"""Extract per-NUCLEOTIDE conservation (phyloP) for every TF isoform.

This produces the real values behind the "ConSway" tracks on the TF page,
replacing the placeholder `mockConswayScore`. Two hg38 phyloP tracks are
scored ("serían dos scores"):

  - zoonomia241 : UCSC cactus241way.phyloP.bw  (241-mammal Zoonomia alignment)
  - mammal470   : UCSC hg38.phyloP470way.bw    (470-mammal alignment)

Resolution is per *nucleotide* (3 values per codon) so the viewer can show the
1st/2nd-vs-3rd codon-position (wobble) signal, not just a per-residue mean.
Every isoform that resolves to an Ensembl transcript is scored in its OWN
frame, so selecting an isoform on the page reframes the track to that
transcript's conservation.

Pipeline per isoform transcript:

  1. prot2exon maps the transcript (ENST) to its CDS exons, giving for each
     exon: chrom, strand, genomic span, and the 1-based CDS-nt span in
     *translation* order (`cds_nt_start/end`).
  2. pyBigWig reads per-base phyloP over each exon's genomic span.
  3. Bases are placed into a CDS-ordered array via `cds_nt_*` (strand-safe:
     minus strand maps genomic-ascending to cds-descending). The array is the
     coding sequence's per-base phyloP, in translation order.
  4. Clamp to [0, cap]/cap per track (cap = p99 of positive per-base scores
     across all isoforms) -> [0, 1]; negatives (acceleration) -> 0; bases with
     no aligned data -> null.

Output: one file per isoform, keyed by isoform NAME (matches the page's
`requestIsoformName$` selection signal):
    public/mock/conway/{isoformName}.json
      { "name": str, "uniprot": str|null, "length": int (aa),
        "z241": [float|null, ...],   # len == 3 * aa length (codon order)
        "m470": [float|null, ...] }
plus `_meta.json` with caps + provenance.

The page's ConswayPanel fetches {canonicalName}.json on mount and the
selected isoform's file on change, falling back to a per-residue mock when a
file is missing (isoforms with no coding transcript).

Run (from repo root), after the bigWigs have downloaded + verified:
    protein2genomic/.venv/bin/python scripts/extract_conway_phylop.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pyBigWig

warnings.filterwarnings("ignore", message="Mean of empty slice")

REPO = Path(__file__).resolve().parent.parent
TF_DIR = REPO / "public" / "mock" / "tfs"
OUT_DIR = REPO / "public" / "mock" / "conway"

DATA = REPO.parent / "protein2genomic_data"
P2E_BIN = REPO.parent / "protein2genomic" / "build" / "prot2exon"
P2E_IDX = DATA / "human.idx"           # GENCODE v49, chr-prefixed contigs
TRACKS = {
    "z241": DATA / "phylop" / "cactus241way.phyloP.bw",
    "m470": DATA / "phylop" / "hg38.phyloP470way.bw",
}
CAP_PERCENTILE = 99.0


def load_isoforms() -> list[dict]:
    """Every (isoform -> Ensembl transcript) pair across all TFs.

    ENSTs are sourced primarily from ``genomicStructures[].transcriptId`` —
    that covers every genomically-mapped isoform the viewer can show (~7.4k),
    not just the subset that also appears in the ``transcripts`` array (~3.4k,
    which is what this used to read, and why ~4k isoforms had no phyloP). Falls
    back to the transcripts array for anything not in genomicStructures.
    Returns [{name, enst, length(aa), uniprot}]; isoforms with no resolvable
    ENST are dropped.
    """
    out = []
    for f in sorted(TF_DIR.glob("*.json")):
        d = json.loads(f.read_text())
        iso_len = {
            iso["name"]: int(iso["length"])
            for iso in (d.get("isoforms") or [])
            if iso.get("name") and iso.get("length")
        }
        iso_uni = {
            iso["name"]: iso.get("uniprotIso")
            for iso in (d.get("isoforms") or [])
            if iso.get("name")
        }
        enst_by_name = {
            t["name"]: t["ensemblId"]
            for t in (d.get("transcripts") or [])
            if t.get("name") and t.get("ensemblId")
        }
        seen: set[str] = set()
        # 1) genomic structures — the full set of isoforms the page renders.
        for s in (d.get("genomicStructures") or []):
            name = s.get("isoformName")
            enst = s.get("transcriptId")
            if not name or name in seen:
                continue
            if not enst or not str(enst).startswith("ENST"):
                continue
            length = iso_len.get(name) or s.get("isoformLength")
            if not length:
                continue
            seen.add(name)
            out.append({"name": name, "enst": enst,
                        "length": int(length), "uniprot": iso_uni.get(name)})
        # 2) transcripts-array isoforms not covered by a genomic structure.
        for iso in (d.get("isoforms") or []):
            name = iso.get("name")
            if not name or name in seen or not iso.get("length"):
                continue
            enst = enst_by_name.get(name)
            if not enst:
                continue
            seen.add(name)
            out.append({"name": name, "enst": enst,
                        "length": int(iso["length"]), "uniprot": iso.get("uniprotIso")})
    return out


def run_prot2exon(isos: list[dict], workdir: Path) -> Path:
    """Map every transcript to its CDS structure in one batch call.

    input_id (BED col 4) is the isoform name so outputs key back directly.
    Duplicate (name, enst) rows are de-duplicated. No aa range -> whole
    transcript (no-domain) mode.
    """
    bed = workdir / "queries.bed"
    seen = set()
    with bed.open("w") as fh:
        for t in isos:
            key = (t["name"], t["enst"])
            if key in seen:
                continue
            seen.add(key)
            fh.write(f"{t['enst']}\t0\t0\t{t['name']}\n")
    out = workdir / "p2e"
    subprocess.run(
        [str(P2E_BIN), "--index", str(P2E_IDX), "--bed", str(bed),
         "--out-dir", str(out), "--output", "isoform"],
        check=True, capture_output=True, text=True,
    )
    return out / "isoform_structure.tsv"


def parse_cds(tsv: Path) -> dict[str, dict]:
    """Group CDS exons by input_id (isoform name).

    -> {name: {chrom, strand, blocks:[(g0,g1,c0,c1)]}}, all 1-based inclusive.
    """
    cds: dict[str, dict] = {}
    with tsv.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if r[col["feature_type"]] != "CDS":
                continue
            name = r[col["input_id"]]
            e = cds.setdefault(name, {
                "chrom": r[col["chrom"]],
                "strand": r[col["strand"]],
                "blocks": [],
            })
            e["blocks"].append((
                int(r[col["feature_genomic_start"]]),
                int(r[col["feature_genomic_end"]]),
                int(r[col["cds_nt_start"]]),
                int(r[col["cds_nt_end"]]),
            ))
    return cds


def per_nt_raw(entry: dict, bw: "pyBigWig.pyBigWig") -> np.ndarray:
    """Per-base phyloP over the CDS, in translation (codon) order. NaN = no data."""
    chrom = entry["chrom"]
    if chrom not in bw.chroms():
        return np.array([])
    minus = entry["strand"] == "-"
    total_nt = max(b[3] for b in entry["blocks"])
    cds_vals = np.full(total_nt + 1, np.nan, dtype=np.float64)  # 1-indexed
    for g0, g1, c0, c1 in entry["blocks"]:
        vals = bw.values(chrom, g0 - 1, g1)  # genomic-ascending, len g1-g0+1
        for k, v in enumerate(vals):
            cpos = (c1 - k) if minus else (c0 + k)
            cds_vals[cpos] = v
    return cds_vals[1:]  # drop the 1-indexing pad; length == total_nt


def fit_to_length(raw: np.ndarray, aa_len: int) -> np.ndarray:
    """Trim/pad the per-nt array to exactly 3*aa_len (drops stop codon; pads short)."""
    target = 3 * aa_len
    if raw.shape[0] == target:
        return raw
    if raw.shape[0] > target:
        return raw[:target]
    return np.concatenate([raw, np.full(target - raw.shape[0], np.nan)])


def to_track(raw: np.ndarray, cap: float) -> list:
    out = []
    for v in raw:
        out.append(None if np.isnan(v)
                   else round(float(min(max(v, 0.0), cap)) / cap, 3))
    return out


def main() -> None:
    for p in (P2E_BIN, P2E_IDX, *TRACKS.values()):
        if not p.exists():
            sys.exit(f"missing required input: {p}")

    isos = load_isoforms()
    print(f"[conway] {len(isos)} isoform transcripts; mapping CDS ...", flush=True)
    with tempfile.TemporaryDirectory() as td:
        cds = parse_cds(run_prot2exon(isos, Path(td)))
    print(f"[conway] CDS structures for {len(cds)} transcripts", flush=True)

    bws = {name: pyBigWig.open(str(path)) for name, path in TRACKS.items()}

    # Pass 1: raw per-nt per isoform/track; collect positives for caps.
    raw_by_iso: dict[str, dict] = {}
    pools: dict[str, list[np.ndarray]] = {name: [] for name in TRACKS}
    skipped = []
    for t in isos:
        name = t["name"]
        if name in raw_by_iso:
            continue
        entry = cds.get(name)
        if not entry:
            skipped.append(name)
            continue
        rec = {"uniprot": t["uniprot"], "length": t["length"], "raw": {}}
        for tname, bw in bws.items():
            raw = per_nt_raw(entry, bw)
            rec["raw"][tname] = raw
            if raw.size:
                pos = raw[~np.isnan(raw) & (raw > 0)]
                if pos.size:
                    pools[tname].append(pos)
        raw_by_iso[name] = rec

    caps = {n: float(np.percentile(np.concatenate(pools[n]), CAP_PERCENTILE))
            for n in TRACKS if pools[n]}
    print(f"[conway] caps (p{CAP_PERCENTILE:g}): "
          + ", ".join(f"{k}={v:.3f}" for k, v in caps.items()), flush=True)

    # Pass 2: normalise, fit to 3*aa, write one file per isoform.
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for old in OUT_DIR.glob("*.json"):       # clear stale (per-uniprot) files
        old.unlink()
    written = 0
    for name, rec in raw_by_iso.items():
        aa_len = rec["length"]
        obj = {"name": name, "uniprot": rec["uniprot"], "length": aa_len}
        ok = False
        for tname in TRACKS:
            raw = rec["raw"][tname]
            if not raw.size:
                obj[tname] = [None] * (3 * aa_len)
                continue
            ok = True
            obj[tname] = to_track(fit_to_length(raw, aa_len), caps[tname])
        if ok:
            (OUT_DIR / f"{name}.json").write_text(json.dumps(obj, separators=(",", ":")))
            written += 1

    (OUT_DIR / "_meta.json").write_text(json.dumps({
        "description": "Per-nucleotide phyloP conservation per TF isoform "
                       "transcript (3 values/codon), normalised to [0,1].",
        "tracks": {
            "z241": {"label": "phyloP 241-way (Zoonomia)",
                     "source": "UCSC cactus241way.phyloP.bw (hg38)"},
            "m470": {"label": "phyloP 470-way (mammals)",
                     "source": "UCSC hg38.phyloP470way.bw"},
        },
        "normalization": f"clamp(phyloP,0,cap)/cap; cap = p{CAP_PERCENTILE:g} of "
                         "positive per-base scores across all isoforms",
        "caps": caps,
        "resolution": "per nucleotide; array length == 3 * aa length, codon order",
        "mapping": "prot2exon (GENCODE v49) ENST -> CDS bases",
        "n_written": written,
        "n_skipped_no_cds": len(skipped),
    }, indent=2))

    print(f"[conway] wrote {written} isoform files to {OUT_DIR}", flush=True)
    if skipped:
        print(f"[conway] {len(skipped)} transcripts had no CDS mapping "
              f"(non-coding?): {', '.join(skipped[:8])}"
              f"{' ...' if len(skipped) > 8 else ''}")


if __name__ == "__main__":
    main()
