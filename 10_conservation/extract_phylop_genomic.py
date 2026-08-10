"""Per-GENE, genomic-coordinate phyloP — full transcript locus incl. introns.

Supersedes the per-isoform CDS-only `conway/*.json` for the genomic
conservation view. Two key wins over the old extractor:

  • Stored ONCE per gene (keyed by genomic coordinate), not once per isoform.
    Isoforms of one gene overlap ~6×, so per-gene dedup turns a 490 Mbp
    per-isoform footprint into an 83 Mbp per-gene union.
  • Covers the FULL locus (first→last feature base, introns included), so the
    genomic view can show intronic conservation — not just CDS.

Values are the RAW phyloP score (not normalised), quantised to int8 to keep
the payload small:

    stored = clamp(round(phyloP * SCALE), -127, 127)   # de-quant: stored/SCALE
    -128 = NO DATA (base had no aligned phyloP)

At SCALE=6 the range is ±21.16 at 0.167 resolution — phyloP for these tracks
sits well inside that, and display only needs ~0.2 precision.

Two tracks per gene:
    z241 : UCSC cactus241way.phyloP.bw   (241-mammal Zoonomia)
    m470 : UCSC hg38.phyloP470way.bw     (470-mammal)

Output (public/mock/conway_genomic/):
    {SYMBOL}.bin      raw bytes: [z241 int8 × span][m470 int8 × span]
    _index.json       { SYMBOL: {chrom, lo, hi, strand, span} }   (lo/hi 1-based incl.)
    _meta.json        scale, track labels, provenance

Frontend: fetch {SYMBOL}.bin as an ArrayBuffer, slice into two Int8Arrays of
length span; value at genomic position g is arr[g - lo] (de-quant /SCALE).

Run (from repo root), bigWigs already downloaded + verified:
    protein2genomic/.venv/bin/python scripts/extract_phylop_genomic.py
"""

from __future__ import annotations

import json
import struct
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyBigWig

warnings.filterwarnings("ignore", message="Mean of empty slice")

REPO = Path(__file__).resolve().parent.parent
TF_DIR = REPO / "public" / "mock" / "tfs"
OUT_DIR = REPO / "public" / "mock" / "conway_genomic"
DATA = REPO.parent / "protein2genomic_data"
TRACKS = {
    "z241": DATA / "phylop" / "cactus241way.phyloP.bw",
    "m470": DATA / "phylop" / "hg38.phyloP470way.bw",
}
SCALE = 6
NODATA = -128


def gene_loci() -> dict[str, dict]:
    """One genomic span per gene: union of every isoform's feature extent.

    Multi-chrom genes (shouldn't happen for these TFs) keep the chrom with the
    widest span. Returns {SYMBOL: {chrom, lo, hi, strand}} (1-based inclusive).
    """
    loci: dict[str, dict] = {}
    for f in sorted(TF_DIR.glob("*.json")):
        d = json.loads(f.read_text())
        sym = d.get("symbol") or f.stem
        by_chrom: dict[str, list[int]] = defaultdict(list)  # chrom -> [lo,hi]
        strand_by_chrom: dict[str, str] = {}
        for s in (d.get("genomicStructures") or []):
            feats = [ft for ft in (s.get("features") or [])
                     if ft.get("start") and ft.get("end")]
            if not feats:
                continue
            ch = s.get("chrom") or "?"
            lo = min(ft["start"] for ft in feats)
            hi = max(ft["end"] for ft in feats)
            ent = by_chrom[ch]
            if not ent:
                ent.extend([lo, hi])
            else:
                ent[0] = min(ent[0], lo)
                ent[1] = max(ent[1], hi)
            strand_by_chrom.setdefault(ch, s.get("strand") or "+")
        if not by_chrom:
            continue
        # widest-span chrom wins
        ch = max(by_chrom, key=lambda c: by_chrom[c][1] - by_chrom[c][0])
        lo, hi = by_chrom[ch]
        loci[sym] = {"chrom": ch, "lo": int(lo), "hi": int(hi),
                     "strand": strand_by_chrom.get(ch, "+")}
    return loci


def quantise(vals: np.ndarray) -> bytes:
    """phyloP float array -> int8 bytes (NaN -> NODATA sentinel)."""
    q = np.where(np.isnan(vals), float(NODATA),
                 np.clip(np.round(vals * SCALE), -127, 127))
    return q.astype(np.int8).tobytes()


def main() -> None:
    for p in TRACKS.values():
        if not p.exists():
            sys.exit(f"missing bigWig: {p}")

    loci = gene_loci()
    print(f"[phylop] {len(loci)} genes; span "
          f"{sum(l['hi'] - l['lo'] + 1 for l in loci.values()):,} bp", flush=True)

    bws = {name: pyBigWig.open(str(path)) for name, path in TRACKS.items()}
    chrom_sets = {name: set(bw.chroms()) for name, bw in bws.items()}

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for old in OUT_DIR.glob("*.bin"):
        old.unlink()

    index: dict[str, dict] = {}
    written = 0
    total_bytes = 0
    for i, (sym, loc) in enumerate(sorted(loci.items())):
        chrom, lo, hi = loc["chrom"], loc["lo"], loc["hi"]
        span = hi - lo + 1
        blob = bytearray()
        ok = False
        for tname, bw in bws.items():
            if chrom not in chrom_sets[tname]:
                blob.extend(bytes([NODATA & 0xFF]) * span)
                continue
            vals = np.array(bw.values(chrom, lo - 1, hi), dtype=np.float64)
            if vals.size != span:  # clamp to chrom edge edge-cases
                vals = np.resize(vals, span)
            ok = True
            blob.extend(quantise(vals))
        if not ok:
            continue
        (OUT_DIR / f"{sym}.bin").write_bytes(bytes(blob))
        index[sym] = {"chrom": chrom, "lo": lo, "hi": hi,
                      "strand": loc["strand"], "span": span}
        written += 1
        total_bytes += len(blob)
        if (i + 1) % 200 == 0:
            print(f"  {i + 1}/{len(loci)} ({total_bytes / 1e6:.0f} MB)", flush=True)

    (OUT_DIR / "_index.json").write_text(json.dumps(index, separators=(",", ":")))
    (OUT_DIR / "_meta.json").write_text(json.dumps({
        "description": "Per-gene genomic-coordinate phyloP (full locus incl. "
                       "introns), raw scores quantised int8 (stored/SCALE).",
        "format": "{SYMBOL}.bin = [z241 int8 × span][m470 int8 × span]; "
                  "value at genomic g = arr[g-lo]/SCALE; -128 = no data.",
        "scale": SCALE,
        "nodata": NODATA,
        "tracks": {
            "z241": {"label": "phyloP 241-way (Zoonomia)",
                     "source": "UCSC cactus241way.phyloP.bw (hg38)"},
            "m470": {"label": "phyloP 470-way (mammals)",
                     "source": "UCSC hg38.phyloP470way.bw"},
        },
        "n_genes": written,
        "total_bytes": total_bytes,
    }, indent=2))

    print(f"[phylop] wrote {written} gene bins, {total_bytes / 1e6:.0f} MB total "
          f"-> {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
