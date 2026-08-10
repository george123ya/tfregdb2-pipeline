"""Parse locally-provided COSMIC somatic-mutation VCFs and emit a TSV the
variant-DB builder can ingest.

We don't download COSMIC from this script — it requires an academic
licence. The expected workflow:

  1. Log in at https://cancer.sanger.ac.uk/cosmic/download
  2. From the current GRCh38 release, download BOTH:
       • Cosmic_GenomeScreensMutant_v<N>_GRCh38.vcf.bgz   (genome-wide screens)
       • Cosmic_TargetedScreensMutant_v<N>_GRCh38.vcf.bgz (targeted panels)
     Skip `Cosmic_NormalVariants_*` — those are germline-like controls, not
     cancer hotspots, and would dilute the lollipop signal.
  3. Drop both files in `data/variants/source/` (any filename — pass with
     --vcfs if you don't want the default glob).
  4. Run this script.

Slice logic mirrors fetch_gnomad.py: keep only sites that fall inside the
canonical CDS intervals (±9 bp padding) of TFs we already emit per-TF JSONs
for.

We KEEP one row per (chrom, pos, ref, alt, transcript) — i.e. COSMIC's
native per-transcript granularity — so the projector can later read the
iso-specific HGVS.p / consequence directly instead of re-deriving it from
prot2exon. Dedup across the two input files (Genome- + TargetedScreens)
on the same key, summing sample counts.

Output columns:
    chrom, pos, ref, alt, cosmic_id, gene, transcript, hgvs_p, so_term,
    is_canonical, count, consequence
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


PAD = 9
DEFAULT_CHROMS = [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY"]


def _build_intervals(p2x_tsv: Path, mock_dir: Path) -> dict[str, list[tuple[int, int]]]:
    df = pd.read_csv(p2x_tsv, sep="\t", low_memory=False, usecols=[
        "protein_id", "gene_name", "feature_type", "is_ensembl_canonical",
        "chrom", "feature_genomic_start", "feature_genomic_end",
    ])
    df = df[(df["feature_type"] == "CDS") & (df["is_ensembl_canonical"] == True)].copy()  # noqa: E712
    seen_symbols = {p.stem.upper() for p in mock_dir.glob("*.json")}
    df["gene_upper"] = df["gene_name"].astype(str).str.upper()
    df = df[df["gene_upper"].isin(seen_symbols)]
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for _, r in df.iterrows():
        c = str(r["chrom"])
        raw[c].append((
            max(1, int(r["feature_genomic_start"]) - PAD),
            int(r["feature_genomic_end"]) + PAD,
        ))
    merged: dict[str, list[tuple[int, int]]] = {}
    for c, ivs in raw.items():
        ivs.sort()
        out = [ivs[0]]
        for s, e in ivs[1:]:
            if s <= out[-1][1] + 1:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[c] = out
    return merged


def _norm_chrom(c: str) -> str:
    c = c.strip()
    if not c:
        return ""
    return c if c.startswith("chr") else f"chr{c}"


def _in_intervals(chrom: str, pos: int, ivs: dict[str, list[tuple[int, int]]]) -> bool:
    """Linear scan over a chrom's intervals — fine because the list is short
    (TF CDS regions on one chromosome rarely exceed a few hundred intervals
    after merge)."""
    arr = ivs.get(chrom)
    if not arr:
        return False
    for s, e in arr:
        if pos < s:
            return False
        if pos <= e:
            return True
    return False


def _parse_info(info_str: str) -> dict[str, str]:
    """COSMIC VCF info is `KEY=VAL;KEY=VAL;FLAG`."""
    out: dict[str, str] = {}
    for piece in info_str.split(";"):
        if "=" in piece:
            k, v = piece.split("=", 1)
            out[k] = v
        elif piece:
            out[piece] = ""
    return out


def _iter_vcf(vcf_path: Path):
    """Open .vcf, .vcf.gz, or .vcf.bgz transparently."""
    p = str(vcf_path)
    if p.endswith(".gz") or p.endswith(".bgz"):
        return gzip.open(p, "rt")
    return open(p, "r")


def main() -> int:
    # Match both bgzipped (.vcf.bgz) and gzipped (.vcf.gz) — COSMIC v103+
    # ships .vcf.gz inside tar archives, while older releases used .vcf.bgz.
    default_glob = "data/variants/source/Cosmic_*creensMutant*v*GRCh38.vcf.*z"
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcfs", type=str, default=default_glob,
                    help=("Glob OR comma-separated list of COSMIC VCF paths. "
                          "Default picks up both Genome- and Targeted-"
                          "ScreensMutant files dropped under data/variants/source/."))
    ap.add_argument("--prot2exon", type=Path,
                    default=Path("/home/goguxor/Desktop/tfregdb_source/prot2exon/results_isoforms/isoform_structure.tsv"))
    ap.add_argument("--mock-dir", type=Path, default=Path("public/mock/tfs"))
    ap.add_argument("--out", type=Path, default=Path("data/variants/cosmic.tsv"))
    args = ap.parse_args()

    # Resolve --vcfs into a concrete file list. Skip explicit NormalVariants
    # files even if a user passes them — they'd dilute the cancer signal.
    import glob as _glob
    if "," in args.vcfs:
        candidates = [Path(p.strip()) for p in args.vcfs.split(",") if p.strip()]
    else:
        candidates = [Path(p) for p in _glob.glob(args.vcfs)]
    vcfs = [p for p in candidates if p.exists() and "NormalVariants" not in p.name]
    if not vcfs:
        print(
            f"ERROR: no COSMIC VCF found at {args.vcfs}.\n\n"
            "COSMIC requires an academic licence — log in at\n"
            "  https://cancer.sanger.ac.uk/cosmic/download\n"
            "Download Cosmic_GenomeScreensMutant_v<N>_GRCh38.vcf.bgz AND\n"
            "Cosmic_TargetedScreensMutant_v<N>_GRCh38.vcf.bgz from the\n"
            "current GRCh38 release. Skip Cosmic_NormalVariants_* — those\n"
            "are germline-like controls, not cancer hotspots. Drop both\n"
            "files in data/variants/source/ and rerun.",
            file=sys.stderr,
        )
        return 1
    print(f"[cosmic] input VCFs:")
    for v in vcfs:
        print(f"    {v}")

    intervals = _build_intervals(args.prot2exon, args.mock_dir)
    print(f"[cosmic] intervals across {len(intervals)} chroms")

    # Dedup across files by (chrom, pos, ref, alt, transcript) so each
    # COSMIC per-transcript row survives but a mutation observed in both
    # Genome- and Targeted-screen files contributes once with combined
    # sample counts. TARGETED_SCREEN_SAMPLE_COUNT / GENOMIC_SCREEN_SAMPLE_COUNT
    # are kept separate in CNT-summing fallback.
    #
    # Versions on the TRANSCRIPT field (.2 suffix on ENST00000641515.2) are
    # stripped so the projector can match against unversioned ENSTs from
    # the index without needing to know which Ensembl release each side
    # references.
    def _strip_ver(tid: str) -> str:
        return tid.split(".", 1)[0] if tid else tid

    aggregate: dict[tuple[str, int, str, str, str], dict] = {}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_in = n_kept = 0
    for vcf_path in vcfs:
        print(f"[cosmic] scanning {vcf_path}")
        with _iter_vcf(vcf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                n_in += 1
                if n_in % 500000 == 0:
                    print(f"  scanned {n_in:,} ({n_kept:,} kept, {len(aggregate):,} unique)")
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                chrom = _norm_chrom(parts[0])
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue
                if not _in_intervals(chrom, pos, intervals):
                    continue
                cosmic_id = parts[2]
                ref, alt = parts[3], parts[4]
                info = _parse_info(parts[7])
                gene = info.get("GENE", "")
                transcript = _strip_ver(info.get("TRANSCRIPT", ""))
                hgvs_p = info.get("HGVSP", "") or info.get("AA", "")
                so_term = info.get("SO_TERM", "")
                # IS_CANONICAL is "y"/"n" per row. Store the boolean so
                # the projector can prefer canonical rows when iso-aware
                # lookup falls through.
                is_canonical = (info.get("IS_CANONICAL", "n").lower() == "y")
                # COSMIC v103 uses TARGETED_SCREEN_SAMPLE_COUNT for the
                # CompleteTargetedScreensMutant file and a separate
                # GENOMIC_SCREEN_SAMPLE_COUNT for the GenomeScreensMutant.
                # Older releases used CNT or SCORE. We accept all four so
                # this script handles past + future COSMIC schema drift.
                try:
                    cnt = int(
                        info.get("TARGETED_SCREEN_SAMPLE_COUNT")
                        or info.get("GENOMIC_SCREEN_SAMPLE_COUNT")
                        or info.get("CNT")
                        or info.get("SCORE")
                        or "1"
                    )
                except ValueError:
                    cnt = 1
                cons = info.get("AA", "") or info.get("CDS", "")
                key = (chrom, pos, ref, alt, transcript)
                if key in aggregate:
                    aggregate[key]["count"] += cnt
                else:
                    aggregate[key] = {
                        "cosmic_id": cosmic_id, "gene": gene,
                        "transcript": transcript, "hgvs_p": hgvs_p,
                        "so_term": so_term, "is_canonical": is_canonical,
                        "count": cnt, "consequence": cons,
                    }
                n_kept += 1

    with args.out.open("w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "pos", "ref", "alt",
                         "cosmic_id", "gene", "transcript", "hgvs_p",
                         "so_term", "is_canonical", "count", "consequence"])
        for (chrom, pos, ref, alt, transcript), v in sorted(aggregate.items()):
            writer.writerow([chrom, pos, ref, alt,
                             v["cosmic_id"], v["gene"], v["transcript"],
                             v["hgvs_p"], v["so_term"],
                             "1" if v["is_canonical"] else "0",
                             v["count"], v["consequence"]])

    print()
    print(f"  rows scanned:        {n_in:,}")
    print(f"  rows kept (in CDS):  {n_kept:,}")
    print(f"  unique variants:     {len(aggregate):,}")
    print(f"  → {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
