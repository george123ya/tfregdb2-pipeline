"""Download ClinVar's variant_summary.txt.gz and emit a slimmed TSV for the
TFRegDB2 variants build.

variant_summary.txt is NCBI's tab-delimited per-variant view. We keep only:
  - GRCh38 rows (drop GRCh37 duplicates)
  - SNV / indel / dup / del types (skip structural variants and CNVs)
  - Variants that have an explicit chrom/pos/ref/alt (PositionVCF columns)

Output schema (tab-separated):
    chrom, pos, ref, alt,
    clinvar_id, rsid, gene, name (HGVS),
    type, consequence, significance, review_status, origin, pmids
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
import urllib.request
from pathlib import Path

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

# Keep only the variant types we can project to AA position (point mutations
# plus small indels). Structural variants live elsewhere.
KEEP_TYPES = {
    "single nucleotide variant",
    "Deletion",
    "Insertion",
    "Indel",
    "Duplication",
    "Microsatellite",
}


def _download(out_path: Path) -> None:
    if out_path.exists():
        print(f"  already cached: {out_path} ({out_path.stat().st_size / 1e6:.1f} MB)")
        return
    print(f"  downloading {CLINVAR_URL} → {out_path}...")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(CLINVAR_URL, out_path)
    print(f"  downloaded {out_path.stat().st_size / 1e6:.1f} MB")


def _norm_chrom(c: str) -> str:
    """Normalise chromosome to `chr1` / `chrX` / etc. (matches prot2exon)."""
    c = c.strip()
    if not c:
        return ""
    if c.startswith("chr"):
        return c
    return f"chr{c}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path,
                    default=Path("data/variants/source/variant_summary.txt.gz"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/variants/clinvar.tsv"))
    args = ap.parse_args()

    print("[fetch_clinvar] step 1: download")
    _download(args.cache)

    print("[fetch_clinvar] step 2: parse + filter")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_in = n_kept = n_no_pos = n_wrong_assembly = n_wrong_type = 0
    with gzip.open(args.cache, "rt") as fh, args.out.open("w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "chrom", "pos", "ref", "alt",
            "clinvar_id", "rsid", "gene", "name",
            "type", "consequence", "significance", "review_status",
            "origin", "pmids",
        ])
        for row in reader:
            n_in += 1
            if row.get("Assembly") != "GRCh38":
                n_wrong_assembly += 1
                continue
            vtype = row.get("Type", "")
            if vtype not in KEEP_TYPES:
                n_wrong_type += 1
                continue
            pos = row.get("PositionVCF") or ""
            ref = row.get("ReferenceAlleleVCF") or ""
            alt = row.get("AlternateAlleleVCF") or ""
            chrom = _norm_chrom(row.get("Chromosome", ""))
            if not (pos and pos != "-1" and chrom and ref and alt):
                n_no_pos += 1
                continue
            rsid = row.get("RS# (dbSNP)") or ""
            if rsid and rsid != "-1":
                rsid = f"rs{rsid}"
            else:
                rsid = ""
            writer.writerow([
                chrom, pos, ref, alt,
                row.get("VariationID", ""),
                rsid,
                row.get("GeneSymbol", ""),
                row.get("Name", ""),
                vtype,
                # variant_summary.txt doesn't carry consequence directly; we
                # infer from `Name` (HGVS) downstream. Leave blank here.
                "",
                row.get("ClinicalSignificance", ""),
                row.get("ReviewStatus", ""),
                row.get("OriginSimple", ""),
                # citation PMIDs aren't in variant_summary.txt either — they
                # live in var_citations.txt.gz. Add a separate fetch later
                # (skipped for v1 to keep this step single-file).
                "",
            ])
            n_kept += 1
            if n_kept % 500000 == 0:
                print(f"  ... kept {n_kept:,}")

    print()
    print(f"  rows in:              {n_in:,}")
    print(f"  kept:                 {n_kept:,}")
    print(f"  skipped (assembly):   {n_wrong_assembly:,}")
    print(f"  skipped (type):       {n_wrong_type:,}")
    print(f"  skipped (no VCF pos): {n_no_pos:,}")
    print(f"  → {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
