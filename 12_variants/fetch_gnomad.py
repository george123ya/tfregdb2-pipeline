"""Download gnomAD v4.1 per-chromosome VCFs and slice them to the canonical
CDS intervals of our TFs.

Strategy:
  1. Build a per-chromosome interval BED from prot2exon's CDS exon table —
     only canonical ENSPs of TFs that have a per-TF JSON. Padded by 1 codon
     each side so splice-region variants right at the boundary still land
     in the slice.
  2. For each chromosome:
       a. Stream the gnomAD VCF + its .tbi index from the public S3 bucket
          (no auth, ~5 GB chr1 + 2 KB index). We keep the full VCF on disk
          only while extracting; immediately delete after the slice.
       b. Use pysam's tabix-based fetch to iterate only records inside our
          intervals — orders-of-magnitude faster than parsing the whole VCF.
       c. Append every kept site to a single global TSV
          (data/variants/gnomad.tsv). Fields: chrom, pos, ref, alt, AF, AC,
          AN, popmax_af. We skip non-SNV multiallelics that have more than
          one alt at the same site — gnomAD emits one VCF row per (locus,
          alt) so we don't have to split allele lists.

Peak transient disk:
  ~5-8 GB per chromosome (autosomes 1-2 are the biggest). Each VCF is
  removed as soon as the slice completes, so the working set never exceeds
  one chromosome.

Final disk:
  data/variants/gnomad.tsv ≈ 50-150 MB depending on how many variants land
  inside TF CDS regions.

Re-run is idempotent: completed-chromosomes are recorded in a marker file
(`.gnomad_progress.json`) so a partial run resumes where it stopped.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam


# Public, no-auth bucket. v4.1 exomes (807k samples) — the v4 genomes VCF is
# even larger; exomes give us the protein-coding coverage we need.
GNOMAD_URL = (
    "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/exomes/"
    "gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz"
)

# Pad each CDS interval by ±9 bp (one codon worth) so splice-region variants
# right at the exon boundary still land in the slice. Helps catch a few
# splice-disrupting calls without inflating the variant count materially.
PAD = 9

# Standard autosomes + sex; skip chrM and unscaffolded contigs.
DEFAULT_CHROMS = [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY"]


def _build_intervals(p2x_tsv: Path, mock_dir: Path) -> dict[str, list[tuple[int, int]]]:
    """One sorted, non-overlapping interval list per chromosome from the
    canonical ENSPs of TFs we already emitted per-TF JSONs for."""
    print(f"[gnomad] reading prot2exon {p2x_tsv}...")
    df = pd.read_csv(p2x_tsv, sep="\t", low_memory=False, usecols=[
        "protein_id", "feature_type", "chrom", "feature_genomic_start",
        "feature_genomic_end", "is_ensembl_canonical",
    ])
    # Canonical-flag rows AND the CDS sub-set.
    df = df[(df["feature_type"] == "CDS") & (df["is_ensembl_canonical"] == True)].copy()  # noqa: E712
    # Only consider TFs we actually emit JSONs for — skips ~10k irrelevant ENSPs.
    seen_symbols = {p.stem.upper() for p in mock_dir.glob("*.json")}
    df["gene_name_upper"] = df["protein_id"]  # placeholder, see below
    # Re-derive gene→ensp via prot2exon's gene_name column
    df2 = pd.read_csv(p2x_tsv, sep="\t", low_memory=False,
                      usecols=["protein_id", "gene_name", "feature_type", "is_ensembl_canonical"])
    df2 = df2[(df2["feature_type"] == "CDS") & (df2["is_ensembl_canonical"] == True)]  # noqa: E712
    canonical_ensps = (
        df2.assign(gene_upper=df2["gene_name"].astype(str).str.upper())
           .loc[lambda d: d["gene_upper"].isin(seen_symbols), "protein_id"]
           .unique()
    )
    df = df[df["protein_id"].isin(canonical_ensps)]
    print(f"  canonical ENSPs for {len(seen_symbols)} TFs: {len(canonical_ensps)}")

    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for _, r in df.iterrows():
        raw[str(r["chrom"])].append((
            max(1, int(r["feature_genomic_start"]) - PAD),
            int(r["feature_genomic_end"]) + PAD,
        ))

    # Merge overlapping intervals on each chrom — keeps the tabix-fetch list
    # small and avoids double-counting variants on shared exons across iso
    # CDS rows.
    merged: dict[str, list[tuple[int, int]]] = {}
    total_intervals = 0
    for chrom, ivs in raw.items():
        ivs.sort()
        out = [ivs[0]]
        for s, e in ivs[1:]:
            if s <= out[-1][1] + 1:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
        total_intervals += len(out)
    print(f"  merged into {total_intervals} intervals across {len(merged)} chroms")
    return merged


def _download(chrom: str, dest: Path) -> tuple[Path, Path]:
    """Fetch the VCF + its .tbi index for one chromosome. Returns (vcf, tbi)
    paths. Idempotent — skips when both files already exist with non-zero
    size, so a re-run after a partial download just resumes."""
    url = GNOMAD_URL.format(chrom=chrom)
    tbi_url = url + ".tbi"
    vcf_path = dest / f"gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz"
    tbi_path = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
    dest.mkdir(parents=True, exist_ok=True)
    if not tbi_path.exists() or tbi_path.stat().st_size == 0:
        print(f"  fetching .tbi for {chrom}...")
        urllib.request.urlretrieve(tbi_url, tbi_path)
    if not vcf_path.exists() or vcf_path.stat().st_size == 0:
        print(f"  fetching VCF for {chrom} (this is the big one)...")
        # urlretrieve buffers in memory for huge files; use wget if present
        # for streaming progress + retry.
        if shutil.which("wget"):
            subprocess.check_call([
                "wget", "-q", "--show-progress", "-O", str(vcf_path), url
            ])
        else:
            urllib.request.urlretrieve(url, vcf_path)
    return vcf_path, tbi_path


def _slice_one_chrom(
    vcf_path: Path,
    intervals: list[tuple[int, int]],
    out_writer: csv.writer,
    chrom: str,
) -> int:
    """Iterate variants in `intervals` for one chromosome via tabix-backed
    pysam fetch, write to `out_writer`. Returns rows kept."""
    vcf = pysam.VariantFile(str(vcf_path))
    # Pre-detect which AF / AC / AN / popmax-style INFO fields actually
    # exist in this VCF's header. pysam raises ValueError on `.info.get(k)`
    # when k isn't declared, not just when it's missing per-record — so we
    # have to check the header up front rather than try/except per record.
    info_keys = set(vcf.header.info.keys())
    af_k = "AF" if "AF" in info_keys else None
    ac_k = "AC" if "AC" in info_keys else None
    an_k = "AN" if "AN" in info_keys else None
    pop_k = next(
        (k for k in ("AF_grpmax", "popmax_AF", "AF_popmax", "AF_grpmax_joint") if k in info_keys),
        None,
    )
    n_kept = 0
    seen_pos: set[tuple[int, str, str]] = set()  # de-dup overlapping intervals
    for s, e in intervals:
        # pysam.fetch uses 0-based half-open coords internally — pass start-1.
        for rec in vcf.fetch(chrom, max(0, s - 1), e):
            if rec.alts is None or len(rec.alts) != 1:
                continue
            ref = rec.ref or ""
            alt = rec.alts[0] or ""
            if not ref or not alt:
                continue
            key = (rec.pos, ref, alt)
            if key in seen_pos:
                continue
            seen_pos.add(key)
            af = rec.info.get(af_k) if af_k else None
            ac = rec.info.get(ac_k) if ac_k else None
            an = rec.info.get(an_k) if an_k else None
            popmax = rec.info.get(pop_k) if pop_k else None
            if isinstance(af, tuple):
                af = af[0] if af else None
            if isinstance(ac, tuple):
                ac = ac[0] if ac else None
            if isinstance(popmax, tuple):
                popmax = popmax[0] if popmax else None
            out_writer.writerow([
                chrom, rec.pos, ref, alt,
                "" if af is None else f"{af:.8g}",
                "" if ac is None else str(ac),
                "" if an is None else str(an),
                "" if popmax is None else f"{popmax:.8g}",
            ])
            n_kept += 1
    vcf.close()
    return n_kept


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prot2exon", type=Path,
                    default=Path(os.environ.get(
                        "P2E_ISOFORM_TSV",
                        os.path.expanduser(
                            "~/Desktop/tfregdb_source/prot2exon/"
                            "results_isoforms/isoform_structure.tsv"
                        ),
                    )))
    ap.add_argument("--mock-dir", type=Path, default=Path("public/mock/tfs"))
    ap.add_argument("--work-dir", type=Path, default=Path("data/variants/gnomad_work"))
    ap.add_argument("--out", type=Path, default=Path("data/variants/gnomad.tsv"))
    ap.add_argument("--chroms", type=str, default=None,
                    help="Comma-sep chrom subset (defaults to chr1-22+X+Y).")
    ap.add_argument("--keep-vcf", action="store_true",
                    help="Don't delete the per-chrom VCF after slicing.")
    args = ap.parse_args()

    intervals = _build_intervals(args.prot2exon, args.mock_dir)
    chroms = args.chroms.split(",") if args.chroms else DEFAULT_CHROMS

    progress_path = args.work_dir / ".gnomad_progress.json"
    args.work_dir.mkdir(parents=True, exist_ok=True)
    progress: dict[str, int] = {}
    if progress_path.exists():
        try:
            progress = json.loads(progress_path.read_text())
        except json.JSONDecodeError:
            progress = {}

    # Resume mode: append to existing TSV; rewrite header only when starting fresh.
    fresh = not args.out.exists() or len(progress) == 0
    mode = "w" if fresh else "a"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_fh = args.out.open(mode, newline="")
    writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
    if fresh:
        writer.writerow(["chrom", "pos", "ref", "alt", "af", "ac", "an", "popmax_af"])

    total_kept = sum(progress.values())
    for chrom in chroms:
        if chrom in progress:
            print(f"[gnomad] {chrom}: already done ({progress[chrom]:,} kept) — skip")
            continue
        if chrom not in intervals:
            print(f"[gnomad] {chrom}: no intervals — skip")
            continue
        print(f"\n[gnomad] {chrom}: download → slice")
        vcf_path, tbi_path = _download(chrom, args.work_dir)
        try:
            n = _slice_one_chrom(vcf_path, intervals[chrom], writer, chrom)
        except Exception as e:
            out_fh.flush()
            print(f"  ERROR on {chrom}: {e}")
            raise
        out_fh.flush()
        progress[chrom] = n
        progress_path.write_text(json.dumps(progress, indent=2))
        total_kept += n
        print(f"  {chrom}: kept {n:,} variants (running total {total_kept:,})")
        if not args.keep_vcf:
            try:
                vcf_path.unlink()
                tbi_path.unlink()
                print(f"  deleted {vcf_path.name}")
            except OSError as e:
                print(f"  (couldn't delete VCF: {e})")

    out_fh.close()
    print()
    print(f"  total variants:   {total_kept:,}")
    print(f"  → {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")
    print(f"  progress saved at {progress_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
