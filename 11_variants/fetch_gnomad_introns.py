"""Fetch gnomAD v4.1 GENOME variants over the FULL gene locus (introns
included) for every TF — the missing piece behind "variants also introns".

Two changes vs the old fetch_gnomad.py (exomes, CDS±9bp):
  • gnomAD GENOMES, not exomes — exomes only capture coding bait regions, so
    they have no intronic variants. Genomes cover the whole locus.
  • Intervals = each gene's FULL transcript span (from the phyloP per-gene
    index), not just CDS — so introns are swept.

Density reality: gnomAD genomes carry ~300 variants/kb, ~25M over all TF
introns — far too many (mostly singletons). We keep only AF ≥ MIN_AF
(default 0.1%), the population-relevant intronic variation, which bounds the
set to ~1.3M (~250 MB merged). Coding variants are unaffected — they come
from the existing ClinVar/COSMIC/exome pipeline and keep full resolution.

Speed: PARALLEL across chromosomes (the user's ask — not one-at-a-time).
htslib does remote range reads over the public S3 bucket via the .tbi index,
so NO multi-GB VCF download — each worker opens its chromosome's remote VCF
once and tabix-fetches only our gene intervals. GIL is released during IO, so
a thread per chromosome genuinely overlaps.

Output: data/variants/gnomad_introns.tsv
    symbol, chrom, pos, ref, alt, af, ac, an, popmax_af

Run from repo root (env with pysam):
    python \
      scripts/variants/fetch_gnomad_introns.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pysam

GENOMES_URL = (
    "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/"
    "gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz"
)
MIN_AF = 0.001  # keep AF ≥ 0.1% — bounds the intronic set to ~1.3M variants


def gene_intervals(index_path: Path) -> dict[str, list[tuple[str, int, int]]]:
    """chrom -> [(symbol, lo, hi)], from the phyloP per-gene index."""
    idx = json.loads(index_path.read_text())
    by_chrom: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    for sym, e in idx.items():
        by_chrom[e["chrom"]].append((sym, int(e["lo"]), int(e["hi"])))
    for lst in by_chrom.values():
        lst.sort(key=lambda t: t[1])
    return by_chrom


def fetch_chrom(chrom: str, genes: list[tuple[str, int, int]],
                min_af: float) -> list[list]:
    """Open the chrom's remote VCF once; tabix-fetch each gene interval."""
    url = GENOMES_URL.format(chrom=chrom)
    rows: list[list] = []
    try:
        vf = pysam.VariantFile(url)
    except Exception as e:  # noqa: BLE001
        print(f"  [{chrom}] open failed: {e}", flush=True)
        return rows
    info = set(vf.header.info.keys())
    af_k = "AF" if "AF" in info else None
    ac_k = "AC" if "AC" in info else None
    an_k = "AN" if "AN" in info else None
    pop_k = next((k for k in ("AF_grpmax", "AF_grpmax_joint", "AF_popmax")
                  if k in info), None)
    seen: set[tuple[int, str, str]] = set()
    kept = 0
    for sym, lo, hi in genes:
        try:
            it = vf.fetch(chrom, max(0, lo - 1), hi)
        except Exception:  # noqa: BLE001
            continue
        for rec in it:
            if rec.alts is None or len(rec.alts) != 1:
                continue
            ref, alt = rec.ref or "", rec.alts[0] or ""
            if not ref or not alt:
                continue
            af = rec.info.get(af_k) if af_k else None
            if isinstance(af, tuple):
                af = af[0] if af else None
            if af is None or af < min_af:
                continue
            key = (rec.pos, ref, alt)
            if key in seen:
                continue
            seen.add(key)
            ac = rec.info.get(ac_k) if ac_k else None
            an = rec.info.get(an_k) if an_k else None
            pop = rec.info.get(pop_k) if pop_k else None
            if isinstance(ac, tuple):
                ac = ac[0] if ac else None
            if isinstance(pop, tuple):
                pop = pop[0] if pop else None
            rows.append([sym, chrom, rec.pos, ref, alt,
                         f"{af:.8g}",
                         "" if ac is None else str(ac),
                         "" if an is None else str(an),
                         "" if pop is None else f"{pop:.8g}"])
            kept += 1
    vf.close()
    print(f"  [{chrom}] {len(genes)} genes → {kept:,} variants (AF≥{min_af})",
          flush=True)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=Path,
                    default=Path("public/mock/conway_genomic/_index.json"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/variants/gnomad_introns.tsv"))
    ap.add_argument("--min-af", type=float, default=MIN_AF)
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    by_chrom = gene_intervals(args.index)
    chroms = sorted(by_chrom, key=lambda c: -len(by_chrom[c]))
    print(f"[gnomad-introns] {sum(len(v) for v in by_chrom.values())} genes "
          f"across {len(chroms)} chroms; {args.workers} parallel workers", flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    all_rows: list[list] = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(fetch_chrom, c, by_chrom[c], args.min_af): c
                for c in chroms}
        for fut in as_completed(futs):
            all_rows.extend(fut.result())

    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["symbol", "chrom", "pos", "ref", "alt", "af", "ac", "an", "popmax_af"])
        w.writerows(all_rows)
    print(f"[gnomad-introns] {len(all_rows):,} variants in {(time.time()-t0)/60:.1f} min "
          f"→ {args.out} ({args.out.stat().st_size/1e6:.0f} MB)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
