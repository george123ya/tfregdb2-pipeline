"""Add INTRONIC ClinVar + COSMIC variants to the per-TF JSONs.

The original variant pipeline projects every variant onto a protein position
via the CDS, so anything outside coding exons (introns, deep splice, UTR) is
dropped — that's why only gnomAD currently has intron data (added separately).
This does the same for the two clinical/somatic sources, straight from their
local raw files, filtered to each gene's full transcript locus.

Sources (already downloaded under data/variants/source/):
  • ClinVar  variant_summary.txt.gz                 (free, NCBI) — splice /
    deep-intronic pathogenic variants are clinically the point here.
  • COSMIC   Cosmic_GenomeScreensMutant_v103 + Cosmic_CompleteTargetedScreens
             (academic licence) — genome-wide so it carries intronic somatic calls.

Loci come from the per-gene phyloP index (full transcript incl. introns). A
variant is matched to a gene by genomic position. Records already present in
the per-TF file (coding, by genomic pos+ref+alt) are skipped, so the net add is
the non-coding set. Added records get pos=0 (genomic-view only), consequence
"intronic", and the same clinvar{}/cosmicId shape the coding records use.

Idempotent: prior `clinvar-intron:` / `cosmic-intron:` records are dropped and
rebuilt each run. COSMIC intronic is capped to recurrent sites (sample count ≥
COSMIC_MIN_SAMPLES) to bound size.

Run from repo root:
    python3 scripts/variants/add_noncoding_introns.py
"""

from __future__ import annotations

import bisect
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
TF_DIR = REPO / "public" / "mock" / "variants"
SRC = REPO / "data" / "variants" / "source"
INDEX = REPO / "public" / "mock" / "conway_genomic" / "_index.json"

CLINVAR = SRC / "variant_summary.txt.gz"
COSMIC_VCFS = [
    SRC / "Cosmic_GenomeScreensMutant_v103_GRCh38.vcf.gz",
    SRC / "Cosmic_CompleteTargetedScreensMutant_v103_GRCh38.vcf.gz",
]
CLINVAR_KEEP = {"single nucleotide variant", "Deletion", "Insertion",
                "Indel", "Duplication", "Microsatellite"}
COSMIC_MIN_SAMPLES = 2  # recurrent intronic somatic only — bounds size

csv.field_size_limit(1 << 24)


def build_loci() -> dict[str, tuple[list[int], list[tuple[int, str]]]]:
    """chrom -> (sorted starts, [(hi, symbol)]) for position→gene lookup."""
    idx = json.loads(INDEX.read_text())
    by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for sym, e in idx.items():
        by_chrom[e["chrom"]].append((int(e["lo"]), int(e["hi"]), sym))
    out: dict[str, tuple[list[int], list[tuple[int, str]]]] = {}
    for chrom, lst in by_chrom.items():
        lst.sort()
        out[chrom] = ([x[0] for x in lst], [(x[1], x[2]) for x in lst])
    return out


def gene_at(loci, chrom: str, pos: int) -> str | None:
    ent = loci.get(chrom)
    if not ent:
        return None
    starts, hisym = ent
    i = bisect.bisect_right(starts, pos) - 1
    # check the candidate and a couple neighbours (loci may nest slightly)
    for j in (i, i - 1, i + 1):
        if 0 <= j < len(starts) and starts[j] <= pos <= hisym[j][0]:
            return hisym[j][1]
    return None


def norm_chrom(c: str) -> str:
    c = c.strip()
    if not c:
        return ""
    return c if c.startswith("chr") else f"chr{c}"


def sig_tag(s: str) -> str:
    s = s.lower()
    if "conflicting" in s:
        return "conflicting"
    if "drug response" in s:
        return "drug_response"
    if "likely pathogenic" in s:
        return "likely_pathogenic"
    if "pathogenic" in s:
        return "pathogenic"
    if "likely benign" in s:
        return "likely_benign"
    if "benign" in s:
        return "benign"
    if "uncertain" in s:
        return "vus"
    return "other"


def collect_clinvar(loci) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = defaultdict(list)
    n = kept = 0
    with gzip.open(CLINVAR, "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            n += 1
            if n % 1_000_000 == 0:
                print(f"  [clinvar] {n:,} rows · {kept:,} kept", flush=True)
            if row.get("Assembly") != "GRCh38" or row.get("Type") not in CLINVAR_KEEP:
                continue
            pos = row.get("PositionVCF") or ""
            ref = row.get("ReferenceAlleleVCF") or ""
            alt = row.get("AlternateAlleleVCF") or ""
            chrom = norm_chrom(row.get("Chromosome", ""))
            if not (pos and pos != "-1" and chrom and ref and alt and ref != alt):
                continue
            sym = gene_at(loci, chrom, int(pos))
            if not sym:
                continue
            rs = row.get("RS# (dbSNP)") or ""
            sig = row.get("ClinicalSignificance", "")
            out[sym].append({
                "gpos": int(pos), "ref": ref, "alt": alt,
                "rsid": f"rs{rs}" if rs and rs != "-1" else None,
                "clinvar": {
                    "id": row.get("VariationID", ""),
                    "name": row.get("Name", ""),
                    "significance": sig,
                    "significanceTag": sig_tag(sig),
                    "reviewStatus": row.get("ReviewStatus", ""),
                    "origin": row.get("OriginSimple", ""),
                },
            })
            kept += 1
    print(f"  [clinvar] done: {n:,} rows, {kept:,} in-locus", flush=True)
    return out


def parse_info(info: str) -> dict[str, str]:
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def collect_cosmic(loci) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = defaultdict(list)
    for vcf in COSMIC_VCFS:
        if not vcf.exists():
            print(f"  [cosmic] missing {vcf.name} — skip", flush=True)
            continue
        n = kept = 0
        with gzip.open(vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                n += 1
                if n % 2_000_000 == 0:
                    print(f"  [cosmic:{vcf.name[:24]}] {n:,} · {kept:,} kept", flush=True)
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                chrom = norm_chrom(f[0])
                try:
                    pos = int(f[1])
                except ValueError:
                    continue
                ref, alt = f[3], f[4]
                if not ref or not alt or ref == alt:
                    continue
                sym = gene_at(loci, chrom, pos)
                if not sym:
                    continue
                info = parse_info(f[7])
                sc = int(info.get("GENOME_SCREEN_SAMPLE_COUNT")
                         or info.get("TARGETED_SCREEN_SAMPLE_COUNT") or "1")
                if sc < COSMIC_MIN_SAMPLES:
                    continue
                out[sym].append({
                    "gpos": pos, "ref": ref, "alt": alt,
                    "cosmicId": f[2] or info.get("LEGACY_ID"),
                    "count": sc,
                })
                kept += 1
        print(f"  [cosmic] {vcf.name}: {n:,} rows, {kept:,} in-locus (≥{COSMIC_MIN_SAMPLES} samples)", flush=True)
    return out


def merge(clinvar: dict, cosmic: dict) -> None:
    syms = set(clinvar) | set(cosmic)
    n_cv = n_co = n_files = 0
    before = after = 0
    for sym in sorted(syms):
        fp = TF_DIR / f"{sym}.json"
        if not fp.exists():
            continue
        before += fp.stat().st_size
        d = json.loads(fp.read_text())
        variants = [v for v in d.get("variants", [])
                    if not str(v.get("id", "")).startswith(("clinvar-intron:", "cosmic-intron:"))]
        seen = {(v.get("genomicPos"), v.get("ref"), v.get("alt")) for v in variants}
        for rec in clinvar.get(sym, []):
            key = (rec["gpos"], rec["ref"], rec["alt"])
            if key in seen:
                continue
            seen.add(key)
            variants.append({
                "id": f"clinvar-intron:{rec['clinvar']['id']}",
                "pos": 0, "type": "Other", "count": 1, "source": "ClinVar",
                "chrom": None,  # backfilled from the file's canonical chrom below
                "genomicPos": rec["gpos"], "ref": rec["ref"], "alt": rec["alt"],
                "aaRef": None, "aaAlt": None, "consequence": "intronic",
                "rsid": rec["rsid"], "clinvar": rec["clinvar"],
                "gnomad": None, "cosmicId": None, "isoforms": None,
            })
            n_cv += 1
        for rec in cosmic.get(sym, []):
            key = (rec["gpos"], rec["ref"], rec["alt"])
            if key in seen:
                continue
            seen.add(key)
            variants.append({
                "id": f"cosmic-intron:{rec['cosmicId']}",
                "pos": 0, "type": "Other", "count": rec["count"], "source": "COSMIC",
                "chrom": None, "genomicPos": rec["gpos"], "ref": rec["ref"], "alt": rec["alt"],
                "aaRef": None, "aaAlt": None, "consequence": "intronic",
                "rsid": None, "clinvar": None, "gnomad": None,
                "cosmicId": rec["cosmicId"], "isoforms": None,
            })
            n_co += 1
        # backfill chrom from the file's canonical chrom (all variants share it)
        chrom = next((v.get("chrom") for v in variants if v.get("chrom")), None)
        for v in variants:
            if v.get("chrom") is None:
                v["chrom"] = chrom
        d["variants"] = variants
        fp.write_text(json.dumps(d, separators=(",", ":")))
        after += fp.stat().st_size
        n_files += 1
    print(f"[merge] +{n_cv:,} ClinVar + {n_co:,} COSMIC intronic into {n_files} files", flush=True)
    print(f"[merge] size {before/1e6:.0f}MB → {after/1e6:.0f}MB (+{(after-before)/1e6:.0f}MB)", flush=True)


def main() -> int:
    if not INDEX.exists():
        sys.exit("missing conway_genomic/_index.json")
    loci = build_loci()
    print(f"[noncoding] {sum(len(v[0]) for v in loci.values())} gene loci", flush=True)
    cv = collect_clinvar(loci) if CLINVAR.exists() else {}
    if not CLINVAR.exists():
        print("  [clinvar] variant_summary.txt.gz missing — skip", flush=True)
    co = collect_cosmic(loci)
    merge(cv, co)
    return 0


if __name__ == "__main__":
    sys.exit(main())
