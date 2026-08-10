"""
Read the prot2exon `isoform_structure.tsv` + `domain_mapping_summary.tsv` and
fold a compact `genomicStructure` blob into each per-TF JSON.

Shape per TF (when the canonical ENSP appears in the prot2exon outputs):
    {
      "transcriptId": "ENST00000...",
      "proteinId":   "ENSP00000...",
      "chrom":       "chr17",
      "strand":      "+" | "-",
      "isCanonical": bool,
      "isManeSelect": bool,
      "cdsLength":   int,   # total CDS length in nt
      "features": [
          {"type": "5utr"|"cds"|"intron"|"3utr",
           "start": 1-based-incl,
           "end":   1-based-incl,
           "order": 1..N (translation order),
           "aaStart"?, "aaEnd"?   # CDS only
          }, ...
      ],
      "domains": [
          {"aaStart": 1, "aaEnd": 75,
           "genomicStart": ..., "genomicEnd": ...,
           "type": "AD"|"RD"|"DNA"|"Oligo"|"Other",   # cross-referenced from domainReports
           "source": "TFRegDB1" | "Tycko" | ...},
          ...
      ]
    }

We deliberately dedupe features per transcript (the TSV repeats them for every
domain query). One canonical transcript per TF.
"""

from __future__ import annotations
import csv
import json
import os
import re
from collections import defaultdict
from typing import Optional


PROT2EXON_RESULTS = os.environ.get(
    "PROT2EXON_RESULTS",
    os.path.expanduser("~/Desktop/tfregdb_source/prot2exon/results"),
)
TFS_DIR = os.environ.get(
    "TFS_DIR", os.path.expanduser("~/Desktop/tfregdb2/public/mock/tfs")
)

FEATURE_TYPE_NORMALIZE = {
    "five_prime_UTR": "5utr",
    "CDS": "cds",
    "intron": "intron",
    "three_prime_UTR": "3utr",
}


def _na_int(v: str) -> Optional[int]:
    if v == "NA" or v == "":
        return None
    try:
        return int(v)
    except ValueError:
        return None


def load_domains_summary() -> dict[str, list[dict]]:
    """Return gene_name -> list of domain dicts (aa range + genomic envelope)."""
    by_gene: dict[str, list[dict]] = defaultdict(list)
    path = os.path.join(PROT2EXON_RESULTS, "domain_mapping_summary.tsv")
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("status") != "ok":
                continue
            gene = row.get("gene_name") or ""
            if not gene:
                continue
            by_gene[gene].append({
                "id": row["domain_id"],
                "transcriptId": row.get("transcript_id") or "",
                "aaStart": _na_int(row["aa_start"]),
                "aaEnd": _na_int(row["aa_end"]),
                "genomicStart": _na_int(row["domain_genomic_start"]),
                "genomicEnd": _na_int(row["domain_genomic_end"]),
            })
    return by_gene


def load_canonical_structures() -> dict[str, dict]:
    """Return gene_name -> { transcriptId, proteinId, chrom, strand, features[],
       isCanonical, isManeSelect, cdsLength } for the canonical transcript only.

       isoform_structure.tsv repeats the same features for every domain query,
       AND splits CDS exons by each domain's overlap. We want the un-split
       feature list, so we dedupe by feature_id and merge the split parts back
       into a single (start, end) span per feature.
    """
    by_gene: dict[str, dict] = {}
    # gene -> feature_id -> aggregator
    fid_agg: dict[str, dict[str, dict]] = defaultdict(dict)

    path = os.path.join(PROT2EXON_RESULTS, "isoform_structure.tsv")
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            gene = row.get("gene_name") or ""
            if not gene:
                continue
            canon_flag = row.get("is_ensembl_canonical") == "true"
            mane_flag = row.get("is_mane_select") == "true"
            tid = row["transcript_id"]
            if gene not in by_gene:
                by_gene[gene] = {
                    "transcriptId": tid,
                    "proteinId": row["protein_id"],
                    "chrom": row["chrom"],
                    "strand": row["strand"],
                    "isCanonical": canon_flag,
                    "isManeSelect": mane_flag,
                }
            entry = by_gene[gene]
            # Lock the canonical transcript: if we picked a non-canonical one
            # first and later see a canonical row, switch and reset features.
            if tid != entry["transcriptId"]:
                if (canon_flag or mane_flag) and not entry["isCanonical"]:
                    entry["transcriptId"] = tid
                    entry["proteinId"] = row["protein_id"]
                    entry["chrom"] = row["chrom"]
                    entry["strand"] = row["strand"]
                    entry["isCanonical"] = canon_flag
                    entry["isManeSelect"] = mane_flag
                    fid_agg[gene].clear()
                else:
                    continue

            feat_type = FEATURE_TYPE_NORMALIZE.get(row["feature_type"])
            if not feat_type:
                continue
            fid = row["feature_id"]
            s = _na_int(row["feature_genomic_start"])
            e = _na_int(row["feature_genomic_end"])
            order = _na_int(row.get("feature_order_transcript") or "")
            aa_s = _na_int(row.get("aa_start_encoded") or "") if feat_type == "cds" else None
            aa_e = _na_int(row.get("aa_end_encoded") or "") if feat_type == "cds" else None

            cur = fid_agg[gene].get(fid)
            if cur is None:
                fid_agg[gene][fid] = {
                    "type": feat_type,
                    "start": s,
                    "end": e,
                    "order": order,
                    "aaStart": aa_s,
                    "aaEnd": aa_e,
                }
            else:
                # Union the genomic range across CDS-split parts.
                if s is not None and (cur["start"] is None or s < cur["start"]):
                    cur["start"] = s
                if e is not None and (cur["end"] is None or e > cur["end"]):
                    cur["end"] = e
                if aa_s is not None and (cur["aaStart"] is None or aa_s < cur["aaStart"]):
                    cur["aaStart"] = aa_s
                if aa_e is not None and (cur["aaEnd"] is None or aa_e > cur["aaEnd"]):
                    cur["aaEnd"] = aa_e

    # Materialise feature lists + total CDS length per gene.
    for gene, entry in by_gene.values().__iter__() if False else by_gene.items():
        feats = list(fid_agg[gene].values())
        feats.sort(key=lambda f: f["order"] if f["order"] is not None else 1e9)
        # Drop aaStart/aaEnd None for non-CDS rows to keep the payload compact.
        for f in feats:
            if f["type"] != "cds":
                f.pop("aaStart", None)
                f.pop("aaEnd", None)
        entry["features"] = feats
        entry["cdsLength"] = sum(
            (f["end"] - f["start"] + 1)
            for f in feats
            if f["type"] == "cds" and f["start"] is not None and f["end"] is not None
        )
    return by_gene


# Domain id patterns to extract source + type from "{SYM}_{source...}_Row{N}"
# or "{SYM}_{source}_{type}_Row{N}".
DOMAIN_ID_RE = re.compile(
    r"^[A-Za-z0-9]+_(?P<source>[A-Za-z0-9]+(?:[A-Za-z0-9]+)?)(?:_(?P<type>AD|RD|DNA|Oligo))?_Row\d+$"
)


def classify_domain(domain_id: str, tf_reports: list[dict]) -> tuple[str, Optional[str]]:
    """Return (type, source). Falls back to the matching domainReport in our
    per-TF JSON when the ID alone is ambiguous."""
    m = DOMAIN_ID_RE.match(domain_id)
    source = m.group("source") if m else None
    typ = m.group("type") if (m and m.group("type")) else None
    if typ:
        return typ, source
    # Try to match by source name in tf.domainReports (where we know type).
    if source and tf_reports:
        for r in tf_reports:
            r_src = r.get("source") or ""
            if r_src.lower().startswith(source.lower()):
                return r.get("type") or "Other", source
    return "Other", source


def main():
    print(f"Loading prot2exon outputs from {PROT2EXON_RESULTS}")
    summaries = load_domains_summary()
    structures = load_canonical_structures()
    print(f"  genes with mapped domains: {len(summaries)}")
    print(f"  genes with canonical structure: {len(structures)}")

    patched = 0
    no_struct = 0
    for fn in sorted(os.listdir(TFS_DIR)):
        if not fn.endswith(".json"):
            continue
        p = os.path.join(TFS_DIR, fn)
        with open(p) as f:
            tf = json.load(f)
        sym = tf["symbol"]
        struct = structures.get(sym)
        if not struct:
            no_struct += 1
            continue

        # Attach domain overlays (one per BED row that targeted this TF).
        doms_summary = summaries.get(sym, [])
        doms_for_gene: list[dict] = []
        for d in doms_summary:
            typ, src = classify_domain(d["id"], tf.get("domainReports", []))
            doms_for_gene.append({
                "type": typ,
                "source": src,
                "aaStart": d["aaStart"],
                "aaEnd": d["aaEnd"],
                "genomicStart": d["genomicStart"],
                "genomicEnd": d["genomicEnd"],
            })
        # Deduplicate identical (type, aaStart, aaEnd) overlays so the genomic
        # track doesn't render five copies of the same RD claim from different
        # source files — keep one per (type, range).
        seen_keys: set[tuple] = set()
        deduped: list[dict] = []
        for d in doms_for_gene:
            key = (d["type"], d["aaStart"], d["aaEnd"])
            if key in seen_keys:
                continue
            seen_keys.add(key)
            deduped.append(d)
        # Sort by aa start so overlays render in protein order.
        deduped.sort(key=lambda x: (x["aaStart"] or 0, x["aaEnd"] or 0))

        tf["genomicStructure"] = {**struct, "domains": deduped}
        with open(p, "w") as f:
            json.dump(tf, f, separators=(",", ":"))
        patched += 1

    print(f"patched {patched} TF JSONs with genomicStructure")
    print(f"  TFs without structure (probably non-Lambert / missing ENSP): {no_struct}")


if __name__ == "__main__":
    main()
