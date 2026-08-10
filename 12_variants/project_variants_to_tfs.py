"""For each TF, slice the variants SQLite by the canonical ENSP's CDS exons,
project every variant to AA position, and emit a per-TF JSON.

Pipeline:
  1. Load prot2exon's isoform_structure.tsv → per-ENSP CDS map (chrom, strand,
     and a list of (g_start, g_end, aa_start, aa_end) chunks). Apply the
     boundary-codon partitioning (each AA owned by exactly one CDS) so
     variants on the codon junction don't end up double-counted.
  2. Resolve canonical ENSP per gene (Ensembl canonical → APPRIS PRINCIPAL:1
     → longest), matching build_tf_records's chain.
  3. For each TF with a canonical ENSP + a per-TF JSON in public/mock/tfs:
       a. Query SQLite for variants in any of the canonical's CDS genomic
          ranges (chrom + pos BETWEEN g_start AND g_end).
       b. Project each variant: find which CDS exon owns `pos`, compute AA
          position via the same arithmetic as map_canonical_to_isoforms.py.
       c. For every alt isoform of the same gene, attempt the same projection
          onto the iso's CDS map. Variants only "appear" on the iso if its
          CDS includes the codon.
       d. Light HGVS-based consequence inference: pull `p.<X><pos><Y>` from
          the ClinVar `Name` column to label missense vs nonsense vs silent.
  4. Write to public/mock/variants/{SYMBOL}.json AND inject the same list
     into the TF's existing per-TF JSON under `mutations` so the UI loader
     picks it up without a second fetch.

Output JSON schema per TF:
{
  "symbol": "TP53",
  "canonicalEnsp": "ENSP00000269305",
  "variants": [
    {
      "chrom": "chr17", "pos": 7676154, "ref": "G", "alt": "A",
      "aaPos": 175, "aaRef": "R", "aaAlt": "H", "consequence": "missense",
      "rsid": "rs28934578",
      "clinvar": {"id": "12347", "name": "NM_000546.6(TP53):c.524G>A (p.Arg175His)",
                  "significance": "Pathogenic", "reviewStatus": "...",
                  "origin": "germline"},
      "sources": ["clinvar"],
      "isoforms": {"TP53-204": 175, "TP53-205": 175}
    }
  ]
}
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sqlite3
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


# ---- AA3 → AA1 (for HGVS parsing) ---------------------------------------
AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*", "X": "*",
}
# Three-letter HGVS used by ClinVar (`p.Arg175His`, `p.Arg175*`, `p.Arg175del`).
HGVS_P3_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|=|del|fs.*|dup.*)?")
# Single-letter HGVS used by COSMIC (`p.R175H`, `p.L185=`, `p.R209*`).
HGVS_P1_RE = re.compile(r"p\.([A-Z\*])(\d+)([A-Z\*=]|del|fs.*|dup.*)?")


@dataclass
class CdsChunk:
    aa_start: int
    aa_end: int
    g_start: int
    g_end: int


@dataclass
class EnspMap:
    strand: str
    chrom: str
    cds: list[CdsChunk]


def _parse_prot2exon(p2x_tsv: Path) -> dict[str, EnspMap]:
    print(f"[1/4] reading prot2exon {p2x_tsv}...")
    df = pd.read_csv(p2x_tsv, sep="\t", low_memory=False)
    df = df[df["feature_type"] == "CDS"].copy()
    df["aa_start_encoded"] = pd.to_numeric(df["aa_start_encoded"], errors="coerce")
    df["aa_end_encoded"] = pd.to_numeric(df["aa_end_encoded"], errors="coerce")
    df = df.dropna(subset=["aa_start_encoded", "aa_end_encoded"])

    out: dict[str, EnspMap] = {}
    for ensp, sub in df.groupby("protein_id"):
        first_input = sub["input_id"].iloc[0]
        block = sub[sub["input_id"] == first_input].sort_values("feature_order_transcript")
        if block.empty:
            continue
        strand = str(block["strand"].iloc[0])
        chrom = str(block["chrom"].iloc[0])
        # NB: do NOT partition boundary AAs here, unlike the iso-mapping
        # script. prot2exon's `aa_start_encoded` is the first AA *completely
        # encoded* in this CDS exon; the previous exon owns the boundary
        # AA's split codon. Bumping aa_start forward (as the iso-mapping
        # path does to fix chip-rendering bridges) makes the variant→AA
        # formula off-by-one for variants landing at the first nucleotide
        # of an exon. Validated against ClinVar HGVS for AR/RUNX1/...:
        # un-partitioned matches the AA stated in HGVS.
        cds = [
            CdsChunk(
                aa_start=int(r["aa_start_encoded"]),
                aa_end=int(r["aa_end_encoded"]),
                g_start=int(r["feature_genomic_start"]),
                g_end=int(r["feature_genomic_end"]),
            )
            for _, r in block.iterrows()
        ]
        if cds:
            out[ensp] = EnspMap(strand=strand, chrom=chrom, cds=cds)
    print(f"  ENSPs with CDS map: {len(out)}")
    return out


def _resolve_canonical(ids_csv: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Match build_tf_records's canonical picker. Returns (gene→canonical ENSP,
    gene→list of alt ENSPs)."""
    print(f"[2/4] reading TF_completeIDs {ids_csv}...")
    ids = pd.read_csv(ids_csv)
    ids["Translation_Clean"] = ids["Translation ID"].astype(str).str.split(".").str[0]
    ids["Gene_Upper"] = ids["Gene name (HGNC)"].astype(str).str.strip().str.upper()
    ids = ids[ids["Translation_Clean"].str.startswith("ENSP", na=False)].copy()

    canon_ensp: dict[str, str] = {}
    for gene, grp in ids.groupby("Gene_Upper"):
        if not gene or gene == "NAN":
            continue
        pick = None
        ens_canon = grp[grp.get("Ensembl canonical", "").astype(str) == "Canonical"]
        if len(ens_canon):
            pick = ens_canon.iloc[0]
        if pick is None:
            p1 = grp[grp["APPRIS Annotation"].astype(str).str.contains("PRINCIPAL:1", na=False)]
            if len(p1):
                pick = p1.iloc[0]
        if pick is None:
            pany = grp[grp["APPRIS Annotation"].astype(str).str.contains("PRINCIPAL", na=False)]
            if len(pany):
                pick = pany.iloc[0]
        if pick is None and "Protein length" in grp.columns:
            with_len = grp.dropna(subset=["Protein length"])
            if len(with_len):
                pick = with_len.sort_values("Protein length", ascending=False).iloc[0]
        if pick is None:
            pick = grp.iloc[0]
        canon_ensp[gene] = str(pick["Translation_Clean"])

    alt_ensps: dict[str, list[str]] = defaultdict(list)
    for gene, grp in ids.groupby("Gene_Upper"):
        canon = canon_ensp.get(gene)
        for e in grp["Translation_Clean"].dropna().unique():
            if e != canon:
                alt_ensps[gene].append(str(e))
    return canon_ensp, alt_ensps


def _project_genomic_to_aa(
    pos: int, ensp_map: EnspMap, ref_len: int = 1,
) -> int | None:
    """AA whose codon contains the variant, walking CDS in transcript order.

    `pos` is the VCF anchor (leftmost genomic position). `ref_len` is the
    length of the REF allele — for SNVs it's 1; for deletions it spans the
    deleted bases. We scan every position in `[pos, pos + ref_len - 1]`
    looking for a CDS-overlapping base.

    Why scan a range? VCF anchors deletions at the leftmost (lowest)
    genomic position. On minus-strand genes the leftmost position can sit
    just past the stop codon in the 3'UTR even when the deletion *itself*
    overlaps the stop. Looking only at `pos` misses these stop-loss
    variants (FOXP3 c.1293_1294del was a verified case in our audit).

    Per-exon `(aa_start+offset)` arithmetic gets the answer wrong for
    codons that span an exon boundary: the codon "belongs" to the cDNA
    position walked across the splice junction, not to either exon's
    local AA bounds. cDNA-walk handles split codons correctly.
    """
    plus = ensp_map.strand == "+"
    end = pos + max(1, ref_len) - 1
    for scan in range(pos, end + 1):
        cdna = 0
        for c in ensp_map.cds:
            cds_len = c.g_end - c.g_start + 1
            if scan < c.g_start or scan > c.g_end:
                cdna += cds_len
                continue
            offset = (scan - c.g_start) if plus else (c.g_end - scan)
            cdna_pos = cdna + offset + 1  # 1-based cDNA
            aa = (cdna_pos + 2) // 3      # AA 1 = cDNA 1..3, AA 2 = 4..6, …
            return aa
    return None


def _parse_hgvs_consequence(hgvs_name: str) -> tuple[str, str | None, str | None, int | None]:
    """Pull (consequence, aa_ref, aa_alt, aa_pos) from a ClinVar `Name` column
    OR a COSMIC `AA=` INFO field. ClinVar uses three-letter HGVS (`p.Arg175His`)
    while COSMIC uses one-letter (`p.R175H`, `p.L185=`, `p.R209*`). We try
    three-letter first because ClinVar names sometimes contain bare letters
    that look like one-letter codes (e.g. transcript IDs)."""
    if not hgvs_name:
        return "", None, None, None

    # 3-letter (ClinVar) path
    m3 = HGVS_P3_RE.search(hgvs_name)
    ref1: str = ""
    alt: str = ""
    aa_pos: int | None = None
    tail: str = ""  # the slice past the matched alt — exposes fs / Ter / ? suffixes
    if m3:
        ref1 = AA3_TO_AA1.get(m3.group(1).title(), "")
        aa_pos = int(m3.group(2)) if m3.group(2).isdigit() else None
        alt = (m3.group(3) or "").strip()
        tail = hgvs_name[m3.end():m3.end() + 12]
        # Translate 3-letter alt → 1-letter where possible.
        if alt and alt not in ("*", "=", "del") and not alt.startswith(("fs", "dup", "del")):
            alt = AA3_TO_AA1.get(alt.title(), alt)
    else:
        # 1-letter (COSMIC) path
        m1 = HGVS_P1_RE.search(hgvs_name)
        if not m1:
            return "", None, None, None
        ref1 = m1.group(1)
        aa_pos = int(m1.group(2)) if m1.group(2).isdigit() else None
        alt = (m1.group(3) or "").strip()
        tail = hgvs_name[m1.end():m1.end() + 12]

    # "?" anywhere in the substitution slice means the protein effect is
    # unknown / unpredictable (common for start-codon-affecting variants
    # like p.Met1? — the start ATG is destroyed, so we can't say what
    # gets translated). Returning "" makes the UI paint these "Other"
    # rather than "synonymous", which the empty-alt branch below would
    # otherwise do.
    if alt == "?" or tail.startswith("?") or hgvs_name.endswith("p.?"):
        return "", ref1 or None, None, aa_pos
    # Frameshift / nonsense suffix sitting *past* a 3-letter match: the
    # P3_RE regex greedily captures the next 3-letter AA after the
    # position, so `p.Met40AsnfsTer3` becomes group(3)=Asn and the `fs`
    # never reaches the dispatch below. Catch it via the tail.
    if tail.startswith("fs") or tail.startswith("Ter"):
        return "frameshift", ref1, None, aa_pos
    if tail.startswith("delins"):
        return "deletion", ref1, None, aa_pos

    if alt == "=":
        return "synonymous", ref1, ref1, aa_pos
    if alt == "":
        # Empty group(3) on a positioned HGVS = unparseable substitution
        # (e.g. p.Met1, p.Arg40 with no effect captured). Was returning
        # "synonymous" incorrectly; now flagged as Other.
        return "", ref1 or None, None, aa_pos
    if alt in ("*", "X"):
        return "nonsense", ref1, "*", aa_pos
    if alt.startswith("del"):
        return "deletion", ref1, None, aa_pos
    if alt.startswith("dup"):
        return "duplication", ref1, None, aa_pos
    if alt.startswith("fs"):
        return "frameshift", ref1, None, aa_pos
    # By now alt should be a single letter for a missense.
    if alt and ref1 and alt != ref1 and len(alt) == 1 and alt.isalpha():
        return "missense", ref1, alt, aa_pos
    # Same-letter (e.g. M==M) = synonymous, not an unknown effect.
    if alt and ref1 and alt == ref1:
        return "synonymous", ref1, ref1, aa_pos
    return "", ref1 or None, alt or None, aa_pos


def _significance_label(sig: str) -> str:
    """Bucket ClinVar's free-form ClinicalSignificance into a short tag."""
    if not sig:
        return "unknown"
    s = sig.lower()
    if "pathogenic" in s and "likely" not in s and "conflict" not in s:
        return "pathogenic"
    if "likely pathogenic" in s:
        return "likely_pathogenic"
    if "likely benign" in s:
        return "likely_benign"
    if "benign" in s and "likely" not in s:
        return "benign"
    if "conflict" in s:
        return "conflicting"
    if "uncertain" in s or "vus" in s:
        return "vus"
    if "drug response" in s:
        return "drug_response"
    if "risk" in s:
        return "risk_factor"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path, default=Path("data/variants/variants.sqlite"))
    ap.add_argument("--ids-csv", type=Path,
                    default=Path(os.environ.get(
                        "COMPLETE_IDS_CSV",
                        os.path.expanduser(
                            "~/Desktop/tfregdb/data_pipeline/blast_table/"
                            "TF_completeIDs_with_ensembl_canonical.csv"
                        ),
                    )))
    ap.add_argument("--prot2exon", type=Path,
                    default=Path(os.environ.get(
                        "P2E_ISOFORM_TSV",
                        os.path.expanduser(
                            "~/Desktop/tfregdb_source/prot2exon/"
                            "results_isoforms/isoform_structure.tsv"
                        ),
                    )))
    ap.add_argument("--mock-dir", type=Path,
                    default=Path("public/mock/tfs"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("public/mock/variants"))
    ap.add_argument("--limit", type=int, default=None,
                    help="Cap number of TFs (for testing).")
    ap.add_argument("--symbols", type=str, default=None,
                    help="Comma-sep list of TF symbols to process.")
    args = ap.parse_args()

    p2x = _parse_prot2exon(args.prot2exon)
    canon_ensp, alt_ensps = _resolve_canonical(args.ids_csv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    cur = con.cursor()

    # Iterate TFs that have a per-TF JSON in mock/tfs/ (i.e. that build_tf_records
    # produced) so we don't try to write variants for symbols the rest of the
    # build skipped.
    tf_files = sorted(args.mock_dir.glob("*.json"))
    if args.symbols:
        wanted = {s.strip().upper() for s in args.symbols.split(",")}
        tf_files = [p for p in tf_files if p.stem.upper() in wanted]
    if args.limit:
        tf_files = tf_files[: args.limit]
    print(f"[3/4] projecting variants over {len(tf_files)} TFs")

    n_emitted = n_with_variants = n_total_variants = 0
    counts: dict[str, int] = {}
    for fp in tf_files:
        symbol = fp.stem
        gene_u = symbol.upper()
        c_ensp = canon_ensp.get(gene_u)
        if not c_ensp:
            continue
        canon_map = p2x.get(c_ensp)
        if canon_map is None:
            continue

        # SQL: variants in the union of canonical-CDS + every alt-iso-CDS
        # genomic range on the same chromosome. Each interval is padded by
        # INDEL_PAD bp on both sides so VCF-anchored deletions whose
        # leftmost position lands just outside the CDS (e.g. minus-strand
        # stop-loss variants like FOXP3 c.1293_1294del, anchored 2 bp
        # downstream of the last CDS exon) still get fetched. The
        # range-scan in _project_genomic_to_aa then walks the REF span to
        # find the CDS base. 24 bp covers ~99 % of clinically-reported
        # indels. Including alt-iso ranges lets us surface variants that
        # fall inside an alt exon the canonical transcript skips
        # entirely — they get pos=0 (canonical sentinel) + an entry in
        # `isoforms` for whichever alt iso encodes the residue.
        INDEL_PAD = 24
        intervals: list[tuple[int, int]] = []
        for c in canon_map.cds:
            intervals.append((max(1, c.g_start - INDEL_PAD), c.g_end + INDEL_PAD))
        for alt_ensp in alt_ensps.get(gene_u, []):
            a_map = p2x.get(alt_ensp)
            if a_map is None or a_map.chrom != canon_map.chrom:
                continue
            for c in a_map.cds:
                intervals.append((max(1, c.g_start - INDEL_PAD), c.g_end + INDEL_PAD))
        intervals.sort()
        merged_iv: list[tuple[int, int]] = []
        for s, e in intervals:
            if merged_iv and s <= merged_iv[-1][1] + 1:
                merged_iv[-1] = (merged_iv[-1][0], max(merged_iv[-1][1], e))
            else:
                merged_iv.append((s, e))
        params: list = [canon_map.chrom]
        clauses: list[str] = []
        for s, e in merged_iv:
            clauses.append("(pos BETWEEN ? AND ?)")
            params.extend([s, e])
        sql = (
            "SELECT chrom, pos, ref, alt, source, clinvar_id, cosmic_id, "
            "rsid, gene, name, type, consequence, significance, "
            "review_status, origin, af, ac, an, popmax_af, pmids, "
            "transcript_id, hgvs_p, so_term, is_canonical "
            "FROM variants WHERE chrom = ? AND (" + " OR ".join(clauses) + ")"
        )
        rows = cur.execute(sql, params).fetchall()
        if not rows:
            continue

        # Load the TF JSON so we can look up isoform → ENSP + isoform name.
        tf = json.loads(fp.read_text())
        iso_name_by_ensp: dict[str, str] = {}
        for iso in tf.get("isoforms") or []:
            ie = iso.get("proteinId")
            if ie:
                iso_name_by_ensp[ie] = iso.get("name") or ie
        # iso_name → ENST for COSMIC per-transcript matching. COSMIC ships
        # one row per (variant × ENST), so we can attach iso-specific HGVS.p
        # + SO-term by matching this map to the row's `transcript_id`.
        iso_name_by_enst: dict[str, str] = {}
        for t in tf.get("transcripts") or []:
            enst = (t.get("ensemblId") or "").split(".", 1)[0]
            iso_name = t.get("proteinIsoform")
            if enst and iso_name:
                iso_name_by_enst[enst] = iso_name

        # Pre-group COSMIC rows by (chrom,pos,ref,alt). The primary lollipop
        # is driven by the canonical-transcript row when available; the
        # non-canonical rows feed `isoConsequence` + `isoHgvsP` so iso view
        # shows the right HGVS.p for each iso (e.g. T162A on canonical,
        # T141A on TP53-204).
        cosmic_per_variant: dict[tuple[str, int, str, str], list] = defaultdict(list)
        for r in rows:
            if r["source"] == "cosmic":
                cosmic_per_variant[
                    (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                ].append(r)
        # Track which COSMIC rows we've already emitted so we don't double-
        # count when iterating below.
        emitted_cosmic: set[tuple[str, int, str, str]] = set()

        variants: list[dict] = []
        for r in rows:
            # COSMIC dedup: ship one variant entry per (chrom,pos,ref,alt)
            # regardless of how many per-transcript rows we ingested. The
            # primary `pos` + `aaRef` + `aaAlt` come from the canonical
            # COSMIC row (IS_CANONICAL=y) when present, falling back to the
            # current row. The non-canonical rows feed `isoConsequence` +
            # `isoHgvsP` for iso-aware display.
            if r["source"] == "cosmic":
                key = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                if key in emitted_cosmic:
                    continue
                emitted_cosmic.add(key)
                group = cosmic_per_variant.get(key, [r])
                # Prefer the canonical-marked row; falls back to the
                # first row of the group otherwise.
                primary = next(
                    (g for g in group if (g["is_canonical"] or 0)),
                    group[0],
                )
                r = primary
            ref_len = len(r["ref"] or "") or 1
            aa_pos = _project_genomic_to_aa(
                int(r["pos"]), canon_map, ref_len=ref_len,
            )
            cons, aa_ref, aa_alt, hgvs_aa = _parse_hgvs_consequence(r["name"] or "")
            if hgvs_aa is not None and aa_pos is None:
                aa_pos = hgvs_aa
            # When canonical projection fails we may still have a hit in
            # an alt iso (variant lives in an exon canonical skips). We
            # don't continue here — the iso projection below decides.
            iso_only = aa_pos is None
            # gnomAD rows carry no HGVS — they're sliced straight from the
            # VCF without VEP. Infer a coarse consequence from ref/alt
            # length: SNV → likely missense, in-frame indel → "deletion"
            # (other bucket), frameshift otherwise. Better than tagging
            # everything as "Other" which lost the colour signal completely.
            if r["source"] == "gnomad" and not cons:
                ref_s = str(r["ref"] or "")
                alt_s = str(r["alt"] or "")
                if len(ref_s) == 1 and len(alt_s) == 1:
                    cons = "missense"
                elif (len(alt_s) - len(ref_s)) % 3 == 0:
                    cons = "deletion" if len(alt_s) < len(ref_s) else "duplication"
                else:
                    cons = "frameshift"
                # Surface the original AA letter from canonical seq for the
                # SNV case so the side panel shows e.g. "R175?" instead of
                # "?175?". TF JSON's `sequence` is the Ensembl-canonical AA
                # string (set in build_tf_records.py). Skip for iso-only
                # variants (no canonical aa_pos to index into).
                if cons == "missense" and aa_pos is not None:
                    canon_seq = str(tf.get("sequence") or "")
                    if canon_seq and 1 <= aa_pos <= len(canon_seq):
                        aa_ref = canon_seq[aa_pos - 1]
            # Map consequence → existing MutationType union the UI expects.
            mtype = {
                "missense": "Missense",
                "nonsense": "Nonsense",
                "synonymous": "Synonymous",
                "frameshift": "Frameshift",
                "deletion": "Frameshift",
                "duplication": "Frameshift",
            }.get(cons, "Other")
            # Source bucket: match the existing UI toggle names.
            src_bucket = {"clinvar": "ClinVar", "cosmic": "COSMIC", "gnomad": "gnomAD"}[r["source"]]
            # Stable display id: prefer HGVS p.<short> like `R175H` when we
            # have it; else COSMIC id; else rsID; else ClinVar ID.
            if aa_ref and aa_alt and aa_pos:
                vid = f"{aa_ref}{aa_pos}{aa_alt}"
            elif r["cosmic_id"]:
                vid = r["cosmic_id"]
            elif r["rsid"]:
                vid = r["rsid"]
            elif r["clinvar_id"]:
                vid = f"CV:{r['clinvar_id']}"
            else:
                vid = f"{r['chrom']}:{r['pos']}{r['ref']}>{r['alt']}"

            # `count` drives the canvas log-scale lollipop height. For COSMIC
            # we use the stored sample-occurrence count (in the `ac` column).
            # For gnomAD we use allele-count (rare variants get short sticks,
            # common variants tall — visually consistent with the canvas).
            count_val = 1
            if r["source"] == "cosmic" and r["ac"]:
                count_val = max(1, int(r["ac"]))
            elif r["source"] == "gnomad" and r["ac"]:
                count_val = max(1, int(r["ac"]))

            # Rich metadata — optional on the Mutation interface but consumed
            # by the upcoming variant side-panel.
            cv = None
            if r["source"] == "clinvar":
                cv = {
                    "id": r["clinvar_id"],
                    "name": r["name"],
                    "significance": r["significance"],
                    "significanceTag": _significance_label(r["significance"]),
                    "reviewStatus": r["review_status"],
                    "origin": r["origin"],
                }

            # Project genomic pos onto every alt isoform whose CDS encodes it.
            iso_aas: dict[str, int] = {}
            for alt_ensp in alt_ensps.get(gene_u, []):
                a_map = p2x.get(alt_ensp)
                if a_map is None or a_map.chrom != canon_map.chrom:
                    continue
                a_pos = _project_genomic_to_aa(int(r["pos"]), a_map, ref_len=ref_len)
                if a_pos is not None:
                    iso_aas[iso_name_by_ensp.get(alt_ensp, alt_ensp)] = a_pos

            # Per-iso COSMIC annotations: HGVS.p + derived consequence
            # keyed by iso name. Built from the non-primary rows in the
            # COSMIC group (one row per transcript). The UI's iso view
            # consumes this so the side panel shows "p.T141A · missense"
            # on TP53-204 even though canonical reads "p.T162A".
            #
            # We derive consequence from the HGVS.p string itself (via
            # _parse_hgvs_consequence) rather than COSMIC's SO_TERM —
            # SO_TERM is the variant CLASS ("SNV" / "insertion" /
            # "deletion") not the protein consequence. The CLASS gives
            # users a misleading colour signal in the lollipop tooltip.
            iso_cons: dict[str, str] = {}
            iso_hgvs_p: dict[str, str] = {}
            if r["source"] == "cosmic":
                group = cosmic_per_variant.get(
                    (r["chrom"], int(r["pos"]), r["ref"], r["alt"]),
                    [r],
                )
                for g in group:
                    enst = g["transcript_id"]
                    iso_name = iso_name_by_enst.get(enst or "")
                    if not iso_name:
                        continue
                    hp = g["hgvs_p"]
                    if hp:
                        iso_hgvs_p[iso_name] = hp
                        cons_iso, _, _, _ = _parse_hgvs_consequence(hp)
                        if cons_iso:
                            iso_cons[iso_name] = cons_iso

            # Iso-only variants: canonical missed but at least one alt
            # iso encodes this position. Emit with pos=0 sentinel — the
            # canonical lollipop chart filters these out (pos < 1), the
            # iso chart picks them up via `isoforms[isoName]`.
            if iso_only:
                if not iso_aas:
                    continue
                aa_pos = 0

            # gnomAD payload (allele frequency / popmax) for the side panel.
            gn = None
            if r["source"] == "gnomad":
                gn = {
                    "af": r["af"],
                    "ac": r["ac"],
                    "an": r["an"],
                    "popmaxAf": r["popmax_af"],
                }

            variants.append({
                # Required by Mutation interface
                "id": vid,
                "pos": aa_pos,
                "type": mtype,
                "count": count_val,
                "source": src_bucket,
                # Rich extras
                "chrom": r["chrom"], "genomicPos": int(r["pos"]),
                "ref": r["ref"], "alt": r["alt"],
                "aaRef": aa_ref, "aaAlt": aa_alt,
                "consequence": cons or None,
                "rsid": r["rsid"] or None,
                "clinvar": cv,
                "gnomad": gn,
                "cosmicId": r["cosmic_id"] if r["source"] == "cosmic" else None,
                "isoforms": iso_aas,
                # Per-iso COSMIC consequence + HGVS.p (empty for ClinVar/gnomAD).
                **(
                    {"isoConsequence": iso_cons, "isoHgvsP": iso_hgvs_p}
                    if (iso_cons or iso_hgvs_p)
                    else {}
                ),
            })

        # Sort by AA position for stable rendering.
        variants.sort(key=lambda v: (v["pos"], v.get("genomicPos") or 0))

        out = {
            "symbol": symbol,
            "canonicalEnsp": c_ensp,
            "variants": variants,
        }
        out_path = args.out_dir / f"{symbol}.json"
        out_path.write_text(json.dumps(out, separators=(",", ":")))
        counts[symbol] = len(variants)
        n_emitted += 1
        if variants:
            n_with_variants += 1
            n_total_variants += len(variants)

        # Variants are NOT inlined into the TF mock JSON — TP53 has 8k+ rows,
        # which ballooned tfs/TP53.json to 6 MB and made dev-server SSR + React
        # hydration take 40+ seconds per page load. Instead the TF JSON keeps
        # `mutations: []` and MutationTrack fetches /mock/variants/{symbol}.json
        # client-side after first paint.
        if tf.get("mutations") != []:
            tf["mutations"] = []
            fp.write_text(json.dumps(tf, separators=(",", ":")))

    # Compact {SYMBOL: count} index for the Protein Overview stat. The full
    # per-TF variants files are too large to ship in the source repo / build
    # tarball (served from R2 at runtime), so the page can't fs-read them at
    # build time. This index IS small enough to commit, so it's read server-
    # side at build and the count renders into the static HTML.
    counts_path = args.out_dir.parent / "variant_counts.json"
    counts_path.write_text(
        json.dumps(dict(sorted(counts.items())), separators=(",", ":"))
    )

    print()
    print(f"  per-TF variant files emitted: {n_emitted}")
    print(f"  TFs with at least 1 variant:  {n_with_variants}")
    print(f"  total variants projected:     {n_total_variants:,}")
    print(f"  → {args.out_dir}")
    print(f"  variant count index:          {counts_path}")
    con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
