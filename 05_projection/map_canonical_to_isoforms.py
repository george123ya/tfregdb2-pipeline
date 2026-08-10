"""Map canonical-frame curated domains onto every alt isoform's frame.

For each curated effector / DBD call (canonical-frame coords s..e in the
Ensembl-canonical ENSP), emit one row per alt isoform of the same gene with
the equivalent iso-frame coords. Two paths:

  1. **Genomic (prot2exon)** — preferred when both the canonical and alt ENSPs
     have an entry in prot2exon's isoform_structure.tsv. We translate the
     canonical AA range into a genomic interval via the canonical CDS map,
     then intersect that interval with each alt's CDS exons to find the iso
     AA positions encoded by the same genome region. Real biological splicing
     drops out for free — if the alt skips the exon carrying the canonical
     domain, the row gets no projection.

  2. **Sequence alignment fallback** — when prot2exon has no entry for the alt
     ENSP (rare; only when the BED didn't include it), we global-pairwise-
     align the canonical sequence vs the iso sequence (biopython Align) and
     walk the alignment to map canonical (s,e) → iso (s',e'). Identity from
     the alignment goes into `alignment_identity`.

Output schema matches the legacy `TF_effector_data_with_derivadas_and_alignment_ideal.csv`
so build_tf_mocks can swap it in via SRC_EFFECTOR_CSV — but only the columns
the build actually reads downstream:
    Gene name (HGNC), Translation ID, Source, Domain_origin,
    From, To, alignment_identity, alignment_coverage.

`Source` encodes the lineage of the canonical call so downstream type-resolution
keeps working: `{gene}_{Source_Database}_Row{row_idx}` (matches the legacy
format). `Domain_origin` is `ORIGINAL_X` for the canonical row, `DERIVADA_X`
for genomic-mapped alt isos, `ALIENEAMIENTO_X` for sequence-aligned ones —
lineage letter X is unique per (gene × source row).
"""

from __future__ import annotations

import argparse
import string
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices


@dataclass
class CdsExon:
    """One CDS chunk of a protein: encodes aa_start..aa_end at genomic
    [g_start, g_end] on the given strand. + strand: AA increases as genomic
    coord increases; - strand: AA increases as genomic coord decreases.
    All ranges are inclusive."""
    aa_start: int
    aa_end: int
    g_start: int
    g_end: int


@dataclass
class EnspMap:
    strand: str
    cds: list[CdsExon]   # in transcript (5'→3') order
    seq_len: int


def _parse_prot2exon(p2x_tsv: Path) -> dict[str, EnspMap]:
    """Build per-ENSP CDS exon list from prot2exon's isoform_structure.tsv.

    The file has one row per (input_query × feature) — many duplicate rows for
    the same ENSP because each canonical row's BED entry produces a separate
    block. We pick one block per ENSP (the first we see) and keep only CDS
    features ordered by transcript position."""
    print(f"Reading prot2exon table {p2x_tsv}...")
    df = pd.read_csv(p2x_tsv, sep="\t", low_memory=False)
    df = df[df["feature_type"] == "CDS"].copy()
    df["aa_start_encoded"] = pd.to_numeric(df["aa_start_encoded"], errors="coerce")
    df["aa_end_encoded"] = pd.to_numeric(df["aa_end_encoded"], errors="coerce")
    df = df.dropna(subset=["aa_start_encoded", "aa_end_encoded"])

    out: dict[str, EnspMap] = {}
    for ensp, sub in df.groupby("protein_id"):
        # Use the FIRST input_id's block — all blocks for the same ENSP share
        # the same CDS structure (only `overlaps_domain` columns differ).
        first_input = sub["input_id"].iloc[0]
        block = sub[sub["input_id"] == first_input]
        strand = str(block["strand"].iloc[0])
        block = block.sort_values("feature_order_transcript")
        cds_raw = [
            (
                int(r["aa_start_encoded"]),
                int(r["aa_end_encoded"]),
                int(r["feature_genomic_start"]),
                int(r["feature_genomic_end"]),
            )
            for _, r in block.iterrows()
        ]
        if not cds_raw:
            continue
        # prot2exon emits a boundary AA on both sides of a split codon — e.g.
        # CDS_1 reports aa 1-60 and CDS_2 reports aa 60-113 (AA 60's codon
        # spans the junction). Without disambiguating, an iso-frame projection
        # of canonical RD 84-183 onto RARA-206 returns iso aa 60-86 instead
        # of 61-86, and the boundary AA later triggers a chip that visibly
        # bridges across a 20 kb intron in the genomic view. Give each
        # boundary AA to the earlier exon by bumping the later exon's aa_start.
        cds: list[CdsExon] = []
        prev_end = 0
        for a_s, a_e, g_s, g_e in cds_raw:
            adj_s = max(a_s, prev_end + 1)
            if a_e < adj_s:
                continue
            cds.append(CdsExon(aa_start=adj_s, aa_end=a_e, g_start=g_s, g_end=g_e))
            prev_end = a_e
        if not cds:
            continue
        out[ensp] = EnspMap(strand=strand, cds=cds, seq_len=max(c.aa_end for c in cds))
    print(f"  ENSPs with CDS map: {len(out)}")
    return out


def _aa_to_genomic_range(ensp_map: EnspMap, s: int, e: int) -> list[tuple[int, int]]:
    """Genomic intervals (one per CDS chunk that overlaps AA range s..e)."""
    out: list[tuple[int, int]] = []
    plus = ensp_map.strand == "+"
    for cds in ensp_map.cds:
        if e < cds.aa_start or s > cds.aa_end:
            continue
        ovl_aa_s = max(s, cds.aa_start)
        ovl_aa_e = min(e, cds.aa_end)
        # Each AA = 3 nt. AA offset within CDS = ovl - cds.aa_start.
        if plus:
            g_lo = cds.g_start + (ovl_aa_s - cds.aa_start) * 3
            g_hi = cds.g_start + (ovl_aa_e - cds.aa_start + 1) * 3 - 1
        else:
            g_hi = cds.g_end - (ovl_aa_s - cds.aa_start) * 3
            g_lo = cds.g_end - (ovl_aa_e - cds.aa_start + 1) * 3 + 1
        out.append((g_lo, g_hi))
    return out


def _genomic_range_to_aa(ensp_map: EnspMap, g_lo: int, g_hi: int) -> tuple[int, int] | None:
    """Iso AA positions whose codons fall inside genomic [g_lo, g_hi]. Returns
    (min, max) over all iso AAs that overlap, or None if no overlap."""
    plus = ensp_map.strand == "+"
    iso_positions: list[int] = []
    for cds in ensp_map.cds:
        if cds.g_end < g_lo or cds.g_start > g_hi:
            continue
        ovl_g_lo = max(g_lo, cds.g_start)
        ovl_g_hi = min(g_hi, cds.g_end)
        if plus:
            aa_lo = cds.aa_start + (ovl_g_lo - cds.g_start) // 3
            aa_hi = cds.aa_start + (ovl_g_hi - cds.g_start) // 3
        else:
            aa_lo = cds.aa_start + (cds.g_end - ovl_g_hi) // 3
            aa_hi = cds.aa_start + (cds.g_end - ovl_g_lo) // 3
        iso_positions.extend([min(aa_lo, aa_hi), max(aa_lo, aa_hi)])
    if not iso_positions:
        return None
    return min(iso_positions), max(iso_positions)


def _align_canonical_to_iso(
    canon_seq: str, iso_seq: str, s: int, e: int, aligner: Align.PairwiseAligner,
) -> tuple[int, int, float, float] | None:
    """Sequence-alignment fallback: global align, walk alignment to find iso
    positions corresponding to canonical [s..e]. Returns (iso_s, iso_e,
    identity_pct, coverage_pct) or None on no overlap."""
    if not canon_seq or not iso_seq:
        return None
    # The aligner returns an iterator of equally-optimal alignments. Calling
    # len() can overflow when match/mismatch scoring leaves many ties — just
    # take the first one (BLOSUM62 typically produces a unique optimum).
    try:
        aln = next(iter(aligner.align(canon_seq, iso_seq)))
    except (StopIteration, OverflowError, ValueError):
        return None
    # aln.coordinates: shape (2, N+1) — start/end of each aligned block in
    # canonical (row 0) and iso (row 1). Walk through blocks to build a
    # canon-pos → iso-pos lookup just for positions in [s, e].
    coords = aln.coordinates
    canon_to_iso: dict[int, int] = {}
    n_blocks = coords.shape[1] - 1
    for b in range(n_blocks):
        c0, c1 = int(coords[0][b]), int(coords[0][b + 1])
        i0, i1 = int(coords[1][b]), int(coords[1][b + 1])
        if c0 == c1 or i0 == i1:
            continue
        # Both advance: matched block (or substitution). Map 1:1.
        block_len = min(c1 - c0, i1 - i0)
        for k in range(block_len):
            canon_pos = c0 + k + 1  # 1-based
            iso_pos = i0 + k + 1
            if s <= canon_pos <= e:
                canon_to_iso[canon_pos] = iso_pos
    if not canon_to_iso:
        return None
    iso_s = min(canon_to_iso.values())
    iso_e = max(canon_to_iso.values())
    # Identity: matched chars in [s, e] / (e - s + 1).
    matched = sum(
        1
        for cp, ip in canon_to_iso.items()
        if canon_seq[cp - 1] == iso_seq[ip - 1]
    )
    domain_len = e - s + 1
    identity = 100.0 * matched / domain_len if domain_len else 0.0
    coverage = 100.0 * len(canon_to_iso) / domain_len if domain_len else 0.0
    return iso_s, iso_e, round(identity, 1), round(coverage, 1)


def _lineage_letters(n: int) -> list[str]:
    """A, B, ..., Z, AA, AB, ... up to n labels."""
    letters = list(string.ascii_uppercase)
    out: list[str] = []
    i = 0
    while len(out) < n:
        a, b = divmod(i, 26)
        out.append((letters[a - 1] if a > 0 else "") + letters[b])
        i += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unified", type=Path, required=True,
                    help="Unified_Valid_Domains.xlsx (canonical-frame domains).")
    ap.add_argument("--ids-csv", type=Path, required=True,
                    help="TF_completeIDs_with_ensembl_canonical.csv")
    ap.add_argument("--prot2exon", type=Path, required=True,
                    help="prot2exon results_isoforms/isoform_structure.tsv")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output canonical-to-iso mapping CSV.")
    args = ap.parse_args()

    # ---- Load canonical metadata + sequences --------------------------------
    print(f"Loading TF_completeIDs {args.ids_csv}...")
    ids = pd.read_csv(args.ids_csv)
    ids["Translation_Clean"] = ids["Translation ID"].astype(str).str.split(".").str[0]
    ids["Gene_Upper"] = ids["Gene name (HGNC)"].astype(str).str.strip().str.upper()
    ids = ids[ids["Translation_Clean"].str.startswith("ENSP", na=False)].copy()

    # Pick canonical ENSP per gene: Ensembl canonical → APPRIS PRINCIPAL:1 →
    # PRINCIPAL → longest. Match build_tf_mocks's chain.
    canon_ensp: dict[str, str] = {}
    ensp_seq: dict[str, str] = {}
    ensp_gene: dict[str, str] = {}
    for ensp, row in ids.groupby("Translation_Clean").first().iterrows():
        seq = str(row.get("Sequence aa") or "").strip()
        if seq and seq.lower() != "nan":
            ensp_seq[ensp] = seq
        ensp_gene[ensp] = str(row["Gene_Upper"])

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
    print(f"  Genes with canonical ENSP: {len(canon_ensp)}")

    # Per-gene alt-ENSPs: every translation ID for that gene minus the canonical.
    gene_to_alt_ensps: dict[str, list[str]] = defaultdict(list)
    for gene, grp in ids.groupby("Gene_Upper"):
        canon = canon_ensp.get(gene)
        for ensp in grp["Translation_Clean"].dropna().unique():
            if ensp != canon and ensp in ensp_seq:
                gene_to_alt_ensps[gene].append(ensp)

    # ---- Load prot2exon CDS maps --------------------------------------------
    p2x = _parse_prot2exon(args.prot2exon)

    # ---- Load curated canonical-frame domains -------------------------------
    print(f"Loading unified {args.unified}...")
    u = pd.read_excel(args.unified)
    u = u[u["Final_Human_UNIPROT_ID"].notna()].copy()
    # Re-bin every row under the canonical HGNC gene that owns its UniProt —
    # same logic as build_tf_mocks. This moves the TEAD1-misrouted "TEF" rows
    # (Original_Gene=TEF, UniProt=P28347) onto TEAD1 where the projection
    # actually landed, so the canonical→iso mapping doesn't propagate the
    # curator's gene-name error into 7+ iso rows per misrouted entry.
    up_to_gene: dict[str, str] = {}
    for ensp, gene in ensp_gene.items():
        if not gene or gene == "NAN":
            continue
        # Pick canonical UniProt for this gene from the canonical ENSP row.
        c_ensp = canon_ensp.get(gene)
        if not c_ensp:
            continue
        canon_row = ids[ids["Translation_Clean"] == c_ensp]
        if canon_row.empty:
            continue
        up = str(canon_row["UniProt_ID"].iloc[0] or "").split("-")[0].strip().upper()
        if up and up not in up_to_gene:
            up_to_gene[up] = gene
    u["UP_Clean"] = u["Final_Human_UNIPROT_ID"].astype(str).str.split("-").str[0].str.strip().str.upper()
    u["_Gene_Norm"] = u["Original_Gene"].map(_normalize_for_lookup)
    n_moved = 0
    def _rebind(row):
        nonlocal n_moved
        up = row["UP_Clean"]
        if up in up_to_gene:
            new_g = up_to_gene[up]
            if new_g != row["_Gene_Norm"]:
                n_moved += 1
            return new_g
        return row["_Gene_Norm"]
    u["Gene_Upper"] = u.apply(_rebind, axis=1)
    print(f"  unified rows: {len(u)} (re-binned {n_moved} rows by UniProt)")

    # ---- Build the mapping ---------------------------------------------------
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    out_rows: list[dict] = []
    n_p2x = n_aln = n_skip = 0
    # Generate one lineage letter per source row, per gene.
    per_gene_counter: dict[str, int] = defaultdict(int)

    for idx, r in u.iterrows():
        gene = str(r.get("Gene_Upper") or "").strip().upper()
        if not gene or gene == "NAN":
            continue
        c_ensp = canon_ensp.get(gene)
        if not c_ensp:
            n_skip += 1
            continue
        coords = str(r.get("Final_Human_DomainCoords") or "").strip()
        if "-" not in coords:
            continue
        try:
            s, e = [int(x) for x in coords.split("-")[:2]]
        except ValueError:
            continue

        src_db = str(r.get("Source_Database") or "")
        row_idx = r.get("Original_Row_Index") or idx
        # Source string matches the legacy format the effector csv uses so
        # downstream type-resolution lookups still work.
        source = f"{gene}_{src_db}_Row{int(row_idx)}"

        # Lineage letter — one per (gene × source row).
        lin = _lineage_letters(per_gene_counter[gene] + 1)[-1]
        per_gene_counter[gene] += 1

        # ORIGINAL row for the canonical iso (coords as-is).
        out_rows.append({
            "Gene name (HGNC)": gene,
            "Translation ID": c_ensp,
            "From": s,
            "To": e,
            "Source": source,
            "Domain_origin": f"ORIGINAL_{lin}",
            "alignment_identity": None,
            "alignment_coverage": None,
        })

        canon_seq = ensp_seq.get(c_ensp, "")
        canon_p2x = p2x.get(c_ensp)

        for alt in gene_to_alt_ensps.get(gene, []):
            alt_p2x = p2x.get(alt)
            iso_coords: tuple[int, int] | None = None
            origin = None
            identity = None
            coverage = None
            # 1) Genomic path: both ENSPs in prot2exon.
            if canon_p2x and alt_p2x and canon_p2x.strand == alt_p2x.strand:
                g_ranges = _aa_to_genomic_range(canon_p2x, s, e)
                iso_aas: list[int] = []
                for g_lo, g_hi in g_ranges:
                    hit = _genomic_range_to_aa(alt_p2x, g_lo, g_hi)
                    if hit:
                        iso_aas.extend(hit)
                if iso_aas:
                    iso_coords = (min(iso_aas), max(iso_aas))
                    origin = f"DERIVADA_{lin}"
                    n_p2x += 1
            # 2) Sequence-alignment fallback.
            if iso_coords is None:
                iso_seq = ensp_seq.get(alt, "")
                ali = _align_canonical_to_iso(canon_seq, iso_seq, s, e, aligner)
                if ali is not None:
                    iso_s, iso_e, identity, coverage = ali
                    iso_coords = (iso_s, iso_e)
                    origin = f"ALIENEAMIENTO_{lin}"
                    n_aln += 1
            if iso_coords is None:
                continue

            out_rows.append({
                "Gene name (HGNC)": gene,
                "Translation ID": alt,
                "From": iso_coords[0],
                "To": iso_coords[1],
                "Source": source,
                "Domain_origin": origin,
                "alignment_identity": identity,
                "alignment_coverage": coverage,
            })

    out_df = pd.DataFrame(out_rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False)
    print()
    print(f"Rows written: {len(out_df)}")
    print(f"  ORIGINAL:      {(out_df['Domain_origin'].str.startswith('ORIGINAL')).sum()}")
    print(f"  DERIVADA:      {n_p2x}")
    print(f"  ALIENEAMIENTO:  {n_aln}")
    print(f"  unified rows skipped (no canonical ENSP for gene): {n_skip}")
    print(f"Output: {args.out}")
    return 0


def _normalize_for_lookup(g) -> str:
    """Strip _HUMAN suffix + upper-case (matches build_tf_mocks normalisation)."""
    s = str(g or "").strip().upper()
    if s.endswith("_HUMAN"):
        s = s[: -len("_HUMAN")]
    return s


if __name__ == "__main__":
    sys.exit(main())
