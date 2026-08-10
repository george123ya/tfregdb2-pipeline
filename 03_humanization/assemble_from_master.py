"""Assemble Unified_Valid_Domains.xlsx from the user's curated UpdatedTable.xlsx.

The master workbook already contains every BLAST projection done upstream:

  • HumanizedTable        — 1,053 humanized rows from CleanTableFiltered → 838
                            with Final_Human_UNIPROT_ID after BLAST.
  • TFRegDB1_Humanized    — 924 legacy TFRegDB1 records with BLAST projections.
  • Studies_Dataframe     — 3,548 rows pulled from the 10 study-specific
                            sheets (Tycko 2020/2024, DelRosso, Klaus, Arnold,
                            Alerasool); 3,136 have Final_Human_UNIPROT_ID.

We keep only rows that successfully projected onto a human canonical
(Final_Human_UNIPROT_ID present + non-failure Mapping_Status), normalize
column names to the schema build_tf_mocks.py expects, and tag each row's
Source_Database with the fine-grained study label (Provenance) — this is
what makes SOURCE_TYPE_HINT lookups resolve downstream.

Output columns required by build_tf_mocks.py:
    Final_Human_UNIPROT_ID, Final_Human_DomainCoords, Final_Human_ENSP_ID,
    Final_Human_Domain_Sequence, Final_Identity, Mapping_Status,
    Original_Gene, Source_Database, PUBMED_ID,
    Domain type / DomainType, DomainCoordinates / Coordinates,
    BLAST_UNIPROT_ID, BLAST_DomainCoordinates, BLAST_Identity_Pct.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def _normalize_gene(g) -> str:
    """Trim, uppercase, strip `_HUMAN` suffix (UniProt-style entry names)."""
    s = str(g or "").strip().upper()
    if s.endswith("_HUMAN"):
        s = s[: -len("_HUMAN")]
    return s


def _is_success(status) -> bool:
    s = str(status or "").strip().lower()
    return s.startswith("successful") or s.startswith("confirmed human")


def _from_humanized(ht: pd.DataFrame) -> pd.DataFrame:
    df = ht.copy()
    df["Source_Database"] = "HumanizedTable"
    df["Original_Gene"] = df["Gene"].map(_normalize_gene)
    # PUBMED_ID is already the column name on this sheet.
    df["UNIPROT ID"] = df.get("UNIPROT ID")
    df["UNIPROT SEQ"] = df.get("UNIPROT SEQ")
    df["UNIPROT LENGTH"] = df.get("UNIPROT LENGTH")
    return df


def _from_tfregdb1(tf: pd.DataFrame) -> pd.DataFrame:
    df = tf.copy()
    df["Source_Database"] = "TFRegDB1"
    df["Original_Gene"] = df["TF name"].map(_normalize_gene)
    # Build a PUBMED_ID column from Reference (PMID) so the downstream PMID
    # coalesce (which reads Reference (PMID) OR PUBMED_ID) gets a hit either way.
    df["PUBMED_ID"] = df.get("Reference (PMID)")
    return df


def _from_studies(sd: pd.DataFrame) -> pd.DataFrame:
    df = sd.copy()
    # Provenance values like "DelRosso_AD" become Source_Database — matches the
    # legacy fine-grained labels that SOURCE_TYPE_HINT keys off.
    df["Source_Database"] = df["Provenance"].astype(str)
    df["Original_Gene"] = df["Original_Gene"].map(_normalize_gene)
    # Studies have no PMID per row — leave blank; build_tf_mocks falls back to
    # PAPER_TITLE_FALLBACK for the study-paper titles.
    df["PUBMED_ID"] = ""
    # Studies don't carry the curator's free-form DomainCoordinates; the
    # human-frame coords are the only ones — use them under DomainCoordinates so
    # the coord-based effector-type lookup in build_tf_mocks finds them.
    df["DomainCoordinates"] = df["Final_Human_DomainCoords"]
    df["DomainType"] = df["Domain type"]
    return df


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", type=Path, required=True,
                    help="Path to the curated UpdatedTable.xlsx.")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output Unified_Valid_Domains.xlsx path.")
    args = ap.parse_args()

    print(f"Loading {args.master}...")
    ht = pd.read_excel(args.master, sheet_name="HumanizedTable")
    tf = pd.read_excel(args.master, sheet_name="TFRegDB1_Humanized")
    sd = pd.read_excel(args.master, sheet_name="Studies_Dataframe")
    print(f"  HumanizedTable:        {len(ht):>5} rows")
    print(f"  TFRegDB1_Humanized:    {len(tf):>5} rows")
    print(f"  Studies_Dataframe:     {len(sd):>5} rows")

    parts = []
    for label, frame in (("HumanizedTable", _from_humanized(ht)),
                         ("TFRegDB1", _from_tfregdb1(tf)),
                         ("Studies", _from_studies(sd))):
        before = len(frame)
        # Filter to rows that actually projected onto a human canonical.
        frame = frame[frame["Final_Human_UNIPROT_ID"].notna()].copy()
        # Drop the "Failed:" rows even if they have a UniProt — keep only
        # Successful BLAST Hit (>95%) and Confirmed Human variants.
        if "Mapping_Status" in frame.columns:
            frame = frame[frame["Mapping_Status"].map(_is_success)].copy()
        after = len(frame)
        print(f"  {label:>16} kept {after}/{before}")
        parts.append(frame)

    unified = pd.concat(parts, ignore_index=True, sort=False)
    print(f"Unified rows: {len(unified)}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    unified.to_excel(args.out, index=False)
    print(f"Written → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
