"""Stage 0: hydrate HumanizedTable rows with UniProt + sequence by joining
on gene symbol → canonical Ensembl protein.

The source xlsx ships ~3,165 of 3,428 HumanizedTable rows with empty
`UNIPROT ID` / `UNIPROT SEQ` cells. Without those, stage 1 can't slice the
domain sequence out, BLAST has nothing to query, and the projection drops
the row.

This stage uses `TF_completeIDs_with_ensembl_canonical.csv` as the
source-of-truth: every gene that appears in the curated tables is keyed to
its canonical ENSP (per Ensembl's `Ensembl canonical` flag, with APPRIS
PRINCIPAL:1 + longest as fallback ties), and that ENSP carries the
UniProt accession + protein sequence we need.

Behaviour:
  • Human rows → looked up by `Gene` (HGNC symbol). The row's UNIPROT ID +
    UNIPROT SEQ get filled from the canonical ENSP. DomainCoordinates are
    untouched.
  • Non-human rows (Mouse, Rat, Xenopus, …) → can't be projected without
    fetching the animal sequence from UniProt. We tag them as
    `NEEDS_UNIPROT_FETCH` in the failure log and skip enrichment — they
    can be re-run later through an API-aware step.

The output is `hydrated_HumanizedTable.xlsx` with the same schema as the
input but with rescued sequences. TFRegDB1 is mostly hydrated already so we
leave it alone.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# Treat any Species string containing one of these substrings as human.
# A row tagged "Human and Mouse" still gets the human canonical mapping —
# the domain coords are validated against the human reference.
HUMAN_TOKENS = ("human", "homo sapiens")


def _is_human(species: str) -> bool:
    if not species:
        return False
    s = str(species).lower()
    return any(tok in s for tok in HUMAN_TOKENS)


def _build_canonical_map(csv_path: Path) -> dict[str, dict]:
    """Gene symbol → canonical ENSP record. Priority within a gene:
       1. Row with `Ensembl canonical == "Canonical"`
       2. Row whose APPRIS Annotation starts with "PRINCIPAL:1"
       3. Longest protein
    """
    df = pd.read_csv(csv_path)
    df["Translation_Clean"] = df["Translation ID"].astype(str).str.split(".").str[0]
    df["Gene_Upper"] = df["Gene name (HGNC)"].astype(str).str.strip().str.upper()
    df["UniProt_Clean"] = df["UniProt_ID"].astype(str).str.split("-").str[0].str.strip()

    def pick(group: pd.DataFrame) -> pd.Series:
        canon = group[group["Ensembl canonical"].astype(str).str.strip() == "Canonical"]
        if len(canon):
            return canon.iloc[0]
        p1 = group[group["APPRIS Annotation"].astype(str).str.contains("PRINCIPAL:1", na=False)]
        if len(p1):
            return p1.iloc[0]
        return group.sort_values("Protein length", ascending=False).iloc[0]

    out: dict[str, dict] = {}
    for g, grp in df.dropna(subset=["Translation_Clean", "Sequence aa"]).groupby("Gene_Upper"):
        if not g or g == "NAN":
            continue
        row = pick(grp)
        out[g] = {
            "ensp": row["Translation_Clean"],
            "uniprot": row["UniProt_Clean"],
            "sequence": str(row["Sequence aa"]).strip().upper(),
            "length": int(row.get("Protein length") or 0),
        }
    return out


def _coords_in_range(coords: str, length: int) -> bool:
    """Quick check: do the DomainCoordinates fit inside `length` aa?"""
    if not coords or pd.isna(coords):
        return False
    raw = str(coords).strip().replace("–", "-").replace("—", "-").lstrip(":")
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", raw)
    if not m:
        return False
    s, e = int(m.group(1)), int(m.group(2))
    return 1 <= s <= e <= length


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--humanized", type=Path, required=True,
                    help="UpdatedTable.xlsx — uses Hoja 1 sheet.")
    ap.add_argument("--canonical", type=Path, required=True,
                    help="TF_completeIDs_with_ensembl_canonical.csv")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output hydrated_HumanizedTable.xlsx path.")
    ap.add_argument("--failures", type=Path, required=True,
                    help="TSV log of rows that still couldn't be hydrated.")
    args = ap.parse_args()

    print(f"Loading canonical map from {args.canonical}...")
    canon = _build_canonical_map(args.canonical)
    print(f"  {len(canon)} genes mapped to canonical ENSP + sequence")

    print(f"Loading HumanizedTable from {args.humanized}...")
    ht = pd.read_excel(args.humanized, sheet_name="Hoja 1")
    print(f"  {len(ht)} rows total")

    n_already = 0
    n_hydrated_human = 0
    n_human_no_gene = 0
    n_human_coord_overflow = 0
    n_non_human = 0
    failures: list[dict] = []

    for idx, row in ht.iterrows():
        existing_up = str(row.get("UNIPROT ID", "") or "").strip()
        existing_seq = str(row.get("UNIPROT SEQ", "") or "").strip()
        if existing_up and existing_seq and existing_seq.lower() != "nan":
            n_already += 1
            continue

        species = str(row.get("Species", "") or "")
        gene = str(row.get("Gene", "") or "").strip().upper()
        coords = row.get("DomainCoordinates")

        # Always try the canonical-ENSP lookup first — works for rows
        # tagged Human AND for the ~540 "Not specified." rows that turn
        # out to be human anyway (their Gene matches an HGNC symbol).
        # Only when this fails do we fall back to species-based routing
        # for stage 00b's UniProt API fetch.
        rec = canon.get(gene)

        if rec is None:
            # Gene name doesn't match any human canonical record. If the
            # row is also non-human, fall through to stage 00b. If it
            # *claims* human but we can't find it, log and skip.
            if _is_human(species):
                n_human_no_gene += 1
                failures.append({
                    "row_idx": idx,
                    "gene": gene,
                    "species": species,
                    "reason": "gene not in TF_completeIDs canonical map",
                })
            else:
                n_non_human += 1
                failures.append({
                    "row_idx": idx,
                    "gene": gene,
                    "species": species,
                    "reason": "NEEDS_UNIPROT_FETCH (non-human)",
                })
            continue

        if not _coords_in_range(coords, rec["length"]):
            # Gene resolves but coords overflow canonical length — likely
            # curated against a longer isoform.
            n_human_coord_overflow += 1
            failures.append({
                "row_idx": idx,
                "gene": gene,
                "species": species,
                "reason": (
                    f"DomainCoordinates '{coords}' don't fit canonical "
                    f"{rec['ensp']} ({rec['length']} aa) — likely curated "
                    f"against a longer isoform"
                ),
            })
            continue

        ht.at[idx, "UNIPROT ID"] = rec["uniprot"]
        ht.at[idx, "UNIPROT SEQ"] = rec["sequence"]
        ht.at[idx, "UNIPROT LENGTH"] = rec["length"]
        n_hydrated_human += 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.failures.parent.mkdir(parents=True, exist_ok=True)
    ht.to_excel(args.out, sheet_name="Hoja 1", index=False)
    pd.DataFrame(failures).to_csv(args.failures, sep="\t", index=False)

    print()
    print(f"Already hydrated:                   {n_already}")
    print(f"Hydrated from canonical map (human): {n_hydrated_human}")
    print(f"Skipped (non-human, needs API):     {n_non_human}")
    print(f"Skipped (gene not in canonical map): {n_human_no_gene}")
    print(f"Skipped (coord overflow):            {n_human_coord_overflow}")
    print(f"Output:   {args.out}")
    print(f"Failures: {args.failures}  ({len(failures)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
