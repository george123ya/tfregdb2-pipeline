"""Stage 3: project BLAST hits onto the canonical human Ensembl protein and
emit `Unified_Valid_Domains.xlsx`.

For each query (one curated effector domain row from HumanizedTable or
TFRegDB1):

  1. Find the best BLAST hit (max bitscore — already constrained to the
     single best human hit by `-max_target_seqs 1` upstream).
  2. Decide whether the hit clears the >95% identity threshold; if not, log
     it to `blast_failures.tsv` and keep going.
  3. Use the canonical CSV to look up the canonical Ensembl protein for
     the hit's `sseqid` (which is the UniProt of the matched human gene).
  4. Compute the canonical-frame `[sstart, send]` coordinates — those become
     `Final_Human_DomainCoords`.
  5. Write back the source row + the new Final_* columns.

The output is one combined `Unified_Valid_Domains.xlsx` with a stable column
schema downstream readers (build_tf_records.py) already consume.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]

IDENTITY_THRESHOLD = 95.0


def _build_canonical_dicts(canonical_csv: Path) -> tuple[dict, dict, dict]:
    """Return (gene→ENSP, uniprot→ENSP, ENSP→seq) maps. Within a group, the
    row tagged `Ensembl canonical == "Canonical"` wins; otherwise the first
    APPRIS PRINCIPAL row; otherwise the longest."""
    df = pd.read_csv(canonical_csv)
    df["Translation_Clean"] = df["Translation ID"].astype(str).str.split(".").str[0]
    df["Gene_Upper"] = df["Gene name (HGNC)"].astype(str).str.strip().str.upper()
    df["UniProt_Clean"] = df["UniProt_ID"].astype(str).str.split("-").str[0].str.strip()

    def pick(group: pd.DataFrame) -> str:
        canon = group[group["Ensembl canonical"].astype(str).str.strip() == "Canonical"]
        if len(canon):
            return canon["Translation_Clean"].iloc[0]
        p1 = group[group["APPRIS Annotation"].astype(str).str.contains("PRINCIPAL:1", na=False)]
        if len(p1):
            return p1["Translation_Clean"].iloc[0]
        return group.sort_values("Protein length", ascending=False)["Translation_Clean"].iloc[0]

    gene_to_ensp: dict[str, str] = {}
    up_to_ensp: dict[str, str] = {}
    for g, grp in df.dropna(subset=["Translation_Clean"]).groupby("Gene_Upper"):
        if g and g != "NAN":
            gene_to_ensp[g] = pick(grp)
    for u, grp in df.dropna(subset=["Translation_Clean"]).groupby("UniProt_Clean"):
        if u and u != "NAN":
            up_to_ensp[u] = pick(grp)
    ensp_to_seq = dict(zip(df["Translation_Clean"], df["Sequence aa"]))
    return gene_to_ensp, up_to_ensp, ensp_to_seq


def _parse_query_id(qseqid: str) -> dict:
    """`HumanizedTable_row42|TP53|P04637` → {source, row_idx, gene, uniprot}."""
    parts = qseqid.split("|")
    head = parts[0] if parts else qseqid
    src, _, row = head.partition("_row")
    return {
        "source": src,
        "row_idx": int(row) if row.isdigit() else None,
        "gene": parts[1] if len(parts) > 1 else "",
        "uniprot": parts[2] if len(parts) > 2 else "",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast", type=Path, required=True)
    ap.add_argument("--canonical", type=Path, required=True)
    ap.add_argument("--humanized", type=Path, required=True)
    ap.add_argument("--tfregdb1", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--failures", type=Path, required=True)
    args = ap.parse_args()

    blast = pd.read_csv(args.blast, sep="\t", names=BLAST_COLS, header=None)
    # Keep best hit per query (highest bitscore).
    blast = blast.sort_values(["qseqid", "bitscore"], ascending=[True, False])
    best_hit = blast.drop_duplicates("qseqid", keep="first").set_index("qseqid")

    gene_to_ensp, up_to_ensp, ensp_to_seq = _build_canonical_dicts(args.canonical)

    # Re-read source tables so we can stitch back every original column.
    src_rows: list[dict] = []
    for kind, path, sheet in (
        ("HumanizedTable", args.humanized, "Hoja 1"),
        ("TFRegDB1", args.tfregdb1, "Hoja 2"),
    ):
        df = pd.read_excel(path, sheet_name=sheet)
        df["__source__"] = kind
        df["__row_idx__"] = df.index
        src_rows.extend(df.to_dict("records"))

    out_rows: list[dict] = []
    failures: list[dict] = []

    for r in src_rows:
        qid = f"{r['__source__']}_row{r['__row_idx__']}"
        hit = best_hit[best_hit.index.str.startswith(qid + "|")]
        if hit.empty:
            failures.append({**r, "reason": "no BLAST hit"})
            continue
        h = hit.iloc[0]
        pident = float(h["pident"])
        sseqid = str(h["sseqid"])
        # sseqid format depends on how the proteome FASTA was built. The
        # human UniProt FASTA convention is `sp|<accession>|<entry-name>`,
        # so we strip to the bare accession.
        up_acc = sseqid.split("|")[1] if "|" in sseqid else sseqid
        canon_ensp = up_to_ensp.get(up_acc)
        canon_seq = ensp_to_seq.get(canon_ensp, "") if canon_ensp else ""

        # >95% identity = trust the canonical projection; otherwise log.
        if pident < IDENTITY_THRESHOLD:
            failures.append({
                **r,
                "reason": f"identity {pident:.1f}% < {IDENTITY_THRESHOLD}%",
                "BLAST_UNIPROT_ID": up_acc,
                "BLAST_Identity_Pct": pident,
            })
            continue

        s_start = int(h["sstart"])
        s_end = int(h["send"])
        coords = f"{s_start}-{s_end}"
        # Legacy-pipeline column aliases — downstream build_tf_records.py reads
        # `Original_Gene`, `Source_Database`, and `PUBMED_ID`. We map them from
        # whichever source table this row came from so the schema is stable.
        src = str(r.get("__source__") or "")
        orig_gene = (
            str(r.get("Gene") or "").strip() if src == "HumanizedTable"
            else str(r.get("TF name") or "").strip()
        )
        pubmed = str(r.get("PMID") or "").strip()
        out_rows.append({
            **r,
            "Original_Gene": orig_gene,
            "Source_Database": src,
            "PUBMED_ID": pubmed,
            "BLAST_UNIPROT_ID": up_acc,
            "BLAST_DomainCoordinates": coords,
            "BLAST_Identity_Pct": pident,
            "Final_Human_UNIPROT_ID": up_acc,
            "Final_Human_DomainCoords": coords,
            "Final_Identity": pident,
            "Mapping_Status": "Successful BLAST Hit (>95%)",
            "Final_Human_ENSP_ID": canon_ensp,
            "Final_Human_Domain_Sequence": (
                canon_seq[s_start - 1 : s_end] if canon_seq else ""
            ),
        })

    df_out = pd.DataFrame(out_rows)
    df_fail = pd.DataFrame(failures)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_excel(args.out, index=False)
    df_fail.to_csv(args.failures, sep="\t", index=False)

    print(f"Projected: {len(df_out)} rows → {args.out}")
    print(f"Failures:  {len(df_fail)} rows → {args.failures}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
