"""Stage 1: build queries.fasta from the curated input tables.

One FASTA record per row of HumanizedTable (Hoja 1) and TFRegDB1 (Hoja 2)
in `UpdatedTable.xlsx`. Header encodes (source, row index, gene, original
UniProt) so the projector can rejoin BLAST hits to their source row.

Header format: `>{source}_row{N}|{gene}|{uniprot}`

For HumanizedTable, the domain sequence is sliced live from the row's
`UNIPROT SEQ` using the `DomainCoordinates` column (e.g. "555-701"). For
TFRegDB1 the `Sequence` column carries the domain sequence directly. Rows
that can't yield a sequence ≥ 5 aa are skipped and counted in the summary.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Some xlsx cells use en/em dashes — normalise to plain hyphen before parsing.
DASHES = ("–", "—", "−")


def _slice_domain(uniprot_seq: str, coords: str) -> str:
    """Slice the domain out of `uniprot_seq` according to `coords`
    ("123-456" / "1-100" / "-50:" / "100-"). Returns "" if anything is off."""
    if not uniprot_seq or not coords:
        return ""
    seq = str(uniprot_seq).strip().upper()
    raw = str(coords).strip()
    for d in DASHES:
        raw = raw.replace(d, "-")
    if not seq or seq.lower() == "nan" or not raw or raw.lower() == "nan":
        return ""
    # Accept "N-M", "N-", "-N", ":N-M", etc. Be permissive.
    try:
        # Strip a leading ":" if present (rare formatting).
        raw = raw.lstrip(":")
        a, _, b = raw.partition("-")
        a, b = a.strip(), b.strip()
        if not a and b:  # "-50" → last 50 aa
            n = int(b)
            return seq[-n:] if 0 < n <= len(seq) else ""
        if a and not b:  # "100-" → from 100 to end
            start = int(a)
            return seq[max(0, start - 1):] if start >= 1 else ""
        start, end = int(a), int(b)
        if start < 1 or end < start or end > len(seq):
            return ""
        return seq[start - 1:end]
    except ValueError:
        return ""


def _emit(rows, out_fh, source, hdr) -> tuple[int, int]:
    n_w = n_s = 0
    for idx, row in rows.iterrows():
        seq = hdr["seq_fn"](row).strip().upper() if "seq_fn" in hdr else \
            str(row.get(hdr["seq_col"], "")).strip().upper()
        if not seq or seq.lower() in {"nan", "none"} or len(seq) < 5:
            n_s += 1
            continue
        gene = str(row.get(hdr["gene_col"], "")).strip() or "UNK"
        up_raw = str(row.get(hdr["uniprot_col"], "")).strip()
        up = up_raw.split("-")[0] if up_raw else "UNK"
        out_fh.write(f">{source}_row{idx}|{gene}|{up}\n{seq}\n")
        n_w += 1
    return n_w, n_s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--humanized", type=Path, required=True,
                    help="Path to UpdatedTable.xlsx (Hoja 1 sheet).")
    ap.add_argument("--tfregdb1", type=Path, required=True,
                    help="Path to UpdatedTable.xlsx (Hoja 2 sheet). May be the same file.")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    total_w = total_s = 0
    with args.out.open("w") as fh:
        # ── HumanizedTable (Hoja 1) ──────────────────────────────────────
        ht = pd.read_excel(args.humanized, sheet_name="Hoja 1")
        w, s = _emit(
            ht, fh, source="HumanizedTable",
            hdr={
                "gene_col": "Gene",
                "uniprot_col": "UNIPROT ID",
                "seq_fn": lambda r: _slice_domain(r.get("UNIPROT SEQ"),
                                                  r.get("DomainCoordinates")),
            },
        )
        print(f"HumanizedTable: {w} written / {s} skipped (no sequence)")
        total_w += w; total_s += s

        # ── TFRegDB1 (Hoja 2) — the older effector-domain list ───────────
        tf1 = pd.read_excel(args.tfregdb1, sheet_name="Hoja 2")
        w, s = _emit(
            tf1, fh, source="TFRegDB1",
            hdr={
                "gene_col": "TF name",
                "uniprot_col": "Uniprot ID",
                "seq_col": "Sequence",
            },
        )
        print(f"TFRegDB1:       {w} written / {s} skipped (no sequence)")
        total_w += w; total_s += s

    print(f"\nTotal: {total_w} queries → {args.out}  (skipped {total_s})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
