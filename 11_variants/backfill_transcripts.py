"""Backfill `transcripts[]` on every per-TF JSON from TF_completeIDs.

Pulls (UniProt, Transcript ID, Translation ID, Protein length, Isoform name)
per row and attaches a Transcript[] cross-link (ENST → iso name) to each
public/mock/tfs/{SYMBOL}.json. Lets the variant projector match COSMIC's
per-transcript rows to isoforms without re-running the full
build_tf_records pipeline.

Idempotent — re-running just overwrites the field.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pandas as pd


IDS_CSV = Path(os.environ.get(
    "COMPLETE_IDS_CSV",
    os.path.expanduser(
        "~/Desktop/tfregdb/data_pipeline/blast_table/"
        "TF_completeIDs_with_ensembl_canonical.csv"
    ),
))
TFS_DIR = Path("public/mock/tfs")


def _clean(s: object) -> str:
    return str(s or "").split(".")[0].strip()


def main() -> int:
    ids = pd.read_csv(IDS_CSV, low_memory=False)
    ids["UniProt_Clean"] = ids["UniProt_ID"].astype(str).str.strip()
    ids["ENST"] = ids["Transcript ID"].astype(str).str.split(".").str[0].str.strip()
    ids["ENSP"] = ids["Translation ID"].astype(str).str.split(".").str[0].str.strip()

    n_ok = n_skip = 0
    for fp in sorted(TFS_DIR.glob("*.json")):
        tf = json.loads(fp.read_text())
        upro = str(tf.get("uniprot") or "").strip()
        if not upro:
            n_skip += 1
            continue
        # Match by UniProt accession (drops the -2 / -3 isoform suffix).
        upro_base = upro.split("-", 1)[0]
        rows = ids[ids["UniProt_Clean"].str.startswith(upro_base, na=False)]
        if rows.empty:
            # Fall back to gene name.
            sym = str(tf.get("symbol") or "").upper()
            rows = ids[ids["Gene name (HGNC)"].astype(str).str.upper() == sym]
        if rows.empty:
            n_skip += 1
            continue

        # The CSV's `Ensembl canonical` column carries the literal string
        # "Canonical" on the one canonical row + NaN elsewhere — not a
        # boolean. We also accept MANE_Select as a tiebreaker, then
        # APPRIS PRINCIPAL:1, then row order. Picking the wrong row here
        # silently mis-labels the canonical transcript for the whole TF.
        canon_mask = rows["Ensembl canonical"].astype(str).str.strip().str.lower() == "canonical"
        if canon_mask.any():
            canon_row = rows.loc[canon_mask].iloc[0]
        else:
            mane_mask = rows["MANE"].astype(str).str.contains("MANE_Select", case=False, na=False)
            if mane_mask.any():
                canon_row = rows.loc[mane_mask].iloc[0]
            else:
                appris_mask = rows["APPRIS Annotation"].astype(str).str.startswith("PRINCIPAL:1", na=False)
                canon_row = rows.loc[appris_mask].iloc[0] if appris_mask.any() else rows.iloc[0]
        canon_enst = _clean(canon_row["Transcript ID"])
        canon_ensp_str = _clean(canon_row["Translation ID"])

        # Map iso name → ENST by joining on Translation ID against the
        # existing isoforms[].proteinId. Falls back to APPRIS / row order
        # for isoforms whose Translation ID is missing.
        ensp_to_enst: dict[str, str] = {}
        for _, r in rows.iterrows():
            ensp = _clean(r["Translation ID"])
            enst = _clean(r["Transcript ID"])
            if ensp and enst:
                ensp_to_enst[ensp] = enst

        transcripts: list[dict] = []
        seen_enst: set[str] = set()
        for iso in tf.get("isoforms") or []:
            iso_name = iso.get("name")
            iso_ensp = iso.get("proteinId") or ""
            iso_len = int(iso.get("length") or 0)
            enst = ensp_to_enst.get(iso_ensp)
            if not enst or enst in seen_enst:
                continue
            transcripts.append({
                "name": iso_name,
                "ensemblId": enst,
                "length": iso_len * 3 + 3,
                "exons": [{"start": 1, "end": iso_len * 3 + 3}],
                "proteinIsoform": iso_name,
                "canonical": enst == canon_enst,
            })
            seen_enst.add(enst)

        if not transcripts:
            n_skip += 1
            continue

        tf["transcripts"] = transcripts
        fp.write_text(json.dumps(tf, separators=(",", ":")))
        n_ok += 1

    print(f"[backfill_transcripts] updated {n_ok} TF JSONs · skipped {n_skip}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
