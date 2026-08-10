"""Per-isoform PTM files from PTM_aggregated.xlsx.

dbPTM gives a 21-mer peptide + UniProt (NO isoform info). The upstream pipeline
(05_obtain_ptm.py) locates that exact peptide in EACH isoform's sequence and
recomputes the site position in that isoform's frame (`Real Position`) — i.e.
per-isoform PTMs are a *peptide-context presence* mapping (the modification site
exists in the isoform), NOT an isoform-specific observation. We surface them on
the PTM track with that framing.

The canonical PTMs are already baked into tfs/{symbol}.json by build_tf_records
(canonical translation bucket). This emits the SAME record shape for every
NON-canonical isoform, keyed by isoform NAME so the PTM track can fetch
mock/ptms/{isoformName}.json on isoform selection (matching the
conservation / adpred per-isoform convention).

Output: public/mock/ptms/{isoformName}.json
    [{ "pos", "type", "residue", "pmid"? }, ...]   (same as tfs ptms[])

Run from repo root (env with pandas+openpyxl):
    python scripts/extract_ptms_isoform.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SRC_PTMS = REPO.parent / "tfregdb_source" / "tables" / "PTM_aggregated.xlsx"
OUT_DIR = REPO / "public" / "mock" / "ptms"

# Same map build_tf_records uses; anything else is dropped.
PTM_TYPE_MAP = {
    "Phosphorylation": "Phospho",
    "Dephosphorylation": "Phospho",
    "Acetylation": "Acetyl",
    "Methylation": "Methyl",
    "Ubiquitination": "Ubiquitin",
    "Sumoylation": "Sumo",
    "Neddylation": "Sumo",
    "O-linked Glycosylation": "Glyco",
    "N-linked Glycosylation": "Glyco",
}


def main() -> None:
    if not SRC_PTMS.exists():
        raise SystemExit(f"missing {SRC_PTMS}")
    print(f"[ptm-iso] loading {SRC_PTMS} …", flush=True)
    df = pd.read_excel(SRC_PTMS)
    df["pmid"] = df["PMID"].astype(str).str.split(";").str[0].str.split(".").str[0].str.strip()
    df["ptype"] = df["PTM"].astype(str).map(PTM_TYPE_MAP).fillna("Other")
    df = df[df["ptype"] != "Other"]
    df["pos_i"] = pd.to_numeric(df["Real Position"], errors="coerce")
    df = df[df["pos_i"].notna()]
    df["pos_i"] = df["pos_i"].astype(int)
    df = df[df["pos_i"] >= 1]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for old in OUT_DIR.glob("*.json"):
        old.unlink()

    n_files = n_sites = 0
    for iso_name, grp in df.groupby("Isoform name"):
        iso_name = str(iso_name).strip()
        if not iso_name or iso_name.lower() == "nan":
            continue
        seen: set[tuple] = set()
        recs = []
        for _, r in grp.iterrows():
            pos = int(r["pos_i"])
            ptype = r["ptype"]
            key = (pos, ptype)
            if key in seen:
                continue
            seen.add(key)
            aa = str(r.get("Modified aa") or "X").strip()[:1] or "X"
            pmid = r["pmid"]
            rec = {"pos": pos, "type": ptype, "residue": f"{aa}{pos}"}
            if pmid and pmid != "nan":
                rec["pmid"] = pmid
            recs.append(rec)
        recs.sort(key=lambda p: (p["pos"], p["type"]))
        (OUT_DIR / f"{iso_name}.json").write_text(json.dumps(recs, separators=(",", ":")))
        n_files += 1
        n_sites += len(recs)

    print(f"[ptm-iso] wrote {n_files} isoform PTM files ({n_sites:,} sites) → {OUT_DIR}",
          flush=True)


if __name__ == "__main__":
    main()
