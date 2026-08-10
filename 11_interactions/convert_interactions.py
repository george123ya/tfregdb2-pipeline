"""Convert the finches epsilon SUMMARY into per-TF interaction JSON.

Summary-first build (the 28 GB matrix pkl is deferred to a workstation pass):
reads `epsilon_values_complete.csv` (one row per TF-IDR x partner with epsilon
stats) and writes one menu file per TF:

    public/interactions/{symbol}.json
    {
      "tf": "HOXA11",
      "idrs": [{
        "id": "242-270", "start": 242, "end": 270,
        "partners": [{                       # sorted most-attractive first
          "uniprot": "Q12778", "name": "FOXO1"|null, "ens": "ENSP..."|null,
          "mean": -0.907, "max": 2.706, "sum": -1202.4, "std": 1.146,
          "nValid": 1326, "shape": [91, 17]  # partner x TF-IDR (for later heatmap)
        }, ...]
      }]
    }

epsilon < 0 = attractive (favourable), > 0 = repulsive; partners are ranked by
mean epsilon ascending (strongest predicted attraction first).

Partner gene names are resolved best-effort from local sources (the broad
STRING interactome reaches well beyond our TF tables, so many fall back to the
UniProt accession). Pass --idmap FILE (csv/tsv/pkl of uniprot->name) — e.g. the
finches s2 `ensembl_to_uniprot` output — to complete names with no rebuild.

Usage:
    protein2genomic/.venv/bin/python scripts/convert_interactions.py [--idmap FILE]
"""

from __future__ import annotations

import argparse
import glob
import json
import pickle
import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SRC = REPO.parent / "tfregdb_source" / "interaction" / "_zip" / "epsilon_values_complete.csv"
OUT = REPO / "public" / "interactions"
TF_DIR = REPO / "public" / "mock" / "tfs"
OLD_CSV = REPO.parent / "tfregdb" / "data_pipeline" / "raw_files" / "TF_interactions_PDB_filtered.csv"
COMPLETE_IDS = REPO.parent / "tfregdb_source" / "tables" / "TF_completeIDs.xlsx"


def build_name_map(idmap: str | None) -> tuple[dict, dict]:
    name: dict[str, str] = {}
    ens: dict[str, str] = {}
    # 1) Partners that are TFs in our DB -> symbol (strongest source).
    for f in glob.glob(str(TF_DIR / "*.json")):
        d = json.loads(Path(f).read_text())
        if d.get("uniprot"):
            name.setdefault(str(d["uniprot"]), d["symbol"])
            if d.get("ensemblProtein"):
                ens.setdefault(str(d["uniprot"]), str(d["ensemblProtein"]))
    # 2) Old curated interaction CSV (uniprot -> partner name + ENSP).
    if OLD_CSV.exists():
        o = pd.read_csv(OLD_CSV)
        for u, n, e in zip(o.UniProtID.astype(str), o.Partner_Name.astype(str),
                           o.Partner_EnsemblID.astype(str)):
            name.setdefault(u, n)
            ens.setdefault(u, e)
    # 3) TF_completeIDs.xlsx (UniProt_ID -> HGNC name).
    if COMPLETE_IDS.exists():
        try:
            tc = pd.read_excel(COMPLETE_IDS)
            cols = {c.lower(): c for c in tc.columns}
            uc = next((cols[c] for c in cols if "uniprot" in c), None)
            nc = next((cols[c] for c in cols if "hgnc" in c or c in ("gene name", "gene_name")), None)
            if uc and nc:
                for u, n in zip(tc[uc].astype(str), tc[nc].astype(str)):
                    if u and u != "nan":
                        name.setdefault(u, n)
        except Exception as e:  # noqa: BLE001
            print(f"[interactions] xlsx name map skipped: {e}")
    # 4) Optional authoritative idmap (uniprot -> name).
    if idmap:
        p = Path(idmap)
        if p.suffix == ".pkl":
            m = pickle.load(open(p, "rb"))
            for u, n in (m.items() if isinstance(m, dict) else []):
                name[str(u)] = str(n)
        else:
            sep = "\t" if p.suffix in (".tsv", ".txt") else ","
            m = pd.read_csv(p, sep=sep)
            uc = next((c for c in m.columns if "uniprot" in c.lower()), m.columns[0])
            nc = next((c for c in m.columns if "name" in c.lower() or "gene" in c.lower()), m.columns[-1])
            for u, n in zip(m[uc].astype(str), m[nc].astype(str)):
                name[u] = n
    return name, ens


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--idmap", default=None, help="uniprot->name csv/tsv/pkl")
    args = ap.parse_args()

    df = pd.read_csv(SRC, low_memory=False)
    df = df[df["valid"] == True]  # noqa: E712
    df["tf"] = df.idr_name.str.rsplit("_", n=1).str[0]
    name, ens = build_name_map(args.idmap)

    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob("*.json"):
        old.unlink()

    site_tfs = {Path(f).stem for f in glob.glob(str(TF_DIR / "*.json"))}
    written = no_site = 0
    named = total_partners = 0
    shape_re = re.compile(r"\((\d+),\s*(\d+)\)")

    for tf, g in df.groupby("tf"):
        if tf not in site_tfs:
            no_site += 1
            continue
        idrs = []
        for idr_name, gi in g.groupby("idr_name"):
            rng = idr_name.rsplit("_", 1)[1]
            try:
                s, e = (int(x) for x in rng.split("-"))
            except ValueError:
                continue
            partners = []
            for _, r in gi.sort_values("epsilon_mean").iterrows():
                u = str(r.interactor_uniprot)
                nm = name.get(u)
                total_partners += 1
                if nm:
                    named += 1
                m = shape_re.search(str(r.matrix_shape))
                partners.append({
                    "uniprot": u,
                    "name": nm,
                    "ens": ens.get(u),
                    "mean": round(float(r.epsilon_mean), 3),
                    "max": round(float(r.epsilon_max), 3),
                    "sum": round(float(r.epsilon_sum), 1),
                    "std": round(float(r.epsilon_std), 3),
                    "nValid": int(r.n_valid_cells),
                    "shape": [int(m.group(1)), int(m.group(2))] if m else None,
                })
            idrs.append({"id": rng, "start": s, "end": e, "partners": partners})
        (OUT / f"{tf}.json").write_text(
            json.dumps({"tf": tf, "idrs": idrs}, separators=(",", ":")))
        written += 1

    print(f"[interactions] wrote {written} TF files to {OUT}")
    print(f"[interactions] {no_site} finches TFs had no matching site TF (skipped)")
    print(f"[interactions] partner names resolved: {named}/{total_partners} "
          f"({100 * named // max(total_partners, 1)}%) — rest show UniProt ID")


if __name__ == "__main__":
    main()
