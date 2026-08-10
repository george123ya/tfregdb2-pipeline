"""Inject HGNC ID + CCDS IDs into each per-TF JSON, from the HGNC complete set.

The HGNC "complete set" TSV (genenames.org) maps every approved gene symbol to
its stable HGNC accession and its Consensus CDS (CCDS) IDs — so one file gives
both. We match by symbol (the TF JSON's `symbol`) and add:
    "hgnc": "HGNC:9864"
    "ccds": ["CCDS11537", "CCDS45693"]   # [] when none
to public/mock/tfs/{SYMBOL}.json. The TF page renders these as external-link
chips under the name (HGNC → genenames.org, CCDS → NCBI CCDS browser).

Idempotent: re-running just overwrites the two fields.

Download first (cached at data/hgnc/hgnc_complete_set.txt):
    curl -sL https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt -o data/hgnc/hgnc_complete_set.txt
Then run from repo root:
    python3 scripts/add_hgnc_ccds.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TF_DIR = REPO / "public" / "mock" / "tfs"
HGNC_TSV = REPO / "data" / "hgnc" / "hgnc_complete_set.txt"


def load_lookup() -> dict[str, dict]:
    """symbol(upper) -> {hgnc, ccds[]}, incl. previous symbols/aliases as fallbacks."""
    primary: dict[str, dict] = {}
    alias: dict[str, dict] = {}
    with HGNC_TSV.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            hgnc = (row.get("hgnc_id") or "").strip()
            sym = (row.get("symbol") or "").strip()
            if not hgnc or not sym:
                continue
            ccds = [c.strip() for c in (row.get("ccds_id") or "").split("|") if c.strip()]
            entry = {"hgnc": hgnc, "ccds": ccds}
            primary[sym.upper()] = entry
            for field in ("prev_symbol", "alias_symbol"):
                for a in (row.get(field) or "").split("|"):
                    a = a.strip().upper()
                    if a and a not in alias:
                        alias[a] = entry
    # primary wins over alias
    return {**alias, **primary}


def main() -> int:
    if not HGNC_TSV.exists():
        sys.exit(f"missing {HGNC_TSV} — download the HGNC complete set first")
    lut = load_lookup()
    print(f"[hgnc] {len(lut)} symbol/alias keys loaded", flush=True)

    n_ok = n_ccds = n_miss = 0
    missing = []
    for fp in sorted(TF_DIR.glob("*.json")):
        d = json.loads(fp.read_text())
        sym = (d.get("symbol") or fp.stem).upper()
        ent = lut.get(sym)
        if not ent:
            n_miss += 1
            missing.append(sym)
            d["hgnc"] = None
            d["ccds"] = []
        else:
            d["hgnc"] = ent["hgnc"]
            d["ccds"] = ent["ccds"]
            n_ok += 1
            if ent["ccds"]:
                n_ccds += 1
        fp.write_text(json.dumps(d, separators=(",", ":")))

    print(f"[hgnc] matched {n_ok} TFs ({n_ccds} with CCDS) · {n_miss} unmatched", flush=True)
    if missing:
        print(f"[hgnc] unmatched: {', '.join(missing[:15])}"
              f"{' …' if len(missing) > 15 else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
