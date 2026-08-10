"""Generate the bulk download files for the Downloads page.

Lightweight, flat tables + a per-TF JSON bundle — the things worth pre-bundling.
Big per-residue / structural data (AlphaFold PDB ~1 GB, PAE ~2 GB, conservation)
is NOT bundled here: it already lives on R2 as individual files and is reached
via the API (/structures, /biophysics, /conservation) — bundling it would blow
the 10 GB R2 free tier.

Outputs → public/downloads/ (then synced to R2 v1/downloads/ by upload_r2.sh):
  effector_domains.tsv        all curated AD/RD/Bif + DBD domains
  ptms.tsv                    all PTM sites (+ which curated domain they hit)
  variants_in_domains.csv.gz  in-domain ClinVar + COSMIC variants
  tfs_json.tgz                every per-TF JSON profile
  README.txt
"""

from __future__ import annotations

import csv
import gzip
import json
import shutil
import subprocess
from pathlib import Path

OUT = Path("public/downloads")
API = Path("public/mock/api")
TFS = Path("public/mock/tfs")
VAR_CSV = Path("data/variants/api_in_domains.csv")


def write_tsv(path: Path, header: list[str], rows):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    # 1. Effector + DNA-binding domains.
    doms = json.loads((API / "domains.json").read_text())
    write_tsv(OUT / "effector_domains.tsv",
              ["domain_id", "symbol", "uniprot", "ensp", "type", "subtype", "family",
               "tf_family", "start", "end", "length", "confidence", "pmid", "year"],
              ([d["id"], d["symbol"], d["uniprot"], d.get("ensp") or "", d["type"], d["subtype"],
                d.get("family") or "", d.get("tf_family") or "", d["start"], d["end"], d["length"],
                d.get("confidence") or "", d.get("pmid") or "", d.get("year") or ""] for d in doms))
    print(f"  effector_domains.tsv · {len(doms):,} domains")

    # 2. PTMs.
    ptms = json.loads((API / "ptms.json").read_text())
    write_tsv(OUT / "ptms.tsv",
              ["gene", "uniprot", "pos", "residue", "type", "pmid", "in_domain", "domain_types"],
              ([p["gene"], p["uniprot"], p["pos"], p.get("residue") or "", p.get("type") or "",
                p.get("pmid") or "", "yes" if p.get("domains") else "no",
                "|".join(sorted({d["type"] for d in p.get("domains") or []}))] for p in ptms))
    print(f"  ptms.tsv · {len(ptms):,} sites")

    # 3. In-domain variants (ClinVar + COSMIC) — gzip the existing CSV.
    if VAR_CSV.exists():
        with VAR_CSV.open("rb") as src, gzip.open(OUT / "variants_in_domains.csv.gz", "wb") as dst:
            shutil.copyfileobj(src, dst)
        print(f"  variants_in_domains.csv.gz · {(OUT / 'variants_in_domains.csv.gz').stat().st_size/1e6:.1f} MB")

    # 4. Per-TF JSON bundle.
    subprocess.run(["tar", "czf", str(OUT / "tfs_json.tgz"), "-C", str(TFS.parent), TFS.name], check=True)
    print(f"  tfs_json.tgz · {(OUT / 'tfs_json.tgz').stat().st_size/1e6:.1f} MB")

    (OUT / "README.txt").write_text(
        "TFRegDB2 bulk downloads (CC-BY-4.0)\n"
        "===================================\n\n"
        "effector_domains.tsv        Curated AD/RD/Bifunctional effector + CIS-BP 3.0 DNA-binding domains.\n"
        "                            Columns: domain_id, symbol, uniprot, ensp, type, subtype, family,\n"
        "                            tf_family, start, end, length, confidence, pmid, year.\n"
        "ptms.tsv                    PTM sites; in_domain flags whether the site lies in a curated domain.\n"
        "variants_in_domains.csv.gz  ClinVar + COSMIC variants that fall inside a curated domain\n"
        "                            (gnomAD excluded — see the TF pages / API). Columns: gene, uniprot,\n"
        "                            pos, aa_ref, aa_alt, consequence, source, clin_sig, af, count,\n"
        "                            domain_id, domain_type, domain_subtype, family, tf_family.\n"
        "tfs_json.tgz                Every per-TF JSON profile (full payload behind each /tf/{SYMBOL} page).\n\n"
        "For custom slices, query the API: https://tfregdb2.kamluwantan.com/api/v1/docs\n"
        "AlphaFold PDB / PAE / per-residue conservation are large; fetch individually via the API\n"
        "(/api/v1/structures, /biophysics) rather than as bundles.\n\n"
        "Cite v1: Soto et al., Compendium of human TF effector domains, Mol Cell 2021,\n"
        "doi:10.1016/j.molcel.2021.11.007\n"
    )
    print(f"  → {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
