"""Fetch gene names for all finches interaction partners from UniProtKB.

The partners are the broad STRING interactome (mostly non-TF), so local tables
only name ~31%. This pulls accession -> gene name for every distinct partner
UniProt via the UniProt ID-mapping REST API (one bulk job), writing a TSV that
convert_interactions.py consumes via --idmap.

Output: tfregdb_source/interaction/partner_names.tsv  (uniprot<TAB>name)
"""

from __future__ import annotations

import io
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SRC = REPO.parent / "tfregdb_source" / "interaction" / "_zip" / "epsilon_values_complete.csv"
OUT = REPO.parent / "tfregdb_source" / "interaction" / "partner_names.tsv"
API = "https://rest.uniprot.org"


def _post(url: str, data: dict) -> dict:
    body = urllib.parse.urlencode(data).encode()
    with urllib.request.urlopen(urllib.request.Request(url, data=body)) as r:
        import json
        return json.load(r)


def _get(url: str) -> str:
    with urllib.request.urlopen(url) as r:
        return r.read().decode()


def main() -> None:
    ids = sorted(set(pd.read_csv(SRC, low_memory=False).interactor_uniprot.astype(str)))
    print(f"[names] {len(ids)} distinct partner accessions; submitting idmapping job...")
    job = _post(f"{API}/idmapping/run",
                {"ids": ",".join(ids), "from": "UniProtKB_AC-ID", "to": "UniProtKB"})
    jid = job["jobId"]
    print(f"[names] jobId={jid}; polling...")

    import json
    for _ in range(120):  # up to ~10 min
        st = json.loads(_get(f"{API}/idmapping/status/{jid}"))
        if st.get("jobStatus") in (None, "FINISHED") or "results" in st:
            break
        if st.get("jobStatus") == "ERROR":
            raise SystemExit(f"[names] idmapping job errored: {st}")
        time.sleep(5)

    # Stream all results as TSV (accession + primary gene name).
    url = (f"{API}/idmapping/uniprotkb/results/stream/{jid}"
           f"?fields=accession,gene_primary,protein_name&format=tsv")
    tsv = _get(url)
    df = pd.read_csv(io.StringIO(tsv), sep="\t")
    gene_col = next((c for c in df.columns if "gene" in c.lower()), df.columns[1])
    acc_col = df.columns[0]
    # Fall back to protein name when no gene symbol, else leave blank.
    prot_col = next((c for c in df.columns if "protein" in c.lower()), None)
    rows = []
    for _, r in df.iterrows():
        nm = str(r[gene_col]) if pd.notna(r[gene_col]) and str(r[gene_col]).strip() else ""
        if not nm and prot_col and pd.notna(r[prot_col]):
            nm = str(r[prot_col])
        rows.append((str(r[acc_col]), nm))
    out = pd.DataFrame(rows, columns=["uniprot", "name"])
    out = out[out.name.str.len() > 0]
    out.to_csv(OUT, sep="\t", index=False)
    print(f"[names] resolved {len(out)}/{len(ids)} -> {OUT}")


if __name__ == "__main__":
    main()
