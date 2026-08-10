"""Stage 0b: hydrate non-human HumanizedTable rows via UniProt REST API.

Stage 00 fills in sequences for human rows by joining on TF_completeIDs.
That leaves ~1,400 mouse / rat / xenopus / etc. rows with empty UNIPROT
columns. Those rows are essential for the BLAST projection — they're the
inputs that get mapped back onto the human canonical via cross-species
sequence homology.

This stage:
  1. Reads the species column on every row that's still missing a sequence.
  2. Resolves Species name → NCBI taxonomy ID (small in-code map).
  3. For each (Gene, tax_id) pair, queries the UniProt REST API to find
     the reviewed accession (or first unreviewed if no Swiss-Prot entry
     exists) plus its canonical sequence.
  4. Hydrates the row's UNIPROT ID + UNIPROT SEQ + UNIPROT LENGTH columns.
  5. Caches resolutions in a JSON file so reruns are free.

The downstream BLAST step then projects each non-human sequence onto the
human canonical proteome — same flow as the original notebook.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests

# Common species labels in the HumanizedTable → NCBI taxonomy IDs. Add more
# entries as needed; rows whose Species string doesn't match a key here are
# logged as "unknown_species" and skipped.
SPECIES_TO_TAXON: dict[str, int] = {
    # Mammals
    "mouse": 10090,
    "mus musculus": 10090,
    "murine": 10090,
    "rat": 10116,
    "rattus": 10116,
    "bovine": 9913,
    "bos taurus": 9913,
    "cow": 9913,
    "pig": 9823,
    "sus scrofa": 9823,
    "rabbit": 9986,
    "oryctolagus": 9986,
    "dog": 9615,
    "canis": 9615,
    # Amphibians
    "xenopus": 8355,
    "xenopus laevis": 8355,
    "xenopus tropicalis": 8364,
    # Other vertebrates
    "chicken": 9031,
    "gallus": 9031,
    "zebrafish": 7955,
    "danio rerio": 7955,
    "atlantic cod": 8049,
    "gadus morhua": 8049,
    # Invertebrates
    "drosophila": 7227,
    "drosophila melanogaster": 7227,
    "c. elegans": 6239,
    "caenorhabditis elegans": 6239,
    # Yeasts
    "saccharomyces cerevisiae": 4932,
    "yeast": 4932,
    "schizosaccharomyces pombe": 4896,
    # Plants
    "arabidopsis thaliana": 3702,
    "arabidopsis": 3702,
    "oryza sativa": 39947,
    "rice": 39947,
    "wheat": 4565,
    "triticum aestivum": 4565,
    "maize": 4577,
    "zea mays": 4577,
    "capsicum annuum": 4072,
    "solanum tuberosum": 4113,
    "nicotiana tabacum": 4097,
    "medicago truncatula": 3880,
    "hordeum vulgare": 4513,
    "setaria italica": 4555,
    "ipomoea batatas": 4120,
    "thinopyrum intermedium": 87166,
}

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"


def _species_to_taxon(species: str) -> int | None:
    s = (species or "").strip().lower()
    if not s:
        return None
    # "Human" species rows shouldn't land here (stage 00 handles them), but
    # keep the fast-path return so we don't waste an API call.
    if s == "human" or "homo sapiens" in s:
        return 9606
    # Try exact match first, then any keyword from the map appearing inside.
    # Sort keys by length descending so "saccharomyces cerevisiae" wins over
    # "rice" when the substring is just "saccharomyces…".
    if s in SPECIES_TO_TAXON:
        return SPECIES_TO_TAXON[s]
    for key in sorted(SPECIES_TO_TAXON, key=len, reverse=True):
        if key in s:
            return SPECIES_TO_TAXON[key]
    return None


def _query_uniprot(gene: str, taxon: int, session: requests.Session) -> dict | None:
    """Return {'accession', 'sequence', 'length'} for the best hit or None.
    Prefers reviewed (Swiss-Prot) entries; falls back to unreviewed."""
    # `gene_exact` matches the gene name exactly (case-insensitive); reviewed
    # filter prefers Swiss-Prot. Fall back to unreviewed if no SP hit.
    for reviewed_filter in ("AND+reviewed:true", ""):
        q = f"gene_exact:{gene}+AND+organism_id:{taxon}{reviewed_filter}"
        url = (
            f"{UNIPROT_API}?query={q}"
            "&fields=accession,sequence,length&format=tsv&size=1"
        )
        try:
            r = session.get(url, timeout=20)
        except requests.RequestException:
            continue
        if r.status_code != 200:
            continue
        lines = r.text.strip().splitlines()
        if len(lines) < 2:
            continue  # header only
        cols = lines[1].split("\t")
        if len(cols) < 3:
            continue
        return {"accession": cols[0], "sequence": cols[1].strip().upper(),
                "length": int(cols[2]) if cols[2].isdigit() else len(cols[1])}
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hydrated", type=Path, required=True,
                    help="Output of 00_enrich_canonical (hydrated_HumanizedTable.xlsx).")
    ap.add_argument("--out", type=Path, required=True,
                    help="Final hydrated xlsx with animal rows filled in.")
    ap.add_argument("--cache", type=Path, required=True,
                    help="JSON cache of (gene, tax_id) → accession+seq.")
    ap.add_argument("--failures", type=Path, required=True)
    ap.add_argument("--sleep", type=float, default=0.1,
                    help="Sleep seconds between API calls (default 0.1).")
    ap.add_argument("--max-rows", type=int, default=None,
                    help="Cap total API rows for testing — leave unset for full run.")
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.cache.parent.mkdir(parents=True, exist_ok=True)
    args.failures.parent.mkdir(parents=True, exist_ok=True)

    cache: dict[str, dict] = {}
    if args.cache.exists():
        try:
            cache = json.loads(args.cache.read_text())
            print(f"Loaded {len(cache)} cached entries from {args.cache}")
        except json.JSONDecodeError:
            print(f"Warning: {args.cache} unreadable, starting fresh")

    print(f"Loading {args.hydrated}...")
    df = pd.read_excel(args.hydrated, sheet_name="Hoja 1")

    # Rows still missing both UNIPROT ID and UNIPROT SEQ. Stage 00 set
    # UNIPROT SEQ on rows it could hydrate; remaining are the animals plus
    # the handful of edge cases.
    def _empty(v) -> bool:
        s = str(v or "").strip().lower()
        return s in {"", "nan", "none"}

    todo_idx = [
        i for i, r in df.iterrows()
        if _empty(r.get("UNIPROT ID")) and _empty(r.get("UNIPROT SEQ"))
    ]
    if args.max_rows is not None:
        todo_idx = todo_idx[: args.max_rows]
    print(f"Rows to fetch: {len(todo_idx)}")

    session = requests.Session()
    session.headers.update({"User-Agent": "tfregdb-blast/0.1"})

    n_hit = n_cached = n_no_taxon = n_no_match = 0
    failures: list[dict] = []

    for n, idx in enumerate(todo_idx, 1):
        row = df.loc[idx]
        gene = str(row.get("Gene", "") or "").strip()
        species = str(row.get("Species", "") or "").strip()
        if not gene:
            failures.append({"row_idx": idx, "gene": gene, "species": species,
                             "reason": "empty gene"})
            continue
        taxon = _species_to_taxon(species)
        if taxon is None:
            n_no_taxon += 1
            failures.append({"row_idx": idx, "gene": gene, "species": species,
                             "reason": "unknown species (add to SPECIES_TO_TAXON)"})
            continue
        key = f"{gene.upper()}|{taxon}"
        if key in cache:
            hit = cache[key]
            n_cached += 1
        else:
            hit = _query_uniprot(gene, taxon, session)
            cache[key] = hit
            time.sleep(args.sleep)
            if n % 50 == 0:
                args.cache.write_text(json.dumps(cache))
        if not hit:
            n_no_match += 1
            failures.append({"row_idx": idx, "gene": gene, "species": species,
                             "tax_id": taxon, "reason": "no UniProt match"})
            continue
        df.at[idx, "UNIPROT ID"] = hit["accession"]
        df.at[idx, "UNIPROT SEQ"] = hit["sequence"]
        df.at[idx, "UNIPROT LENGTH"] = hit["length"]
        n_hit += 1
        if n % 100 == 0:
            print(f"  {n}/{len(todo_idx)} processed (hits={n_hit}, "
                  f"cached={n_cached}, no_taxon={n_no_taxon}, "
                  f"no_match={n_no_match})")

    args.cache.write_text(json.dumps(cache))
    df.to_excel(args.out, sheet_name="Hoja 1", index=False)
    pd.DataFrame(failures).to_csv(args.failures, sep="\t", index=False)

    print()
    print(f"API hits:               {n_hit}")
    print(f"Cache hits:             {n_cached}")
    print(f"Unknown species:        {n_no_taxon}")
    print(f"No UniProt match:       {n_no_match}")
    print(f"Output:   {args.out}")
    print(f"Cache:    {args.cache}  ({len(cache)} keys)")
    print(f"Failures: {args.failures}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
