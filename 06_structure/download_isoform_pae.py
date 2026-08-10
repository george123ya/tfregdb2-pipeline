#!/usr/bin/env python3
"""Download per-isoform PAE matrices from AlphaFold-DB into public/pae/{acc}.json
in the same {"length": N, "matrix": [[int,...],...]} shape the PAEViewer + the
canonical files already use. Resolves the current file version via the API
(AFDB's PAE filename drifts: …-predicted_aligned_error_v6.json), rounds floats
to ints to match the existing format and keep files small. Resumable.

Set of accessions = same distinct isoform accessions as download_isoform_pdbs.py.
Run:  python3 scripts/download_isoform_pae.py        (all)
      python3 scripts/download_isoform_pae.py ACC ..  (subset, for validation)
"""
import json, glob, os, sys, urllib.request, concurrent.futures as cf, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAE_DIR = os.path.join(ROOT, "public", "pae")
MANIFEST_GLOB = os.path.join(ROOT, "public", "mock", "structures", "*.json")

def collect_accessions():
    accs = set()
    for f in glob.glob(MANIFEST_GLOB):
        d = json.load(open(f))
        canon = d["uniprot"]
        for i in d.get("isoforms", []):
            a = i.get("uniprotIso")
            if i.get("afdb") and not i.get("canonical") and a and a != canon:
                accs.add(a)
    return sorted(accs)

def pae_url(acc):
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.load(r)
            return d[0].get("paeDocUrl") if d else None
        except urllib.error.HTTPError as e:
            if e.code in (404, 422):
                return None
            time.sleep(0.5)
        except Exception:
            time.sleep(0.5)
    return None

def fetch_one(acc):
    dest = os.path.join(PAE_DIR, f"{acc}.json")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return ("skip", acc)
    u = pae_url(acc)
    if not u:
        return ("missing", acc)
    for _ in range(3):
        try:
            with urllib.request.urlopen(u, timeout=90) as r:
                raw = json.load(r)
            obj = raw[0] if isinstance(raw, list) else raw
            pae = obj.get("predicted_aligned_error")
            if not pae:
                return ("badfmt", acc)
            matrix = [[int(round(v)) for v in row] for row in pae]
            out = {"length": len(matrix), "matrix": matrix}
            tmp = dest + ".part"
            with open(tmp, "w") as fh:
                json.dump(out, fh, separators=(",", ":"))
            os.replace(tmp, dest)
            return ("ok", acc)
        except Exception:
            time.sleep(0.7)
    return ("err", acc)

def main():
    os.makedirs(PAE_DIR, exist_ok=True)
    accs = collect_accessions()
    if len(sys.argv) > 1:
        wanted = set(sys.argv[1:])
        accs = [a for a in accs if a in wanted]
    print(f"distinct isoform accessions to ensure PAE for: {len(accs)}", flush=True)
    counts = {"ok": 0, "skip": 0, "missing": 0, "err": 0, "badfmt": 0}
    with cf.ThreadPoolExecutor(max_workers=10) as ex:
        for i, (status, acc) in enumerate(ex.map(fetch_one, accs), 1):
            counts[status] += 1
            if i % 200 == 0:
                print(f"  {i}/{len(accs)}  {counts}", flush=True)
    print(f"DONE {counts}", flush=True)

if __name__ == "__main__":
    main()
