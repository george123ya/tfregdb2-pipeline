#!/usr/bin/env python3
"""Download the distinct isoform AlphaFold-DB models referenced by the
structure manifests into public/pdb/AF-{acc}-F1.pdb — self-hosted alongside
the canonical models so the viewer never depends on AFDB's drifting file
version suffix (…-model_v4 → _v6 → …).

Set of accessions: every `uniprotIso` in public/mock/structures/*.json that is
afdb-available, non-canonical, and != the canonical accession (deduped). Files
already present (canonical downloads or a prior run) are skipped.

Run:  python3 scripts/download_isoform_pdbs.py
"""
import json, glob, os, sys, urllib.request, concurrent.futures as cf, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDB_DIR = os.path.join(ROOT, "public", "pdb")
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

def afdb_pdb_url(acc):
    """Resolve the current model file URL from the AFDB API (version-proof)."""
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.load(r)
            return d[0].get("pdbUrl") if d else None
        except urllib.error.HTTPError as e:
            if e.code in (404, 422):
                return None
            time.sleep(0.5)
        except Exception:
            time.sleep(0.5)
    return None

def fetch_one(acc):
    dest = os.path.join(PDB_DIR, f"AF-{acc}-F1.pdb")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return ("skip", acc)
    pdb_url = afdb_pdb_url(acc)
    if not pdb_url:
        return ("missing", acc)
    for _ in range(3):
        try:
            with urllib.request.urlopen(pdb_url, timeout=60) as r:
                data = r.read()
            tmp = dest + ".part"
            with open(tmp, "wb") as fh:
                fh.write(data)
            os.replace(tmp, dest)
            return ("ok", acc)
        except Exception:
            time.sleep(0.7)
    return ("err", acc)

def main():
    os.makedirs(PDB_DIR, exist_ok=True)
    accs = collect_accessions()
    # Optional CLI subset: pass accessions to fetch just those (validation).
    if len(sys.argv) > 1:
        wanted = set(sys.argv[1:])
        accs = [a for a in accs if a in wanted]
    print(f"distinct isoform accessions to ensure: {len(accs)}", flush=True)
    counts = {"ok": 0, "skip": 0, "missing": 0, "err": 0}
    missing, errs = [], []
    with cf.ThreadPoolExecutor(max_workers=12) as ex:
        for i, (status, acc) in enumerate(ex.map(fetch_one, accs), 1):
            counts[status] += 1
            if status == "missing":
                missing.append(acc)
            elif status == "err":
                errs.append(acc)
            if i % 200 == 0:
                print(f"  {i}/{len(accs)}  {counts}", flush=True)
    print(f"DONE {counts}", flush=True)
    if missing:
        open(os.path.join(ROOT, "scripts", "isoform_pdb_missing.txt"), "w").write("\n".join(missing))
        print(f"  {len(missing)} not in AFDB → scripts/isoform_pdb_missing.txt", flush=True)
    if errs:
        print(f"  {len(errs)} download errors (re-run to retry): {errs[:10]}", flush=True)

if __name__ == "__main__":
    main()
