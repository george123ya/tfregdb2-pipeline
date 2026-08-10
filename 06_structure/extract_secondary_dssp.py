#!/usr/bin/env python3
"""Isoform secondary structure via DSSP — the SAME program as the canonical set.
mkdssp 4.x only reads mmCIF (it chokes on AlphaFold's PDB header) and gemmi's
re-serialised CIF loses fields mkdssp needs, but mkdssp parses AFDB's own CIF
perfectly. So per isoform we fetch the AFDB mmCIF, run mkdssp, parse its classic
SS column, and write public/mock/secondary/{acc}.json = {"ss": "..."} — same
shape/3-state mapping as scripts/extract_secondary.py. CIFs are streamed to a
temp file and deleted; only the small SS JSON is kept. Resumable (skips existing).

Run with the conda base python (urllib only):
    python scripts/extract_secondary_dssp.py
"""
import json, glob, os, subprocess, sys, tempfile, urllib.request, time
import concurrent.futures as cf

MKDSSP = os.environ.get("MKDSSP_BIN", "mkdssp")
OUT_DIR = "public/mock/secondary"
SS_MAP = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E",
          "T": "C", "S": "C", " ": "C", "-": "C", "P": "C", "": "C"}


def distinct_iso_accessions():
    accs = set()
    for f in glob.glob("public/mock/structures/*.json"):
        d = json.load(open(f))
        canon = d["uniprot"]
        for i in d.get("isoforms", []):
            a = i.get("uniprotIso")
            if i.get("afdb") and not i.get("canonical") and a and a != canon:
                accs.add(a)
    return sorted(accs)


def cif_url(acc):
    for _ in range(3):
        try:
            with urllib.request.urlopen(
                f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=30
            ) as r:
                d = json.load(r)
            return d[0].get("cifUrl") if d else None
        except urllib.error.HTTPError as e:
            if e.code in (404, 422):
                return None
            time.sleep(0.5)
        except Exception:
            time.sleep(0.5)
    return None


def parse_dssp(text):
    lines = text.splitlines()
    start = next((i for i, ln in enumerate(lines) if ln.startswith("  #  RESIDUE")), None)
    if start is None:
        return None
    out = [SS_MAP.get(ln[16], "C") for ln in lines[start + 1:]
           if len(ln) > 16 and ln[13] != "!"]
    return "".join(out) if out else None


def one(acc):
    dest = os.path.join(OUT_DIR, f"{acc}.json")
    if os.path.exists(dest):
        return "skip"
    url = cif_url(acc)
    if not url:
        return "missing"
    try:
        with tempfile.TemporaryDirectory() as t:
            cif = os.path.join(t, "m.cif")
            with urllib.request.urlopen(url, timeout=90) as r:
                data = r.read()
            with open(cif, "wb") as fh:
                fh.write(data)
            dssp = os.path.join(t, "m.dssp")
            r = subprocess.run([MKDSSP, "--output-format", "dssp", cif, dssp],
                               capture_output=True, text=True, timeout=180)
            if r.returncode != 0 or not os.path.exists(dssp):
                return "fail"
            ss = parse_dssp(open(dssp).read())
            if not ss:
                return "fail"
            with open(dest, "w") as fh:
                json.dump({"ss": ss}, fh, separators=(",", ":"))
            return "ok"
    except Exception:
        return "fail"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    accs = distinct_iso_accessions()
    if len(sys.argv) > 1:
        accs = [a for a in accs if a in set(sys.argv[1:])]
    print(f"[extract_secondary_dssp] {len(accs)} isoform accessions", flush=True)
    counts = {"ok": 0, "skip": 0, "missing": 0, "fail": 0}
    with cf.ThreadPoolExecutor(max_workers=8) as ex:
        for i, status in enumerate(ex.map(one, accs), 1):
            counts[status] += 1
            if i % 200 == 0:
                print(f"  {i}/{len(accs)} {counts}", flush=True)
    print(f"DONE {counts} → {OUT_DIR}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
