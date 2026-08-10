"""Remove the DBD-only TFs (no AD/RD/Bif effector domain) so the DB counts
effector-present TFs only (team scope: AD/RD/Bif). Deletes their per-TF JSONs
and updates index.json: recompute ad/rd/dna/bifunctional/ptms/tfs exactly from
the remaining JSONs; PATCH papers (subtract removed TFs' contribution, since the
per-family papers metric is per-entry and not cleanly recomputable)."""
import json, os, re, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TFDIR = os.path.join(ROOT, "public/mock/tfs")
INDEX = os.path.join(ROOT, "public/mock/index.json")
REMOVE = ["NOBOX", "PRDM10", "ZNF189", "GZF1", "THAP11", "ZKSCAN4"]
EFF = {"AD", "RD", "Bifunctional"}

def numeric_pmids(papers):
    out = set()
    for p in papers or []:
        for pm in re.split(r"[;,]", str(p.get("pmid", "") or "")):
            pm = pm.strip()
            if pm.isdigit():
                out.add(pm)
    return out

# 1. read the TFs to remove; refuse if any actually has an effector (safety)
rem_by_fam = defaultdict(lambda: dict(paper_entries=0))
removed_pmids = set()
for sym in REMOVE:
    path = os.path.join(TFDIR, sym + ".json")
    d = json.load(open(path))
    if any(r.get("type") in EFF for r in d.get("domainReports", [])):
        sys.exit(f"ABORT: {sym} has an effector domain — not DBD-only!")
    fam = d.get("family") or "Unclassified"
    rem_by_fam[fam]["paper_entries"] += len(numeric_pmids(d.get("papers", [])) or [])
    removed_pmids |= numeric_pmids(d.get("papers", []))

# 2. delete the files
for sym in REMOVE:
    os.remove(os.path.join(TFDIR, sym + ".json"))
print(f"deleted {len(REMOVE)} DBD-only TF JSONs: {REMOVE}")

# 3. recompute ad/rd/dna/bif/ptms/tfs exactly from remaining; collect all pmids
fam = defaultdict(lambda: dict(ad=0, rd=0, dna=0, bifunctional=0, ptms=0, tfs=0))
all_pmids = set()
for fn in os.listdir(TFDIR):
    if not fn.endswith(".json"):
        continue
    d = json.load(open(os.path.join(TFDIR, fn)))
    f = d.get("family") or "Unclassified"
    S = fam[f]; S["tfs"] += 1; S["ptms"] += len(d.get("ptms", []))
    for r in d.get("domainReports", []):
        t = (r.get("type") or "").lower()
        if t == "ad": S["ad"] += 1
        elif t == "rd": S["rd"] += 1
        elif t in ("dna", "dbd"): S["dna"] += 1
        elif t == "bifunctional": S["bifunctional"] += 1
    all_pmids |= numeric_pmids(d.get("papers", []))

# 4. update index.json
idx = json.load(open(INDEX))
for f in idx["families"]:
    nm = f["name"]; S = fam[nm]
    f["tfs"] = [t for t in f.get("tfs", []) if t.get("symbol") not in REMOVE]
    c = f["counts"]
    c["ad"] = S["ad"]; c["rd"] = S["rd"]; c["dna"] = S["dna"]
    c["bifunctional"] = S["bifunctional"]; c["ptms"] = S["ptms"]; c["tfs"] = S["tfs"]
    c["papers"] = max(0, c.get("papers", 0) - rem_by_fam[nm]["paper_entries"])
    f["curated"] = S["tfs"]
T = idx["totals"]
T["tfs"] = sum(f["counts"]["tfs"] for f in idx["families"])
T["ad"] = sum(f["counts"]["ad"] for f in idx["families"])
T["rd"] = sum(f["counts"]["rd"] for f in idx["families"])
T["ptms"] = sum(f["counts"]["ptms"] for f in idx["families"])
# totals.papers = distinct: subtract only removed PMIDs that survive nowhere else
orphaned = removed_pmids - all_pmids
T["papers"] = T.get("papers", 0) - len(orphaned)
json.dump(idx, open(INDEX, "w"), ensure_ascii=False, separators=(",", ":"))
print(f"index.json updated: totals tfs={T['tfs']} ad={T['ad']} rd={T['rd']} "
      f"papers={T['papers']} (orphaned pmids removed: {len(orphaned)})")
