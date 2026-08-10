"""Recompute index.json aggregates from the per-TF JSONs. Preserves each family's
totalInFamily (Lambert) and the per-TF metadata list from the existing index;
recomputes curated/tfs + counts (ad/rd/dna/bifunctional/ptms/papers) from the
live JSONs. Run with --validate to check it reproduces the current index without
writing; run with --write to write index.json.

papers = distinct PMIDs (union of each TF's papers[].pmid) per family / overall.
"""
import json, os, sys, re
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TFDIR = os.path.join(ROOT, "public/mock/tfs")
INDEX = os.path.join(ROOT, "public/mock/index.json")

def pmid_norm(v):
    out = []
    for p in re.split(r"[;,]", str(v or "")):
        p = p.strip()
        if p.isdigit():
            out.append(p)
    return out

def compute(tfdir):
    fam = defaultdict(lambda: dict(ad=0, rd=0, dna=0, bifunctional=0, ptms=0,
                                   papers=set(), tfs=set()))
    for fn in os.listdir(tfdir):
        if not fn.endswith(".json"):
            continue
        d = json.load(open(os.path.join(tfdir, fn)))
        f = d.get("family") or "Unclassified"
        sym = d.get("symbol", fn[:-5])
        S = fam[f]
        S["tfs"].add(sym)
        S["ptms"] += len(d.get("ptms", []))
        for r in d.get("domainReports", []):
            t = (r.get("type") or "").lower()
            if t == "ad": S["ad"] += 1
            elif t == "rd": S["rd"] += 1
            elif t in ("dna", "dbd"): S["dna"] += 1
            elif t == "bifunctional": S["bifunctional"] += 1
        for p in d.get("papers", []):
            for pm in pmid_norm(p.get("pmid")):
                S["papers"].add(pm)
    return fam

def fam_counts(S):
    return dict(ad=S["ad"], rd=S["rd"], dna=S["dna"], bifunctional=S["bifunctional"],
                ptms=S["ptms"], papers=len(S["papers"]), tfs=len(S["tfs"]))

if __name__ == "__main__":
    idx = json.load(open(INDEX))
    fam = compute(TFDIR)
    if "--validate" in sys.argv:
        mism = 0
        for f in idx["families"]:
            got = fam_counts(fam[f["name"]]); tgt = f["counts"]
            for k in ("ad", "rd", "dna", "bifunctional", "ptms", "papers", "tfs"):
                if got.get(k) != tgt.get(k):
                    print(f"  MISMATCH {f['name']} {k}: got {got.get(k)} vs index {tgt.get(k)}")
                    mism += 1
        print(f"validation mismatches: {mism}")
        sys.exit(0)
    if "--write" in sys.argv:
        removed = set(a for a in sys.argv[1:] if not a.startswith("--"))
        allpap = set()
        tot = dict(ad=0, rd=0, ptms=0)
        for f in idx["families"]:
            S = fam[f["name"]]
            # drop removed TFs from the per-TF metadata list
            if removed:
                f["tfs"] = [t for t in f.get("tfs", []) if t.get("symbol") not in removed]
            fc = fam_counts(S)
            f["counts"] = {**f.get("counts", {}), **fc}
            f["curated"] = fc["tfs"]
            allpap |= S["papers"]
            tot["ad"] += fc["ad"]; tot["rd"] += fc["rd"]; tot["ptms"] += fc["ptms"]
        T = idx.setdefault("totals", {})
        T["tfs"] = sum(f["counts"]["tfs"] for f in idx["families"])
        T["ad"] = tot["ad"]; T["rd"] = tot["rd"]; T["ptms"] = tot["ptms"]
        T["papers"] = len(allpap)
        json.dump(idx, open(INDEX, "w"), ensure_ascii=False, separators=(",", ":"))
        print(f"wrote index.json: totals tfs={T['tfs']} ad={T['ad']} rd={T['rd']} papers={T['papers']}")
