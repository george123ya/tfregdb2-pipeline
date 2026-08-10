#!/usr/bin/env python3
"""Measure how 'novel' each isoform is vs its canonical sequence, to decide how
redundant per-isoform ConSurf is. A residue is 'mappable' (canonical ConSurf score
could be transferred) if it lies within a 15-mer that occurs exactly in the
canonical; otherwise it's 'novel' (alt exon / frameshift / junction) and only
isoform-specific ConSurf can score it. Pure deletions => novel_frac ~ 0 (redundant).
Stdlib only."""
import json, glob, os, sys
from collections import Counter

K = 15
TF_DIR = "public/mock/tfs"
ISO_DIR = "consurf_isoforms/inputs"

# 1) canonical sequences by gene symbol
canon = {}
for f in glob.glob(os.path.join(TF_DIR, "*.json")):
    try:
        d = json.load(open(f))
    except Exception:
        continue
    s = d.get("symbol"); seq = d.get("sequence")
    if s and seq:
        canon[s] = seq

# 2) read isoform fastas
def read_fasta(path):
    hdr = None; seq = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            hdr = line[1:]
        else:
            seq.append(line.strip())
    return hdr, "".join(seq)

def novel_fraction(iso, can):
    """fraction of iso residues NOT covered by any exact K-mer shared with canonical."""
    n = len(iso)
    if n == 0:
        return 0.0, 0
    kk = min(K, n)
    can_kmers = set(can[i:i+kk] for i in range(0, max(1, len(can)-kk+1))) if len(can) >= kk else set()
    covered = bytearray(n)
    for i in range(0, n-kk+1):
        if iso[i:i+kk] in can_kmers:
            for j in range(i, i+kk):
                covered[j] = 1
    cov = sum(covered)
    return (n - cov) / n, (n - cov)

rows = []
missing_canon = 0
identical = 0
substring = 0
for fa in glob.glob(os.path.join(ISO_DIR, "*.fasta")):
    hdr, iso = read_fasta(fa)
    gene = (hdr or "").split("|")[0]
    can = canon.get(gene)
    base = os.path.basename(fa)[:-6]
    if not can:
        missing_canon += 1
        rows.append((base, gene, len(iso), None, None))
        continue
    if iso == can:
        identical += 1
    elif iso in can:
        substring += 1
    nf, nres = novel_fraction(iso, can)
    rows.append((base, gene, len(iso), nf, nres))

scored = [r for r in rows if r[3] is not None]
N = len(scored)
print(f"=== Isoform novelty vs canonical (k={K}, {N} scored, {missing_canon} no-canonical) ===")
print(f"identical-to-canonical: {identical}   pure-substring(truncation): {substring}")

# histogram of novel fraction
buckets = [(0.0,0.001,"0% (pure deletion/trunc - REDUNDANT)"),
           (0.001,0.02,"0-2% (junction noise - mappable)"),
           (0.02,0.05,"2-5% (minor novel)"),
           (0.05,0.10,"5-10% (some novel)"),
           (0.10,0.25,"10-25% (substantial novel)"),
           (0.25,0.50,"25-50% (heavily novel)"),
           (0.50,1.0001,"50-100% (mostly novel)")]
print("\nnovel_fraction distribution:")
for lo,hi,lab in buckets:
    c = sum(1 for r in scored if lo <= r[3] < hi)
    bar = "#"*int(60*c/max(1,N))
    print(f"  {lab:42s} {c:5d}  {100*c/N:5.1f}%  {bar}")

# residue-level totals (what fraction of ALL isoform residues is novel)
tot_res = sum(r[2] for r in scored)
tot_novel = sum(r[4] for r in scored)
print(f"\ntotal isoform residues: {tot_res:,}")
print(f"novel (un-mappable) residues: {tot_novel:,}  ({100*tot_novel/max(1,tot_res):.1f}%)")
print(f"mappable-from-canonical residues: {tot_res-tot_novel:,}  ({100*(tot_res-tot_novel)/max(1,tot_res):.1f}%)")

# actionable: how many isoforms are essentially redundant at thresholds
for thr in (0.001, 0.02, 0.05):
    c = sum(1 for r in scored if r[3] < thr)
    print(f"\n  isoforms with novel_frac < {thr*100:.1f}%: {c} ({100*c/N:.1f}%) "
          f"-> could MAP canonical instead of running ConSurf")

# dump per-isoform CSV for downstream filtering
with open("consurf_isoforms/novelty.csv","w") as out:
    out.write("isoform,gene,length,novel_fraction,novel_residues\n")
    for b,g,l,nf,nr in rows:
        out.write(f"{b},{g},{l},{'' if nf is None else f'{nf:.4f}'},{'' if nr is None else nr}\n")
print("\nwrote consurf_isoforms/novelty.csv (per-isoform)")
