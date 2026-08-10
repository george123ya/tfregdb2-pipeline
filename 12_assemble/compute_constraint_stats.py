"""Aggregate ED-vs-DBD variant-constraint + Figure-1 stats across all TFs.

Reads the per-TF mocks (domains, variants) + per-protein pLDDT and writes a
single small JSON the Statistics page renders. This is the plan's "money
figure": are EDs / DBDs depleted of common variation (gnomAD), enriched for
pathogenic variants (ClinVar), or recurrent in cancer (COSMIC)? Plus the
positional bias of EDs vs DBDs and the fraction of each that sits in
low-confidence / disordered (pLDDT < 50) regions.

Region classes per residue (domains may overlap, so a residue can be in both):
  ED  = AD / RD / Bifunctional curated effector domains
  DBD = Pfam DNA-binding domains (type "DNA")
  Other = everything else (the protein backbone outside curated domains)

Output: public/mock/stats_constraint.json  (read by src/pages/statistics.astro)
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

TFS_DIR = Path("public/mock/tfs")
VAR_DIR = Path("public/mock/variants")
PLDDT_DIR = Path("public/mock/plddt")
OUT = Path("public/mock/stats_constraint.json")

ED_TYPES = {"AD", "RD", "Bifunctional", "Bif"}
DBD_TYPES = {"DNA"}
PATHOGENIC = {"pathogenic", "likely_pathogenic"}
DISORDER_PLDDT = 50.0
N_BINS = 10  # positional-bias deciles


def blank():
    return {
        "residues": 0,
        "disordered": 0,
        "clinvarTotal": 0,
        "clinvarPath": 0,
        "gnomadTotal": 0,
        "gnomadMissense": 0,
        "cosmicVariants": 0,
        "cosmicSamples": 0,
        "missenseAll": 0,
    }


def main() -> int:
    classes = {"AD": blank(), "RD": blank(), "ED": blank(), "DBD": blank(), "Other": blank()}
    positional = {"ED": [0] * N_BINS, "DBD": [0] * N_BINS}
    domain_lengths = {"ED": [], "DBD": []}
    cov = {"edTFs": 0, "dbdTFs": 0, "bothTFs": 0}
    n_tf = n_tf_with_var = 0
    seen_types: Counter = Counter()

    # Fig-8 (variant constraint): per-TF mutation density by source + region.
    SRCS = ("gnomAD", "ClinVar", "COSMIC")
    curve_vals = {s: {r: [] for r in ("full", "AD", "RD", "DBD")} for s in SRCS}
    scatter = {s: [] for s in SRCS}
    scatter_diag = {s: {"above": 0, "below": 0, "equal": 0} for s in SRCS}

    for fp in sorted(TFS_DIR.glob("*.json")):
        tf = json.loads(fp.read_text())
        L = int(tf.get("length") or 0)
        if L <= 0:
            continue
        n_tf += 1
        symbol = tf.get("symbol")
        upro = (tf.get("uniprot") or "").strip()

        ed = bytearray(L)
        dbd = bytearray(L)
        ad = bytearray(L)
        rd = bytearray(L)
        has_ed = has_dbd = False
        for d in tf.get("domains") or []:
            t = d.get("type")
            seen_types[t] += 1
            s = max(1, int(d.get("start") or 0))
            e = min(L, int(d.get("end") or 0))
            if e < s:
                continue
            mid_frac = ((s + e) / 2 - 1) / max(1, L - 1)
            b = min(N_BINS - 1, max(0, int(mid_frac * N_BINS)))
            if t in ED_TYPES:
                is_ad = t in ("AD", "Bifunctional", "Bif")
                is_rd = t in ("RD", "Bifunctional", "Bif")
                for i in range(s - 1, e):
                    ed[i] = 1
                    if is_ad:
                        ad[i] = 1
                    if is_rd:
                        rd[i] = 1
                has_ed = True
                positional["ED"][b] += 1
                domain_lengths["ED"].append(e - s + 1)
            elif t in DBD_TYPES:
                for i in range(s - 1, e):
                    dbd[i] = 1
                has_dbd = True
                positional["DBD"][b] += 1
                domain_lengths["DBD"].append(e - s + 1)

        if has_ed:
            cov["edTFs"] += 1
        if has_dbd:
            cov["dbdTFs"] += 1
        if has_ed and has_dbd:
            cov["bothTFs"] += 1

        # pLDDT for the disorder overlap (canonical, by accession).
        plddt = None
        pf = PLDDT_DIR / f"{upro}.json"
        if pf.exists():
            try:
                arr = json.loads(pf.read_text())
                if isinstance(arr, list) and len(arr) >= L:
                    plddt = arr
            except Exception:
                plddt = None

        def cls_of(i: int):
            out = []
            if ad[i]:
                out.append("AD")
            if rd[i]:
                out.append("RD")
            if ed[i]:
                out.append("ED")
            if dbd[i]:
                out.append("DBD")
            if not ed[i] and not dbd[i]:
                out.append("Other")
            return out

        for i in range(L):
            for c in cls_of(i):
                classes[c]["residues"] += 1
                if plddt is not None and plddt[i] < DISORDER_PLDDT:
                    classes[c]["disordered"] += 1

        # Variants — pos is 1-based canonical AA (0 = unmapped / non-coding).
        vfp = VAR_DIR / f"{symbol}.json"
        if not vfp.exists():
            continue
        n_tf_with_var += 1
        # Fig-8 "mutation density" — count of variants / region length in
        # NUCLEOTIDES (aa × 3), per the 2021 paper (regions concatenated).
        # NOTE: the AF filter + disjoint-region tweaks were removed 2026-06-04 —
        # the figure is hidden pending the paper; revisit the metric then.
        tfc = {s: {r: 0 for r in ("full", "AD", "RD", "DBD", "ED")} for s in SRCS}
        for v in json.loads(vfp.read_text()).get("variants") or []:
            p = int(v.get("pos") or 0)
            if p < 1 or p > L:
                continue
            src = v.get("source")
            cons = v.get("consequence")
            tag = ((v.get("clinvar") or {}).get("significanceTag") or "").lower()
            cnt = int(v.get("count") or 1)
            for c in cls_of(p - 1):
                k = classes[c]
                if cons == "missense":
                    k["missenseAll"] += 1
                if src == "ClinVar":
                    k["clinvarTotal"] += 1
                    if tag in PATHOGENIC:
                        k["clinvarPath"] += 1
                elif src == "gnomAD":
                    k["gnomadTotal"] += 1
                    if cons == "missense":
                        k["gnomadMissense"] += 1
                elif src == "COSMIC":
                    k["cosmicVariants"] += 1
                    k["cosmicSamples"] += cnt
            # Fig-8 variant subsets, per the paper: gnomAD = non-synonymous SNVs
            # (missense + nonsense, no frequency filter), ClinVar = Pathogenic +
            # Likely-pathogenic, COSMIC = non-synonymous (proxy — the local dump
            # lacks COSMIC's FATHMM "pathogenic" flag). Count every variant.
            nonsyn = cons in ("missense", "nonsense")
            if src == "gnomAD":
                qual = nonsyn
            elif src == "ClinVar":
                qual = tag in PATHOGENIC
            elif src == "COSMIC":
                qual = nonsyn
            else:
                qual = False
            if qual and src in tfc:
                j = p - 1
                tc = tfc[src]
                tc["full"] += 1
                if ad[j]:
                    tc["AD"] += 1
                if rd[j]:
                    tc["RD"] += 1
                if dbd[j]:
                    tc["DBD"] += 1
                if ed[j]:
                    tc["ED"] += 1
        # Per-TF density = variants / nucleotide length (aa × 3) → curves + scatter.
        lens = {"full": L, "AD": sum(ad), "RD": sum(rd), "DBD": sum(dbd)}
        lenED = sum(ed)
        for s in SRCS:
            for r in ("full", "AD", "RD", "DBD"):
                if lens[r] > 0:
                    curve_vals[s][r].append(tfc[s][r] / (lens[r] * 3))
            if lenED > 0 and lens["DBD"] > 0:
                e_ = tfc[s]["ED"] / (lenED * 3)
                d_ = tfc[s]["DBD"] / (lens["DBD"] * 3)
                scatter[s].append({"s": symbol, "e": round(e_, 5), "d": round(d_, 5)})
                if e_ > d_:
                    scatter_diag[s]["above"] += 1
                elif d_ > e_:
                    scatter_diag[s]["below"] += 1
                else:
                    scatter_diag[s]["equal"] += 1

    # Derived per-residue densities + fractions (the chart inputs).
    derived = {}
    for c, k in classes.items():
        r = max(1, k["residues"])
        derived[c] = {
            "residues": k["residues"],
            "disorderFrac": round(k["disordered"] / r, 4),
            "gnomadMissensePerKres": round(k["gnomadMissense"] / r * 1000, 2),
            "cosmicSamplesPerRes": round(k["cosmicSamples"] / r, 3),
            "missensePerKres": round(k["missenseAll"] / r * 1000, 2),
            "clinvarPathFrac": round(k["clinvarPath"] / max(1, k["clinvarTotal"]), 4),
            "clinvarTotal": k["clinvarTotal"],
            "clinvarPath": k["clinvarPath"],
        }

    # Domain-length histogram (20-aa bins, 0..200, last bin = 200+).
    LBIN_W = 20
    LBIN_N = 11  # bins 0-20 … 180-200, 200+
    def length_hist(vals):
        h = [0] * LBIN_N
        for v in vals:
            h[min(LBIN_N - 1, int(v) // LBIN_W)] += 1
        return h

    # Fig-8 curves: log10-binned distribution of per-TF mutation density.
    import math
    # Gaussian KDE over log10(density) (matches the paper's smooth "density of
    # TFs"). density==0 TFs are placed at a sentinel position just left of the
    # axis (ZEROPOS, labelled "0") so they form a SMOOTH left-edge peak that's
    # part of the same curve — not a detached spike.
    LOG_LO, LOG_HI, NB = -3.5, 0.5, 72
    ZEROPOS = -4.1
    GRID_LO, GRID_HI = -4.6, 0.5
    grid = [GRID_LO + i * (GRID_HI - GRID_LO) / (NB - 1) for i in range(NB)]
    BW = 0.16
    kinv = -0.5 / (BW * BW)
    def log_hist(vals):
        pos = [math.log10(v) if v > 0 else ZEROPOS for v in vals]
        if not pos:
            return [0.0] * NB
        # NORMALISE by TF count → a proper "density of TFs" (area comparable
        # across regions), so AD isn't dwarfed just because fewer TFs have one.
        n = len(pos)
        return [round(sum(math.exp(kinv * (g - x) * (g - x)) for x in pos) / n * 100, 3) for g in grid]
    fig = {
        "sources": list(SRCS),
        "regions": ["full", "AD", "RD", "DBD"],
        "logLo": GRID_LO, "logHi": GRID_HI, "nBins": NB, "zeroPos": ZEROPOS,
        "metric": {"gnomAD": "non-syn variants / nt", "ClinVar": "pathogenic / nt", "COSMIC": "tumour non-syn / nt"},
        "curves": {s: {r: log_hist(curve_vals[s][r]) for r in ("full", "AD", "RD", "DBD")} for s in SRCS},
        "zeros": {s: {r: sum(1 for v in curve_vals[s][r] if v <= 0) for r in ("full", "AD", "RD", "DBD")} for s in SRCS},
        "scatter": {s: scatter[s] for s in SRCS},
        "scatterDiag": scatter_diag,
    }

    out = {
        "nTFs": n_tf,
        "nTFsWithVariants": n_tf_with_var,
        "coverage": cov,
        "fig": fig,
        "classes": classes,
        "derived": derived,
        "positional": positional,
        "domainLengthMedian": {
            c: (sorted(v)[len(v) // 2] if v else 0) for c, v in domain_lengths.items()
        },
        "domainCount": {c: len(v) for c, v in domain_lengths.items()},
        "domainLengthHist": {
            "binWidth": LBIN_W,
            "nBins": LBIN_N,
            "ED": length_hist(domain_lengths["ED"]),
            "DBD": length_hist(domain_lengths["DBD"]),
        },
    }
    OUT.write_text(json.dumps(out, separators=(",", ":")))

    print(f"[constraint] {n_tf} TFs ({n_tf_with_var} with variants) · domain types seen: {dict(seen_types)}")
    for c in ("ED", "DBD", "Other"):
        d = derived[c]
        print(f"  {c:5s}  res={d['residues']:>7}  disorder={d['disorderFrac']*100:5.1f}%  "
              f"gnomAD-mis/kres={d['gnomadMissensePerKres']:6.2f}  "
              f"ClinVar%path={d['clinvarPathFrac']*100:5.1f}  "
              f"COSMIC smp/res={d['cosmicSamplesPerRes']:6.3f}")
    print(f"  → {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
