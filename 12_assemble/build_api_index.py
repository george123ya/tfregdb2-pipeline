"""Build the dynamic-API data layer.

Two outputs power the Cloudflare Pages Functions API (functions/api/v1/):

  public/mock/api/domains.json        — ALL curated domains (≈7.6k rows),
      denormalised with gene / uniprot / type / subtype / family / coords.
      Tiny, so the /domains endpoint filters it in-memory at the edge.

  data/variants/api_in_domains.csv    — every variant that lands inside a
      curated domain (≈1.2M rows), denormalised with the domain's
      type/subtype/family. Too big to scan per-request, so it is imported into
      Cloudflare D1 and the /variants endpoint runs SQL (see functions/README).
      A matching api_in_domains.sqlite is also written for local testing.

Region/subtype mapping:
  AD→activation, RD→repression, Bifunctional→bifunctional  → type ED
  DNA→dna                                                   → type DBD
"""

from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

TFS = Path("public/mock/tfs")
VARS = Path("public/mock/variants")
OUT_DOMAINS = Path("public/mock/api/domains.json")
OUT_PTMS = Path("public/mock/api/ptms.json")
OUT_XREF = Path("public/mock/api/xref.json")
OUT_CSV = Path("data/variants/api_in_domains.csv")
OUT_SQLITE = Path("data/variants/api_in_domains.sqlite")

SUBTYPE = {"AD": "activation", "RD": "repression", "Bifunctional": "bifunctional", "DNA": "dna"}
TYPE = {"AD": "ED", "RD": "ED", "Bifunctional": "ED", "DNA": "DBD"}
PATHOGENIC = {"pathogenic", "likely_pathogenic"}
# Sources kept in the queryable /variants table (D1). gnomAD is now INCLUDED
# (the API supports source=gnomAD + af_min); set API_AF_MIN to drop rare
# gnomAD variants if the D1 row count needs trimming for the free tier.
API_SKIP_SOURCES: set[str] = set()
API_AF_MIN = 0.001  # free-tier D1: keep ClinVar/COSMIC + gnomAD AF>=0.1% (set 0 for all, paid D1)

# PubMed publication years per PMID (scripts inline fetch → data/pmid_years.json),
# used to tag each curated domain with the year(s) it was reported.
try:
    PMID_YEARS: dict[str, int] = {
        str(k): int(v) for k, v in json.loads(Path("data/pmid_years.json").read_text()).items()
    }
except Exception:
    PMID_YEARS = {}

CSV_COLS = [
    "gene", "uniprot", "pos", "aa_ref", "aa_alt", "consequence", "source",
    "clin_sig", "af", "count", "domain_id", "domain_type", "domain_subtype",
    "family", "tf_family",
]


def main() -> int:
    OUT_DOMAINS.parent.mkdir(parents=True, exist_ok=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    domains_out: list[dict] = []
    ptms_out: list[dict] = []
    xref_out: list[dict] = []
    n_var_rows = 0

    with OUT_CSV.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLS)

        for fp in sorted(TFS.glob("*.json")):
            tf = json.loads(fp.read_text())
            L = int(tf.get("length") or 0)
            if L <= 0:
                continue
            sym = tf.get("symbol")
            upro = (tf.get("uniprot") or "").strip()
            ensp = tf.get("ensemblProtein")
            tf_family = tf.get("family")

            # Per-report (start, end, year) so each consensus domain can be
            # tagged with the earliest/latest publication year covering it.
            rep_years: list[tuple[int, int, int]] = []
            for r in tf.get("domainReports") or []:
                yr = PMID_YEARS.get(str(r.get("pmid") or "").split(".")[0])
                rs, re_ = r.get("start"), r.get("end")
                if yr and rs and re_:
                    rep_years.append((int(rs), int(re_), yr))

            # Per-residue → list of domains covering it (for variant tagging).
            covering: list[list[dict]] = [[] for _ in range(L)]
            for d in tf.get("domains") or []:
                t = d.get("type")
                if t not in TYPE:
                    continue
                s = max(1, int(d.get("start") or 0))
                e = min(L, int(d.get("end") or 0))
                if e < s:
                    continue
                # Publication year range = years of reports overlapping this
                # domain (so /domains?year_min&year_max can filter by when a
                # region was reported). DBDs (Pfam) usually have none.
                yrs = [y for (rs, re_, y) in rep_years if rs <= e and re_ >= s]
                dom = {
                    "id": f"{TYPE[t]}:{upro}:{s}-{e}",
                    "symbol": sym,
                    "uniprot": upro,
                    "ensp": ensp,
                    "type": TYPE[t],
                    "subtype": SUBTYPE[t],
                    "family": d.get("name") or d.get("family") or None,
                    "tf_family": tf_family,
                    "start": s,
                    "end": e,
                    "length": e - s + 1,
                    "confidence": d.get("confidence"),
                    "pmid": d.get("pmid"),
                    "year": min(yrs) if yrs else None,
                    "year_max": max(yrs) if yrs else None,
                }
                domains_out.append(dom)
                for i in range(s - 1, e):
                    covering[i].append(dom)

            # ID crosswalk (for /xref and /proteins/{uniprot}).
            xref_out.append({
                "symbol": sym,
                "name": tf.get("name"),
                "uniprot": upro,
                "ensg": tf.get("ensemblGene"),
                "ensp": ensp,
                "length": L,
            })

            # PTMs, tagged with the curated domain(s) they fall inside.
            for pt in tf.get("ptms") or []:
                p = int(pt.get("pos") or 0)
                doms = covering[p - 1] if 1 <= p <= L else []
                ptms_out.append({
                    "gene": sym,
                    "uniprot": upro,
                    "pos": p,
                    "residue": pt.get("residue"),
                    "type": pt.get("type"),
                    "pmid": pt.get("pmid"),
                    "domains": [{"type": d["type"], "subtype": d["subtype"], "family": d["family"]} for d in doms],
                })

            vf = VARS / f"{sym}.json"
            if not vf.exists():
                continue
            for v in json.loads(vf.read_text()).get("variants") or []:
                if v.get("source") in API_SKIP_SOURCES:
                    continue
                p = int(v.get("pos") or 0)
                if p < 1 or p > L:
                    continue
                doms = covering[p - 1]
                if not doms:
                    continue
                tag = ((v.get("clinvar") or {}).get("significanceTag") or "") or None
                af = (v.get("gnomad") or {}).get("af")
                # Optional gnomAD trim for free-tier D1: drop rare gnomAD below
                # API_AF_MIN. ClinVar/COSMIC (clinical/cancer) are always kept.
                if API_AF_MIN and v.get("source") == "gnomAD" and (af or 0) < API_AF_MIN:
                    continue
                base = [
                    sym, upro, p, v.get("aaRef"), v.get("aaAlt"),
                    v.get("consequence"), v.get("source"), tag, af, int(v.get("count") or 1),
                ]
                for dom in doms:
                    w.writerow(base + [dom["id"], dom["type"], dom["subtype"], dom["family"], tf_family])
                    n_var_rows += 1

    OUT_DOMAINS.write_text(json.dumps(domains_out, separators=(",", ":")))
    OUT_PTMS.write_text(json.dumps(ptms_out, separators=(",", ":")))
    OUT_XREF.write_text(json.dumps(xref_out, separators=(",", ":")))

    # Local SQLite mirror for testing the /variants SQL the D1 endpoint runs.
    if OUT_SQLITE.exists():
        OUT_SQLITE.unlink()
    con = sqlite3.connect(OUT_SQLITE)
    con.execute(
        "CREATE TABLE v (gene TEXT, uniprot TEXT, pos INT, aa_ref TEXT, aa_alt TEXT, "
        "consequence TEXT, source TEXT, clin_sig TEXT, af REAL, count INT, "
        "domain_id TEXT, domain_type TEXT, domain_subtype TEXT, family TEXT, tf_family TEXT)"
    )
    AF_COL = CSV_COLS.index("af")
    with OUT_CSV.open() as fh:
        r = csv.reader(fh)
        next(r)
        # Coerce empty af ("" for ClinVar/COSMIC) → NULL so the af_min / af.lte
        # filters compare numerically (in SQLite, text sorts above any number,
        # so "" would otherwise spuriously satisfy `af >= 0.001`).
        def rows():
            for row in r:
                row[AF_COL] = float(row[AF_COL]) if row[AF_COL] else None
                yield row
        con.executemany("INSERT INTO v VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows())
    con.execute("CREATE INDEX ix_dom ON v(domain_type, domain_subtype)")
    con.execute("CREATE INDEX ix_gene ON v(gene)")
    con.execute("CREATE INDEX ix_cons ON v(consequence)")
    con.execute("CREATE INDEX ix_src ON v(source)")
    con.execute("CREATE INDEX ix_tffam ON v(tf_family)")
    con.commit()
    con.close()

    print(f"[api-index] domains: {len(domains_out):,} → {OUT_DOMAINS} "
          f"({OUT_DOMAINS.stat().st_size/1024:.0f} KB)")
    print(f"[api-index] ptms: {len(ptms_out):,} → {OUT_PTMS} "
          f"({OUT_PTMS.stat().st_size/1024:.0f} KB) · xref: {len(xref_out):,} → {OUT_XREF}")
    print(f"[api-index] in-domain variant rows: {n_var_rows:,} → {OUT_CSV} "
          f"({OUT_CSV.stat().st_size/1e6:.1f} MB) + {OUT_SQLITE.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
