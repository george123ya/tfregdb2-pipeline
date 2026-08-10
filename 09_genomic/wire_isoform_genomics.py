"""
Run prot2exon once for every isoform of every curated TF and inline a
per-isoform genomic structure into each per-TF JSON.

Pipeline:
  1. Read TF_completeIDs.xlsx -> for each gene, map "Isoform name" -> ENSP.
  2. Backfill `proteinId` (ENSP) onto each isoform in every per-TF JSON.
  3. Generate a unified BED:
        - canonical AD/RD/DBD claims (already in Unified_Effector_Domains.bed)
        - one row per isoform ENSP with aa range = 1..length (so prot2exon
          returns the structure for the alternate transcripts too)
  4. Run prot2exon --output isoform against that BED.
  5. Parse isoform_structure.tsv:
        - For each (gene, transcript) seen, build a deduped feature list.
        - Carry domain overlays only on the canonical transcript (the only
          one literature reports were mapped against).
  6. Replace each per-TF JSON's `genomicStructure` (canonical) with a
     `genomicStructures: [...]` array — one entry per isoform that prot2exon
     could map. Canonical entry keeps the domain overlays; alternates get
     bare structure.
"""

from __future__ import annotations
import csv
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from typing import Optional

import pandas as pd


# Paths
ROOT = os.environ.get("TFREGDB2_ROOT", os.path.expanduser("~/Desktop/tfregdb2"))
TFS_DIR = os.path.join(ROOT, "public/mock/tfs")
# prot2exon was renamed to fastCDS (Jun 2026) — repo and binary both moved.
# Both binaries are env-overridable and default to the bare command name so
# they resolve via PATH; the legacy prot2exon is kept as a fallback so an
# unmigrated checkout still works.
FASTCDS_BIN = os.environ.get("FASTCDS_BIN", "fastCDS")
LEGACY_PROT2EXON_BIN = os.environ.get("PROT2EXON_BIN", "prot2exon")
PROT2EXON_BIN = next(
    (p for p in (FASTCDS_BIN, LEGACY_PROT2EXON_BIN) if os.path.exists(p)),
    FASTCDS_BIN,
)
P2E_OUT_BASE = os.environ.get(
    "P2E_OUT_BASE", os.path.expanduser("~/Desktop/tfregdb_source/prot2exon")
)
INDEX = os.path.join(P2E_OUT_BASE, "human_primary.idx")
COMBINED_BED = os.path.join(P2E_OUT_BASE, "combined_with_isoforms.bed")
RESULTS_DIR = os.path.join(P2E_OUT_BASE, "results_isoforms")
EXISTING_BED = os.environ.get(
    "EFFECTOR_BED",
    os.path.expanduser("~/Desktop/tfregdb/data_pipeline/Unified_Effector_Domains.bed"),
)
COMPLETE_IDS_XLSX = os.environ.get(
    "COMPLETE_IDS_XLSX",
    os.path.expanduser("~/Desktop/tfregdb_source/tables/TF_completeIDs.xlsx"),
)
CACHE_DIR = os.environ.get(
    "CACHE_DIR", os.path.expanduser("~/Desktop/tfregdb_source/.cache")
)
os.makedirs(CACHE_DIR, exist_ok=True)


FEATURE_TYPE_NORMALIZE = {
    "five_prime_UTR": "5utr",
    "CDS": "cds",
    "intron": "intron",
    "three_prime_UTR": "3utr",
}


def _na_int(v: str) -> Optional[int]:
    if v == "NA" or v == "":
        return None
    try:
        return int(v)
    except ValueError:
        return None


def read_ids_df() -> pd.DataFrame:
    """Read the TF_completeIDs xlsx (no parquet cache — pyarrow isn't always
    available on the host)."""
    return pd.read_excel(COMPLETE_IDS_XLSX)


def build_isoform_name_to_ensp() -> dict[tuple[str, str], str]:
    """Map (gene_symbol_upper, isoform_name) -> clean ENSP (no version)."""
    df = read_ids_df()
    df = df.dropna(subset=["UniProt_ID", "Translation ID"]).copy()
    df["Gene_Upper"] = df["Gene name (HGNC)"].astype(str).str.strip().str.upper()
    df["ENSP_Clean"] = (
        df["Translation ID"].astype(str).str.split(".").str[0].str.strip()
    )
    out: dict[tuple[str, str], str] = {}
    for _, r in df.iterrows():
        iso_name = r.get("Isoform name")
        if pd.isna(iso_name):
            continue
        key = (r["Gene_Upper"], str(iso_name).strip())
        out.setdefault(key, r["ENSP_Clean"])
    return out


def backfill_isoform_ensps() -> dict[str, list[dict]]:
    """Attach `proteinId` to each isoform in every per-TF JSON and return
    {symbol: [isoform_with_ensp, ...]} for the BED writer."""
    iso_map = build_isoform_name_to_ensp()
    print(f"  isoform-name -> ENSP entries: {len(iso_map)}")

    by_tf: dict[str, list[dict]] = {}
    patched = 0
    for fn in sorted(os.listdir(TFS_DIR)):
        if not fn.endswith(".json"):
            continue
        p = os.path.join(TFS_DIR, fn)
        with open(p) as f:
            tf = json.load(f)
        sym = tf["symbol"]
        gene = sym.upper()
        isos = tf.get("isoforms") or []
        if not isos:
            continue
        changed = False
        for i in isos:
            name = i.get("name")
            if not name or i.get("proteinId"):
                continue
            # The canonical isoform's ENSP is always the TF's main UniProt's
            # paired ENSP — already on the TF as ensemblProtein.
            if i.get("canonical") and tf.get("ensemblProtein"):
                i["proteinId"] = tf["ensemblProtein"]
                changed = True
                continue
            # Strip the " (NNN aa)" disambiguation suffix we added earlier
            # to recover the original Lambert isoform name.
            base_name = re.sub(r" \(\d+ aa\)$", "", name)
            ensp = iso_map.get((gene, base_name))
            if ensp:
                i["proteinId"] = ensp
                changed = True
        if changed:
            patched += 1
            with open(p, "w") as f:
                json.dump(tf, f, separators=(",", ":"))
        by_tf[sym] = isos
    print(f"  per-TF JSONs patched with isoform ENSPs: {patched}")
    return by_tf


def write_combined_bed(by_tf: dict[str, list[dict]]):
    """Combine the existing curated-effector BED with one row per isoform
    ENSP using the full 1..length aa range so prot2exon emits a full
    structure for each isoform."""
    seen_ensps_for_full_query: set[str] = set()
    with open(COMBINED_BED, "w") as out:
        # Carry the curated rows verbatim — these drive the canonical domain
        # overlays.
        with open(EXISTING_BED) as src:
            for line in src:
                out.write(line)
        # Add one structure-only row per non-canonical isoform ENSP.
        added = 0
        for sym, isos in by_tf.items():
            for i in isos:
                ensp = i.get("proteinId")
                length = i.get("length")
                if not ensp or not length:
                    continue
                if ensp in seen_ensps_for_full_query:
                    continue
                seen_ensps_for_full_query.add(ensp)
                # Use a synthetic structure query keyed by "{sym}_iso_{ensp}".
                # aa range 1..length asks prot2exon for the full transcript.
                domain_id = f"{sym}_iso_{ensp}"
                out.write(f"{ensp}\t1\t{length}\t{domain_id}\t100.0\t.\n")
                added += 1
        print(f"  isoform full-range BED rows added: {added}")
    print(f"  combined BED written: {COMBINED_BED}")


def run_prot2exon():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    cmd = [
        PROT2EXON_BIN,
        "--index", INDEX,
        "--bed", COMBINED_BED,
        "--out-dir", RESULTS_DIR,
        "--output", "isoform",
    ]
    print(f"  running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def parse_isoform_structures() -> dict[str, dict[str, dict]]:
    """Return {symbol: {ENSP: structure}} from results_isoforms/isoform_structure.tsv.

    Each ENSP→structure entry has the shape consumed by the React component:
        { transcriptId, proteinId, chrom, strand, isCanonical, isManeSelect,
          cdsLength, features }

    Features are deduped by feature_id (the per-domain CDS splits are merged
    back into single feature spans).
    """
    path = os.path.join(RESULTS_DIR, "isoform_structure.tsv")
    # gene_symbol -> ensp -> aggregator
    by_gene_ensp: dict[str, dict[str, dict]] = defaultdict(dict)
    # gene_symbol -> ensp -> feature_id -> agg
    fid_agg: dict[str, dict[str, dict[str, dict]]] = defaultdict(lambda: defaultdict(dict))

    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            sym = (row.get("gene_name") or "").strip()
            ensp = (row.get("protein_id") or "").strip()
            if not sym or not ensp:
                continue
            canon_flag = row.get("is_ensembl_canonical") == "true"
            mane_flag = row.get("is_mane_select") == "true"
            entry = by_gene_ensp[sym].get(ensp)
            if entry is None:
                entry = {
                    "transcriptId": row["transcript_id"],
                    "proteinId": ensp,
                    "chrom": row["chrom"],
                    "strand": row["strand"],
                    "isCanonical": canon_flag,
                    "isManeSelect": mane_flag,
                }
                by_gene_ensp[sym][ensp] = entry
            feat_type = FEATURE_TYPE_NORMALIZE.get(row["feature_type"])
            if not feat_type:
                continue
            fid = row["feature_id"]
            s = _na_int(row["feature_genomic_start"])
            e = _na_int(row["feature_genomic_end"])
            order = _na_int(row.get("feature_order_transcript") or "")
            aa_s = _na_int(row.get("aa_start_encoded") or "") if feat_type == "cds" else None
            aa_e = _na_int(row.get("aa_end_encoded") or "") if feat_type == "cds" else None
            cur = fid_agg[sym][ensp].get(fid)
            if cur is None:
                fid_agg[sym][ensp][fid] = {
                    "type": feat_type,
                    "start": s,
                    "end": e,
                    "order": order,
                    "aaStart": aa_s,
                    "aaEnd": aa_e,
                }
            else:
                if s is not None and (cur["start"] is None or s < cur["start"]):
                    cur["start"] = s
                if e is not None and (cur["end"] is None or e > cur["end"]):
                    cur["end"] = e
                if aa_s is not None and (cur["aaStart"] is None or aa_s < cur["aaStart"]):
                    cur["aaStart"] = aa_s
                if aa_e is not None and (cur["aaEnd"] is None or aa_e > cur["aaEnd"]):
                    cur["aaEnd"] = aa_e

    # Materialise feature lists per gene+ensp.
    for sym, by_ensp in by_gene_ensp.items():
        for ensp, entry in by_ensp.items():
            feats = list(fid_agg[sym][ensp].values())
            feats.sort(key=lambda f: f["order"] if f["order"] is not None else 1e9)
            for f in feats:
                if f["type"] != "cds":
                    f.pop("aaStart", None)
                    f.pop("aaEnd", None)
            entry["features"] = feats
            entry["cdsLength"] = sum(
                (f["end"] - f["start"] + 1)
                for f in feats
                if f["type"] == "cds" and f["start"] is not None and f["end"] is not None
            )
    return by_gene_ensp


def parse_canonical_domains() -> dict[str, list[dict]]:
    """Return {symbol: [domain_overlay, ...]} parsed from the existing
    domain_mapping_summary.tsv (the curated-domain run only)."""
    path = os.path.join(P2E_OUT_BASE, "results", "domain_mapping_summary.tsv")
    by_gene: dict[str, list[dict]] = defaultdict(list)
    if not os.path.exists(path):
        return by_gene
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("status") != "ok":
                continue
            sym = row.get("gene_name") or ""
            if not sym:
                continue
            by_gene[sym].append({
                "id": row["domain_id"],
                "aaStart": _na_int(row["aa_start"]),
                "aaEnd": _na_int(row["aa_end"]),
                "genomicStart": _na_int(row["domain_genomic_start"]),
                "genomicEnd": _na_int(row["domain_genomic_end"]),
            })
    return by_gene


DOMAIN_ID_RE = re.compile(
    r"^[A-Za-z0-9]+_(?P<source>[A-Za-z0-9]+(?:[A-Za-z0-9]+)?)(?:_(?P<type>AD|RD|DNA|Oligo|Bifunctional))?_Row\d+$"
)


# Sources whose name alone implies the domain type — for HT-MPRA / Tycko-style
# screens that don't tag the type separately in the BED id.
SOURCE_TYPE_BY_SUBSTRING: list[tuple[str, str]] = [
    # order matters — more specific first
    ("Tycko2020_NucRep", "RD"),
    ("Tycko2020_TilingRep", "RD"),
    ("Tycko2020_NucAct", "AD"),
    ("Tycko2024_AD", "AD"),
    ("Tycko2024_RD", "RD"),
    ("DelRosso_AD", "AD"),
    ("DelRosso_RD", "RD"),
    ("Alerasool2022_AD", "AD"),
    ("Klaus2023_RD", "RD"),
    ("Arnold2018_AD", "AD"),
]

# Tokens encoded directly in the BED id (after the source prefix).
ID_TOKEN_TYPE: dict[str, str] = {
    "_AD_": "AD",
    "_RD_": "RD",
    "_DNA_": "DNA",
    "_Oligo_": "Oligo",
    "_Bifunctional_": "Bifunctional",
    "_Bif_": "Bifunctional",
}


def classify_domain(
    domain_id: str,
    aa_start: Optional[int],
    aa_end: Optional[int],
    reports: list[dict],
) -> tuple[str, Optional[str]]:
    """Resolve the (type, source) for a curated BED entry.

    Match priority:
      1. The exact `(aaStart, aaEnd)` in `tf.domainReports` — this is the
         strongest match because the build script already classified each
         report with AD/RD/Bifunctional/DNA/Other.
      2. Same (start, end) match with a ±1 aa wobble (Lambert / Ensembl
         numbering occasionally shifts by 1 between sources).
      3. Regex on the BED domain_id — catches sources that embed the type
         (e.g. `_DelRosso_RD_Row123`).
      4. Source-name fallback: first report from the same source.

    Returns ("Other", source) only when none of the above match — used to be
    ~410 false negatives because the regex was the only path.
    """
    m = DOMAIN_ID_RE.match(domain_id or "")
    source = m.group("source") if m else None
    did = domain_id or ""

    # 1 + 2: exact / wobble match by aa coordinates.
    if aa_start is not None and aa_end is not None and reports:
        for r in reports:
            if r.get("start") == aa_start and r.get("end") == aa_end:
                return r.get("type") or "Other", r.get("source") or source
        for r in reports:
            rs = r.get("start")
            re_ = r.get("end")
            if rs is None or re_ is None:
                continue
            if abs(rs - aa_start) <= 1 and abs(re_ - aa_end) <= 1:
                return r.get("type") or "Other", r.get("source") or source

    # 3: explicit type token inside the BED id (e.g. _AD_, _RD_, _Bif_).
    for token, typ in ID_TOKEN_TYPE.items():
        if token in did:
            return typ, source

    # 4: source-name encodes the type (Tycko2020_NucRep → RD, Tycko2024_AD → AD…).
    for src_key, typ in SOURCE_TYPE_BY_SUBSTRING:
        if src_key in did:
            return typ, src_key

    # 5: regex-encoded type from earlier pattern (backstop).
    typ = m.group("type") if (m and m.group("type")) else None
    if typ:
        return typ, source

    # 6: any report from the same source.
    if source and reports:
        for r in reports:
            r_src = r.get("source") or ""
            if r_src.lower().startswith(source.lower()):
                return r.get("type") or "Other", source

    return "Other", source


def derive_aa_exons(
    iso_struct: dict,
    canon_struct: dict,
    canonical_length: int,
) -> list[dict]:
    """Translate an isoform's genomic CDS coordinates into AA-frame ranges
    on the CANONICAL protein, by intersecting genomic CDS spans with the
    canonical's per-CDS aa map. Output drives the AA-tab coverage filter.

    For each isoform CDS feature, find canonical CDS features that overlap
    it genomically; convert the overlap range into canonical aa positions
    via the canonical feature's linear (aa, genomic) mapping; merge the
    resulting aa ranges into a final list. The result is a sequence of
    `{start, end}` aa intervals (1-based, inclusive) on the canonical frame
    showing exactly which canonical residues are encoded by this isoform.
    """
    iso_cds = [
        f for f in iso_struct.get("features", [])
        if f.get("type") == "cds"
    ]
    canon_cds = [
        f for f in canon_struct.get("features", [])
        if f.get("type") == "cds"
        and f.get("aaStart") is not None
        and f.get("aaEnd") is not None
    ]
    if not iso_cds or not canon_cds:
        return []
    strand = canon_struct.get("strand", "+")

    def canon_aa_for(pos: int, c: dict) -> int:
        # Map a genomic position to canonical aa within canonical CDS feature `c`.
        # Each canonical CDS feature has linear (gStart, aaStart) ↔ (gEnd, aaEnd):
        # on +strand: aa = aaStart + (pos - gStart) / 3
        # on −strand: aa = aaStart + (gEnd - pos) / 3
        g_start, g_end = c["start"], c["end"]
        aa_start, aa_end = c["aaStart"], c["aaEnd"]
        # Clamp position to the feature's genomic span.
        pos = max(g_start, min(g_end, pos))
        if strand == "+":
            return aa_start + (pos - g_start) // 3
        return aa_start + (g_end - pos) // 3

    ranges: list[tuple[int, int]] = []
    for ic in iso_cds:
        ig_s, ig_e = ic["start"], ic["end"]
        for cc in canon_cds:
            cg_s, cg_e = cc["start"], cc["end"]
            ov_s = max(ig_s, cg_s)
            ov_e = min(ig_e, cg_e)
            if ov_s > ov_e:
                continue
            a = canon_aa_for(ov_s, cc)
            b = canon_aa_for(ov_e, cc)
            lo, hi = (a, b) if a <= b else (b, a)
            lo = max(1, min(canonical_length, lo))
            hi = max(1, min(canonical_length, hi))
            ranges.append((lo, hi))

    # Sort + merge overlapping / adjacent ranges.
    ranges.sort()
    merged: list[list[int]] = []
    for s, e in ranges:
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [{"start": s, "end": e} for s, e in merged]


def project_iso_effector_to_genomic(
    eff_start: int,
    eff_end: int,
    iso_struct: dict,
) -> tuple[int, int] | None:
    """Map an iso-frame effector aa range [eff_start, eff_end] onto a single
    genomic span using the iso's per-CDS aa map. Returns (gMin, gMax) across
    every iso CDS feature that intersects the aa range, so a domain that
    spans an intron renders as one continuous genomic bar (matching the
    canonical-domain rendering convention)."""
    iso_cds = [
        f for f in iso_struct.get("features", [])
        if f.get("type") == "cds"
        and f.get("aaStart") is not None
        and f.get("aaEnd") is not None
    ]
    if not iso_cds:
        return None
    # prot2exon reports `aa_start_encoded`/`aa_end_encoded` for each CDS
    # exon. When a codon is split across an exon junction, the shared AA
    # appears as the aaEnd of one exon AND the aaStart of the next (e.g.
    # CDS_1 ends at aa 60 and CDS_2 starts at aa 60). Without disambiguating,
    # an effector at aa 1-60 would match BOTH CDS_1 and CDS_2, and the
    # gen_lo..gen_hi span would balloon across the intervening intron
    # (RARA-206: Row1398 aa 1-60 → ended up rendered as aa 1-365 in the
    # genomic view). Award the boundary AA to the earlier exon by bumping
    # later exons' aaStart past the previous aaEnd.
    iso_cds_sorted = sorted(iso_cds, key=lambda c: int(c["aaStart"]))
    iso_cds_eff: list[dict] = []
    prev_end = 0
    for c in iso_cds_sorted:
        adj_as = max(int(c["aaStart"]), prev_end + 1)
        adj_ae = int(c["aaEnd"])
        if adj_ae < adj_as:
            continue
        iso_cds_eff.append({**c, "aaStart": adj_as, "aaEnd": adj_ae})
        prev_end = adj_ae
    iso_cds = iso_cds_eff
    strand = iso_struct.get("strand", "+")
    gen_lo, gen_hi = None, None
    for c in iso_cds:
        c_as, c_ae = int(c["aaStart"]), int(c["aaEnd"])
        ov_s = max(eff_start, c_as)
        ov_e = min(eff_end, c_ae)
        if ov_s > ov_e:
            continue
        g_start = int(c["start"])
        g_end = int(c["end"])
        # Position within CDS feature: iso aa position → genomic position.
        # +strand: g = g_start + (aa - c_as) * 3
        # −strand: g = g_end   − (aa - c_as) * 3
        if strand == "+":
            g_a = g_start + (ov_s - c_as) * 3
            g_b = g_start + (ov_e - c_as) * 3 + 2
        else:
            g_a = g_end - (ov_s - c_as) * 3 - 2
            g_b = g_end - (ov_e - c_as) * 3
        lo, hi = (min(g_a, g_b), max(g_a, g_b))
        gen_lo = lo if gen_lo is None else min(gen_lo, lo)
        gen_hi = hi if gen_hi is None else max(gen_hi, hi)
    if gen_lo is None or gen_hi is None:
        return None
    return gen_lo, gen_hi


def write_back_per_tf(structures: dict[str, dict[str, dict]],
                      doms_summary: dict[str, list[dict]]):
    """For each per-TF JSON, attach `genomicStructures: [...]` (one entry per
    isoform that prot2exon could map). The canonical entry carries domain
    overlays; alternates do not.

    Also: REPLACE each non-canonical isoform's `exons` field — previously a
    synthetic single-exon `[{start: offset+1, end: length}]` assuming N-terminal
    truncation — with a real list of canonical-aa-frame ranges derived from
    the genomic CDS-vs-canonical intersection. This drives the PTM / mutation
    fade in the AA tab correctly even for C-terminally truncated or internally
    spliced isoforms.
    """
    patched = 0
    skipped = 0
    aa_rewritten = 0
    for fn in sorted(os.listdir(TFS_DIR)):
        if not fn.endswith(".json"):
            continue
        p = os.path.join(TFS_DIR, fn)
        with open(p) as f:
            tf = json.load(f)
        sym = tf["symbol"]
        canonical_ensp = tf.get("ensemblProtein")
        per_ensp = structures.get(sym, {})
        if not per_ensp:
            skipped += 1
            continue
        # Classify + dedupe canonical domain overlays.
        raw_doms = doms_summary.get(sym, [])
        deduped: list[dict] = []
        seen_keys: set[tuple] = set()
        for d in raw_doms:
            typ, src = classify_domain(
                d["id"], d.get("aaStart"), d.get("aaEnd"),
                tf.get("domainReports", []),
            )
            key = (typ, d["aaStart"], d["aaEnd"])
            if key in seen_keys:
                continue
            seen_keys.add(key)
            deduped.append({
                "type": typ,
                "source": src,
                "aaStart": d["aaStart"],
                "aaEnd": d["aaEnd"],
                "genomicStart": d["genomicStart"],
                "genomicEnd": d["genomicEnd"],
            })
        deduped.sort(key=lambda x: (x["aaStart"] or 0, x["aaEnd"] or 0))

        # Find the canonical genomic structure so we can map alt-isoform CDS
        # onto canonical aa coords for the AA-frame exons.
        canon_struct = per_ensp.get(canonical_ensp) if canonical_ensp else None

        # Walk this TF's isoforms in the original order. Attach the structure
        # to each isoform that has a known ENSP; carry domain overlays
        # consistent with the AA tab: every iso emits ORIGINAL + DERIVADA
        # effector chips (ALINEAMIENTO is noise filter applied to match
        # the AA-tab rendering), plus per-isoform Pfam DBDs from
        # iso.pfamDBDs. Coords are derived from iso-aa via the iso's CDS
        # feature map.
        gs: list[dict] = []
        for iso in tf.get("isoforms") or []:
            iso_ensp = iso.get("proteinId")
            if not iso_ensp or iso_ensp not in per_ensp:
                continue
            entry = dict(per_ensp[iso_ensp])
            entry["isoformName"] = iso.get("name") or iso_ensp
            entry["isoformLength"] = iso.get("length")
            entry["isoformCanonical"] = bool(iso.get("canonical"))
            # Same source for every row (canonical + alt): project the
            # per-iso effector + Pfam data into genomic coords, then merge
            # same-type overlapping intervals. ALINEAMIENTO filtered for
            # consistency with the AA tab.
            iso_domains: list[dict] = []
            for ed in iso.get("effectorDomains", []) or []:
                if ed.get("origin") == "alieneamiento":
                    continue
                proj = project_iso_effector_to_genomic(
                    int(ed["start"]), int(ed["end"]), entry
                )
                if not proj:
                    continue
                g_lo, g_hi = proj
                iso_domains.append({
                    "type": ed["type"],
                    "source": ed.get("source"),
                    "aaStart": int(ed["start"]),
                    "aaEnd": int(ed["end"]),
                    "genomicStart": g_lo,
                    "genomicEnd": g_hi,
                    "origin": ed.get("origin"),
                })
            for pf in iso.get("pfamDBDs", []) or []:
                proj = project_iso_effector_to_genomic(
                    int(pf["start"]), int(pf["end"]), entry
                )
                if not proj:
                    continue
                g_lo, g_hi = proj
                iso_domains.append({
                    "type": "DNA",
                    "source": pf.get("name", "CIS-BP 2.0"),
                    "aaStart": int(pf["start"]),
                    "aaEnd": int(pf["end"]),
                    "genomicStart": g_lo,
                    "genomicEnd": g_hi,
                    "origin": "pfam",
                })
            # Same-source duplicate dedup only — two distinct curated calls
            # (e.g. canonical AD 1-87 and AD 154-462) must stay as separate
            # chips even when their iso-frame projections touch at an exon
            # boundary, otherwise a 60-aa N-term AD and a 311-aa C-term AD
            # collapse into one chip covering the whole isoform (RARA-206
            # boundary-codon case). Cross-source merging is gone.
            iso_domains.sort(key=lambda d: (d["type"], d["source"] or "", d["genomicStart"]))
            merged: list[dict] = []
            for d in iso_domains:
                if (
                    merged
                    and merged[-1]["type"] == d["type"]
                    and (merged[-1].get("source") or "") == (d.get("source") or "")
                    and d["genomicStart"] <= merged[-1]["genomicEnd"] + 1
                ):
                    merged[-1]["genomicEnd"] = max(merged[-1]["genomicEnd"], d["genomicEnd"])
                    merged[-1]["aaStart"] = min(merged[-1]["aaStart"], d["aaStart"])
                    merged[-1]["aaEnd"] = max(merged[-1]["aaEnd"], d["aaEnd"])
                else:
                    merged.append(dict(d))
            entry["domains"] = merged
            gs.append(entry)

            # Rewrite the AA-frame exons (replaces the synthetic heuristic).
            if iso.get("canonical"):
                # Canonical covers itself end-to-end.
                iso["exons"] = [{"start": 1, "end": tf["length"]}]
            elif canon_struct:
                real_exons = derive_aa_exons(per_ensp[iso_ensp], canon_struct, tf["length"])
                if real_exons:
                    iso["exons"] = real_exons
                    aa_rewritten += 1
                # else: keep the existing synthetic exons as a fallback so the
                # coverage mask doesn't accidentally hide every PTM.

        if not gs:
            skipped += 1
            continue
        tf["genomicStructures"] = gs
        # Keep the old single-structure field for backwards compatibility but
        # have it point at the canonical one in the new array.
        canon_entry = next((g for g in gs if g.get("isoformCanonical")), gs[0])
        tf["genomicStructure"] = canon_entry
        with open(p, "w") as f:
            json.dump(tf, f, separators=(",", ":"))
        patched += 1

    print(f"  TF JSONs patched with multi-isoform structures: {patched}")
    print(f"  TF JSONs skipped (no mapped isoforms): {skipped}")
    print(f"  non-canonical isoforms with AA exons rewritten from real CDS: {aa_rewritten}")


def main():
    print("[1/5] Backfilling per-isoform ENSPs from TF_completeIDs.xlsx")
    by_tf = backfill_isoform_ensps()
    print("[2/5] Writing combined BED (curated rows + isoform queries)")
    write_combined_bed(by_tf)
    print("[3/5] Running prot2exon on the combined BED")
    run_prot2exon()
    print("[4/5] Parsing per-isoform structures")
    structures = parse_isoform_structures()
    n_iso = sum(len(v) for v in structures.values())
    print(f"  genes with structures: {len(structures)}, total isoforms mapped: {n_iso}")
    print("[5/5] Reading curated domain overlays + writing back to TF JSONs")
    doms_summary = parse_canonical_domains()
    write_back_per_tf(structures, doms_summary)


if __name__ == "__main__":
    main()
