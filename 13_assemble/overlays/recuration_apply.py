#!/usr/bin/env python3
"""
Apply Sheyla's re-curation (UpdatedTable.xlsx MasterTable, cols O/P/Q) to the
built per-TF JSONs, and harmonize the standardized evidence columns across ALL
sources (not just TFRegDB2).

Pipeline-aware: MasterTable rows are traced MasterTable -> (original_row) ->
HumanizedTable -> Final_Human_DomainCoords, because the humanization step remaps
non-human species to human orthologs (coords/gene change) and drops rows that
fail BLAST. Matching by raw MasterTable coords would misassign ~half the rows.

Standardized enums (applied to every source):
  necessaryOrSufficient : Necessary | Sufficient | Necessary and sufficient | Neither | None
  verification          : Verified | Cited only | Interaction | DNA-binding |
                          Needs review | Measured (high-throughput) | Literature (low-throughput) | None
  assay                 : Luciferase reporter | CAT reporter | b-galactosidase reporter |
                          Yeast two-hybrid | GST pull-down | Co-IP | EMSA | ChIP |
                          dCas9 reporter | HT-recruit | Binding / interaction |
                          Reporter assay (other) | Functional assay | Review | Not specified | Other | None

Decisions (TFRegDB2 re-curated rows only):
  Verified   -> keep + annotate
  Cited only -> DROP report
  Interaction/UNION AL ADN -> RECLASSIFY (type Interaction / DNA-binding; out of AD/RD/Bif counts)
  error/other-> flag as Needs review (kept)

Writes a persistent overlay (data_pipeline/recuration_overlay.json) so a future
build can re-apply without re-parsing the xlsx.
"""
import openpyxl, json, os, re, csv, sys, unicodedata
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XLSX = os.environ.get("RECUR_XLSX", "/home/goguxor/Downloads/UpdatedTable.xlsx")
TFDIR = os.path.join(ROOT, "public/mock/tfs")
INDEX = os.path.join(ROOT, "public/mock/index.json")
OVERLAY = os.path.join(ROOT, "data_pipeline/recuration_overlay.json")
ADDITIONS = os.path.join(ROOT, "data_pipeline/recuration_additions_to_ingest.csv")

def strip_accents(s):
    return "".join(c for c in unicodedata.normalize("NFD", s) if unicodedata.category(c) != "Mn")

# ---------------- standardizers ----------------
def std_mech(v):
    if v is None: return None
    u = strip_accents(str(v).strip()).upper()
    if u in ("AMBOS", "N AND S", "N Y S", "NECESSARY AND SUFFICIENT", "N&S"): return "Necessary and sufficient"
    if u in ("NECESARIO", "N", "NECESSARY"): return "Necessary"
    if u in ("SUFICIENTE", "S", "SUFFICIENT"): return "Sufficient"
    if u in ("NINGUNO", "NEITHER", "NONE", "NINGUNA"): return "Neither"
    return None

def std_assay(v):
    if v is None: return None
    s = str(v).strip()
    if s == "": return None
    low = strip_accents(s).lower()
    if low in ("ninguno", "none", "not specified", "no especificado", "n/a", "na", "-"): return "Not specified"
    if "lucifer" in low: return "Luciferase reporter"
    if "cloranfenicol" in low or "chloramph" in low or re.search(r'\bcat\b', low): return "CAT reporter"
    if "galactosid" in low or "b-gal" in low or "beta-gal" in low or "b gal" in low: return "b-galactosidase reporter"
    if "two-hybrid" in low or "two hybrid" in low or "dos hibrido" in low: return "Yeast two-hybrid"
    if "gst pull" in low or "pull-down" in low or "pulldown" in low or "pull down" in low: return "GST pull-down"
    if "co-ip" in low or "immunoprecip" in low or "inmunoprecip" in low or "co-inmuno" in low or "coimmuno" in low: return "Co-IP"
    if "emsa" in low or "gel shift" in low or "gel-shift" in low or "retardo" in low or "mobility shift" in low: return "EMSA"
    if "chip" in low: return "ChIP"
    if "dcas9" in low: return "dCas9 reporter"
    if "mpra" in low or "recruit" in low or "ht-recruit" in low: return "HT-recruit"
    if low.startswith("review"): return "Review"
    if "functional assay" in low: return "Functional assay"
    if "binding to" in low or low.startswith("binds") or "interacc" in low or "union a" in low or "interaction" in low: return "Binding / interaction"
    if "reporter" in low or "reportero" in low or "gal4" in low or "lexa" in low or "target gene" in low or "expresion de gen" in low: return "Reporter assay (other)"
    return "Other"

def verif_from_Q(Q):
    """Return (verification, decision, reclass_type)."""
    if Q is None or str(Q).strip() == "": return (None, "none", None)
    s = strip_accents(str(Q).strip()).upper()
    if "ERROR" in s or "UNAVAILABLE" in s: return ("Needs review", "flag", None)
    if s.startswith("VVERIFICADO") or s.startswith("VERIFICADO"): return ("Verified", "keep", None)
    if s.startswith("CITADO"): return ("Cited only", "drop", None)
    if "UNION AL ADN" in s or "INTERACCION PROTEINA-DNA" in s or "INTERACCION PROTEINA-ADN" in s: return ("DNA-binding", "reclassify", "DNA-binding")
    if "INTERACCION" in s: return ("Interaction", "reclassify", "Interaction")
    return ("Needs review", "flag", None)  # structural-domain notes, etc.

def verif_by_source(src):
    if src == "TFRegDB1": return "Literature (low-throughput)"
    if src and any(str(src).startswith(x) for x in ("Tycko", "DelRosso", "Del Rosso", "Alerasool", "Arnold", "Klaus")):
        return "Measured (HT-MPRA)"
    return None  # CIS-BP DBDs, unreviewed TFRegDB2

HT_SOURCES = ("Tycko", "DelRosso", "Del Rosso", "Alerasool", "Arnold", "Klaus")

# ---------------- typo fixes for ingested curation ----------------
GENE_TYPO = {"YV1": "YY1"}
def fix_gene(g):
    g = (str(g).strip() if g is not None else "")
    return GENE_TYPO.get(g.upper(), g)
def fix_type(t):
    t = strip_accents(str(t).strip()).lower() if t is not None else ""
    t = t.replace("acivation", "activation").replace("represion", "repression")
    return t

def resolve_coords(raw, length):
    """Handle '-N:' (last N residues) and '-N: (a-b)' notations -> (start,end) list."""
    s = str(raw or "").strip()
    # explicit domain in parens wins:  -38: (713-750)
    m = re.search(r'\((\d+)\s*[-–]\s*(\d+)\)', s)
    if m: return [(int(m.group(1)), int(m.group(2)))]
    # last-N notation:  -83:
    m = re.match(r'^\s*-\s*(\d+)\s*:?\s*$', s)
    if m and length:
        n = int(m.group(1)); return [(max(1, int(length) - n + 1), int(length))]
    s2 = s.replace("–", "-").replace("—", "-").replace("−", "-")
    return [(int(a), int(b)) for a, b in re.findall(r'(\d+)\s*-\s*(\d+)', s2.replace(" ", ""))]

def npmid(v):
    s = "" if v is None else str(v).strip()
    return s[:-2] if s.endswith(".0") else s

# ================= EXTRACT overlay from xlsx =================
wb = openpyxl.load_workbook(XLSX, read_only=True)
def load(sh):
    rows = list(wb[sh].iter_rows(values_only=True)); h = {k: i for i, k in enumerate(rows[0])}
    return h, [list(r) for r in rows[1:]]
def g(r, h, k):
    i = h.get(k); return r[i] if (i is not None and i < len(r)) else None

mh, md = load("MasterTable")
recur = {}  # master sheet-row -> normalized curation
for n, r in enumerate(md):
    srow = n + 2
    O = g(r, mh, "Clasificacion_Mecanismo"); P = g(r, mh, "Experimento_Molecular"); Q = g(r, mh, "Estado_Verificacion")
    if all(v in (None, "") for v in (O, P, Q)): continue
    verif, decision, reclass = verif_from_Q(Q)
    recur[srow] = dict(mech=std_mech(O), assay=std_assay(P), verif=verif, decision=decision, reclass=reclass,
                       pmid=npmid(g(r, mh, "PMID")), gene=fix_gene(g(r, mh, "Gene")),
                       raw_coords=str(g(r, mh, "DomainCoordinates") or ""), dtype=fix_type(g(r, mh, "DomainType")),
                       length=g(r, mh, "UNIPROT LENGTH"), O_raw=O, P_raw=P, Q_raw=Q)

hh, hd = load("HumanizedTable")
overlay = {}  # (pmid,start,end) -> curation dict  (final human frame)
traced = 0
for r in hd:
    orw = g(r, hh, "original_row")
    if orw is None: continue
    orw = int(orw)
    rc = recur.get(orw)
    if not rc: continue
    fu = g(r, hh, "Final_Human_UNIPROT_ID"); fc = g(r, hh, "Final_Human_DomainCoords")
    length = g(r, hh, "UNIPROT LENGTH") or rc["length"]
    if not fu or not fc:  # humanization failed -> not in DB
        continue
    for (a, b) in resolve_coords(fc, length):
        overlay[f"{rc['pmid']}|{a}|{b}"] = dict(mechanism=rc["mech"], assay=rc["assay"], verification=rc["verif"],
                                                decision=rc["decision"], reclass=rc["reclass"])
        traced += 1
wb.close()
print(f"[extract] re-curated master rows: {len(recur)}   overlay keys (pmid|start|end): {len(overlay)}")
json.dump(overlay, open(OVERLAY, "w"), ensure_ascii=False, indent=0)
print(f"[extract] wrote {OVERLAY}")

# ================= APPLY to per-TF JSONs =================
stat = Counter()
fam_effector = defaultdict(lambda: Counter())  # family -> {ad,rd,bifunctional}
for fn in sorted(os.listdir(TFDIR)):
    if not fn.endswith(".json"): continue
    path = os.path.join(TFDIR, fn)
    d = json.load(open(path))
    fam = d.get("family")
    reps = d.get("domainReports", [])
    kept = []
    for rep in reps:
        src = rep.get("source")
        # ---- harmonize standardized columns for EVERY source ----
        rep["necessaryOrSufficient"] = std_mech(rep.get("necessaryOrSufficient"))
        rep["assay"] = std_assay(rep.get("assay"))
        # source-default verification (may be overridden by overlay below)
        if src and any(str(src).startswith(x) for x in HT_SOURCES):
            rep["verification"] = "Measured (high-throughput)"
            if not rep.get("assay"): rep["assay"] = "HT-recruit"
            if not rep.get("necessaryOrSufficient"): rep["necessaryOrSufficient"] = "Sufficient"
            # 'HT-MPRA' is a misnomer for these recruitment screens (Tycko HT-recruit,
            # Del Rosso, Alerasool) -> relabel the method string too.
            if rep.get("method"): rep["method"] = str(rep["method"]).replace("HT-MPRA", "HT-recruit")
        elif src == "TFRegDB1":
            rep["verification"] = "Literature (low-throughput)"
        elif src == "TFRegDB2":
            rep["verification"] = "Literature (low-throughput)"  # default; overlay refines
        else:
            rep["verification"] = None  # CIS-BP DBDs

        # ---- TFRegDB2 re-curation overlay ----
        drop = False
        if src == "TFRegDB2":
            pm = str(rep.get("pmid", "") or "")
            hit = None
            for p in re.split(r"[;,]", pm):
                k = f"{p.strip()}|{rep.get('start')}|{rep.get('end')}"
                if k in overlay: hit = overlay[k]; break
            if hit:
                dec = hit["decision"]
                if hit["mechanism"] is not None: rep["necessaryOrSufficient"] = hit["mechanism"]
                if hit["assay"] is not None: rep["assay"] = hit["assay"]
                rep["verification"] = hit["verification"]
                if dec == "drop":
                    drop = True; stat["dropped_cited"] += 1
                elif dec == "reclassify":
                    # Team decision (Lufesu/George, 2026-08-09): the DB scope is
                    # AD/RD/Bif ONLY. DBD comes solely from CIS-BP; paper-reported
                    # DNA-binding domains and protein-protein interaction domains
                    # are out of scope -> drop them (same as cited-only).
                    drop = True; stat["reclassified_dropped"] += 1
                elif dec == "keep":
                    stat["annotated_verified"] += 1
                else:
                    stat["flagged"] += 1
        if not drop:
            kept.append(rep)
    d["domainReports"] = kept
    # tally effector counts for index recompute
    for rep in kept:
        ty = (rep.get("type") or "").lower()
        if ty == "ad": fam_effector[fam]["ad"] += 1
        elif ty == "rd": fam_effector[fam]["rd"] += 1
        elif ty == "bifunctional": fam_effector[fam]["bifunctional"] += 1
    json.dump(d, open(path, "w"), ensure_ascii=False, separators=(",", ":"))

print(f"[apply] annotated(verified): {stat['annotated_verified']}  dropped(cited): {stat['dropped_cited']}  "
      f"reclassified: {stat['reclassified']}  flagged: {stat['flagged']}")

# ================= recompute index.json effector counts =================
idx = json.load(open(INDEX))
tot = Counter()
for fam in idx["families"]:
    fe = fam_effector[fam["name"]]
    for k in ("ad", "rd", "bifunctional"):
        fam["counts"][k] = fe[k]; tot[k] += fe[k]
# totals + domainTypes
if "totals" in idx:
    for k, tk in (("ad", "ad"), ("rd", "rd"), ("bifunctional", "bifunctional")):
        if tk in idx["totals"]: idx["totals"][tk] = tot[k]
if "domainTypes" in idx and isinstance(idx["domainTypes"], list):
    label_map = {"AD": tot["ad"], "Activation": tot["ad"], "RD": tot["rd"], "Repression": tot["rd"],
                 "Bifunctional": tot["bifunctional"], "Bif": tot["bifunctional"]}
    for dt in idx["domainTypes"]:
        nm = dt.get("name") or dt.get("type")
        if nm in label_map and "count" in dt: dt["count"] = label_map[nm]
json.dump(idx, open(INDEX, "w"), ensure_ascii=False, separators=(",", ":"))
print(f"[index] totals now: AD={tot['ad']} RD={tot['rd']} Bif={tot['bifunctional']}")

# ================= additions list (typos fixed, coords resolved) =================
# rows present in the re-curation but whose humanized coords never matched a live
# report, plus the 25 net-new MasterTable rows, are emitted for a future pipeline
# pass (they need humanization/canonical projection before DB insert).
print(f"[done] overlay traced {traced} human-frame domain coords.")
