"""
Build per-TF mock JSON for tfregdb2 from the v1 source bundle living under
~/Desktop/tfregdb_source/ (symlinks; see README).

What it does, in order, for every TF with at least one curated domain report:

  1. Resolve canonical record from TF_completeIDs (UniProt -> ENSP, gene, sequence).
  2. Aggregate domain reports for that UniProt from Unified_Valid_Domains
     (HumanizedTable + TFRegDB1 + Tycko/DelRosso/Klaus/Arnold/Alerasool).
     Carries assay, activity, confidence, N-or-S, PMID, source database.
  3. Aggregate Pfam DBDs from TF_dbd_data_.
  4. Aggregate per-residue PTMs from PTM_aggregated, filtered to canonical
     translation (Translation_ID match) and clipped to the canonical length.
  5. Copy local AlphaFold PDB into public/pdb/ (skip-if-already-there).
  6. Read local PAE matrix, downsample to <=400x400, write to public/pae/.
  7. Read local ConSurf json, trim to {pos, grade} and write to public/conservation/.
  8. Emit public/mock/tfs/{SYMBOL}.json with the union of everything above.

Then it regenerates public/mock/index.json (home table rollups) from whatever
TF JSONs live under public/mock/tfs/ at the end.

Usage:
  python3 scripts/build_tf_records.py                  # full build (~1.2k TFs)
  python3 scripts/build_tf_records.py --limit 20       # fast iteration
  python3 scripts/build_tf_records.py --symbols TP53,MYC,SP1
  python3 scripts/build_tf_records.py --skip-pae       # skip the slow PAE step
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import re
import shutil
import sys
from typing import Optional

import pandas as pd

# ---------------------------------------------------------------------- paths
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.environ.get("TFREGDB_SOURCE", os.path.expanduser("~/Desktop/tfregdb_source"))

SRC_UNIFIED = os.path.join(SRC, "tables", "Unified_Valid_Domains.xlsx")
SRC_PTMS = os.path.join(SRC, "tables", "PTM_aggregated.xlsx")
# New canonical map is the BLAST-table CSV which carries the `Ensembl canonical`
# column. The old xlsx is kept as the fallback when the CSV isn't present, but
# the gene-name-keyed canonical picker reads from this file first.
SRC_IDS_CSV = os.environ.get(
    "COMPLETE_IDS_CSV",
    os.path.expanduser(
        "~/Desktop/tfregdb/data_pipeline/blast_table/"
        "TF_completeIDs_with_ensembl_canonical.csv"
    ),
)
SRC_IDS = os.path.join(SRC, "tables", "TF_completeIDs.xlsx")
SRC_DBD = os.path.join(SRC, "tables", "TF_dbd_data_.xlsx")
# Per-isoform effector domains with SW-alignment quality (one row per
# (ENSP, effector lineage)). Letters in `Domain_origin` (ORIGINAL_A,
# DERIVADA_A, ALINEAMIENTO_A) tie an alt-isoform projection back to the
# canonical/reference "ORIGINAL" call, so types are resolved by lineage —
# every DERIVADA_X / ALINEAMIENTO_X inherits the type of the ORIGINAL_X
# row for the same gene.
# Per-isoform effector projection CSV. Prefer the locally-rebuilt
# `TF_canonical_to_iso_mapping.csv` (built by scripts/map_canonical_to_isoforms.py
# against the current Unified_Valid_Domains.xlsx) when it exists — it tracks
# the user's curated set after re-binning and uses Ensembl-canonical frames.
# Fall back to the legacy `TF_effector_data_with_derivadas_and_alignment_ideal.csv`
# (precomputed upstream against the older unified) when the new one is absent.
_EFFECTOR_NEW = os.path.join(SRC, "tables", "new", "TF_canonical_to_iso_mapping.csv")
_EFFECTOR_LEGACY = os.path.join(
    SRC, "tables", "new", "TF_effector_data_with_derivadas_and_alignment_ideal.csv"
)
# CHOSEN: legacy (2026-05-29 audit). The rebuilt NEW csv (map_canonical_to_isoforms.py)
# drops 5 canonical TFs (ATF7, CLOCK, PRDM10, ZNF189, ZNF746) + ~701 isoform ENSPs
# and inflates ~816 domain ends by +1, whereas LEGACY's coords are 100% consistent
# with its own domain_seq lengths and its gene↔ENSP labels are clean (0 of 27,430
# rows misattributed). Prefer LEGACY until the NEW rebuild is fixed; fall back to
# NEW only if LEGACY is absent. (Effectors are keyed by ENSP, so gene-label drift
# can't misattach; the _sorted_effectors_for / _pfam_dbds_for guards drop any
# domain that overshoots the isoform length rather than clamping it.)
SRC_EFFECTOR_CSV = _EFFECTOR_LEGACY if os.path.exists(_EFFECTOR_LEGACY) else _EFFECTOR_NEW
SRC_LAMBERT = os.path.join(SRC, "tables", "Lambert2018.txt")
SRC_PDB = os.path.join(SRC, "pdb")
SRC_PAE = os.path.join(SRC, "pae")
SRC_CONS = os.path.join(SRC, "conservation")

DEST_MOCK = os.path.join(ROOT, "public", "mock", "tfs")
DEST_PDB = os.path.join(ROOT, "public", "pdb")
DEST_PAE = os.path.join(ROOT, "public", "pae")
DEST_CONS = os.path.join(ROOT, "public", "conservation")
INDEX_OUT = os.path.join(ROOT, "public", "mock", "index.json")

# One-time parquet cache so xlsx parsing (slow) doesn't dominate iteration runs.
CACHE_DIR = os.path.join(SRC, ".cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def read_excel_cached(xlsx_path: str, sheet=0) -> pd.DataFrame:
    """Read an xlsx, caching as parquet on first call. Reuses cache if xlsx is older."""
    cache = os.path.join(CACHE_DIR, os.path.basename(xlsx_path) + f".{sheet}.parquet")
    if os.path.exists(cache) and os.path.getmtime(cache) >= os.path.getmtime(xlsx_path):
        return pd.read_parquet(cache)
    df = pd.read_excel(xlsx_path, sheet_name=sheet)
    try:
        df.to_parquet(cache, index=False)
    except Exception as e:
        print(f"  (warning: failed to cache {cache}: {e})")
    return df

# Gene-name aliases for the BED's shorthand identifiers (HXA10 → HOXA10 etc.).
# Many of the BED rows use 5-character abbreviations that don't match HGNC
# verbatim; without this map we silently drop ~243 BED gene-symbol prefixes.
BED_GENE_ALIASES: dict[str, str] = {
    # NKX (truncated digit-pairs → hyphenated HGNC)
    "NKX11": "NKX1-1", "NKX12": "NKX1-2",
    "NKX21": "NKX2-1", "NKX22": "NKX2-2", "NKX23": "NKX2-3", "NKX24": "NKX2-4",
    "NKX25": "NKX2-5", "NKX26": "NKX2-6", "NKX28": "NKX2-8",
    "NKX31": "NKX3-1", "NKX32": "NKX3-2",
    "NKX61": "NKX6-1", "NKX62": "NKX6-2", "NKX63": "NKX6-3",
    # L3MBTL (L3 moved to front)
    "LMBL1": "L3MBTL1", "LMBL2": "L3MBTL2",
    "LMBL3": "L3MBTL3", "LMBL4": "L3MBTL4",
    # TWIST
    "TWST1": "TWIST1", "TWST2": "TWIST2",
    # PHOX2
    "PHX2A": "PHOX2A", "PHX2B": "PHOX2B",
    # RHOXF
    "RHXF1": "RHOXF1", "RHXF2": "RHOXF2", "RHXFB": "RHOXF2",
    # MYC / MYB family
    "MYCL1": "MYCL",
    "MYBA": "MYBL1",
    # Singletons
    "RX": "RAX",
    "UNC4": "UNCX",
    "YV1": "YY1",
    "BSH": "BSX",
    "NC2B": "DR1",
    "ARNTL": "BMAL1",     # ARNTL is now BMAL1's primary HGNC symbol on some refs
    "HKR1": "ZNF875",
    "PRDM": "PRDM1",
    "TP53-AD1-2": "TP53",
    "C11ORF95": "ZFTA",   # C11orf95 is renamed ZFTA in current HGNC
}


def normalize_bed_gene(name: str) -> str:
    """Apply the alias map + a few regex rules to map BED's shorthand gene
    symbols to current HGNC symbols. Returns the input upper-cased when no
    rule matches.

    Rules covered:
        *_HUMAN suffix → strip (ALX1_HUMAN → ALX1) — Tycko 2020 convention
        HX[ABCD]N      → HOX[ABCD]N      (HXA10 → HOXA10)
        Z\\d+[A-Z]?    → ZNF\\d+[A-Z]?   (Z324A → ZNF324A)
        explicit map   → curated singletons above
    """
    if not name:
        return name
    n = name.strip().upper()
    # Tycko 2020 / 2024 occasionally use the SwissProt mnemonic-style
    # symbol "<HGNC>_HUMAN" (ALX1_HUMAN, ASCL1_HUMAN, …). The audit found
    # ~190 missing TFs all collapsed to this case. Strip the suffix BEFORE
    # the alias-map lookup so e.g. NKX21_HUMAN → NKX21 → NKX2-1.
    if n.endswith("_HUMAN"):
        n = n[:-len("_HUMAN")]
    if n in BED_GENE_ALIASES:
        return BED_GENE_ALIASES[n]
    # HXA10 → HOXA10
    m = re.match(r"^HX([ABCD])(\d+)$", n)
    if m:
        return f"HOX{m.group(1)}{m.group(2)}"
    # Z324A / Z585B / Z587B → ZNF...
    m = re.match(r"^Z(\d{2,})([A-Z]?)$", n)
    if m:
        return f"ZNF{m.group(1)}{m.group(2)}"
    return n


# ---------------------------------------------------------------------- helpers
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def parse_pdb(uniprot: str) -> tuple[Optional[int], Optional[str]]:
    """Return (length, sequence) from chain A CA atoms. None if PDB missing."""
    path = os.path.join(SRC_PDB, f"AF-{uniprot}-F1.pdb")
    if not os.path.exists(path):
        return None, None
    seq_by_res: dict[int, str] = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[21] != "A":
                continue
            try:
                resi = int(line[22:26])
            except ValueError:
                continue
            seq_by_res[resi] = AA3.get(line[17:20].strip(), "X")
    if not seq_by_res:
        return None, None
    length = max(seq_by_res)
    seq = "".join(seq_by_res.get(i, "X") for i in range(1, length + 1))
    return length, seq


# Map v1 long-form PTM names to our PTMType union. Anything not listed here is
# dropped (rare modifications like hydroxylation / nitration / etc.).
PTM_TYPE_MAP = {
    "Phosphorylation": "Phospho",
    "Dephosphorylation": "Phospho",
    "Acetylation": "Acetyl",
    "Methylation": "Methyl",
    "Ubiquitination": "Ubiquitin",
    "Sumoylation": "Sumo",
    "Neddylation": "Sumo",  # functionally similar; group with Sumo for display
    "O-linked Glycosylation": "Glyco",
    "N-linked Glycosylation": "Glyco",
}

# Map source databases -> short method/assay label for the "Browse by method" facet
SOURCE_METHOD = {
    "TFRegDB1": "TFRegDB1",
    "HumanizedTable": "TFRegDB2",
    "Tycko2020_NucRep": "Tycko 2020 (HT-MPRA)",
    "Tycko2020_NucAct": "Tycko 2020 (HT-MPRA)",
    "Tycko2020_TilingRep": "Tycko 2020 tiling (HT-MPRA)",
    "Tycko2024_AD": "Tycko 2024 (HT-MPRA)",
    "Tycko2024_RD": "Tycko 2024 (HT-MPRA)",
    "DelRosso_AD": "Del Rosso 2023 (HT-MPRA)",
    "DelRosso_RD": "Del Rosso 2023 (HT-MPRA)",
    "Alerasool2022_AD": "Alerasool 2022 (HT-MPRA)",
    "Klaus2023_RD": "Klaus 2023",
    "Arnold2018_AD": "Arnold 2018",
}

# Wrong flagship PMIDs that were once hardcoded here and also leaked into a few
# rows of the curation sheet. Each resolves on PubMed to an unrelated paper
# (breast-cancer margins, tropical-forest degradation, mRNA decoding,
# a triglyceride index, a paediatric injury classification, a psychiatry study),
# so they are never a legitimate citation on a TF effector-domain row.
BAD_PMIDS = {
    "33060797", "38961293", "37020024", "35303881", "37935763", "29483352",
}

# The flagship paper for each high-throughput dataset (real PMIDs, each
# verified against PubMed esummary — author, journal and year all match)
SOURCE_PAPER_PMID = {
    "Tycko2020_NucRep":   "33326746",   # Tycko et al., Cell 2020
    "Tycko2020_NucAct":   "33326746",
    "Tycko2020_TilingRep":"33326746",
    "Tycko2024_AD":       "39487265",   # Tycko et al., Nat Biotechnol 2024
    "Tycko2024_RD":       "39487265",
    "DelRosso_AD":        "37020022",   # DelRosso et al., Nature 2023
    "DelRosso_RD":        "37020022",
    "Alerasool2022_AD":   "35016035",   # Alerasool et al., Mol Cell 2022
    "Klaus2023_RD":       "36545802",   # Klaus et al., EMBO J 2023 (Drosophila RDs)
    "Arnold2018_AD":      "30006452",   # Arnold et al., EMBO J 2018
}

# Curated titles for the HT-MPRA flagship papers (we don't have a Hoja 1 entry
# for them, so we keep a small hardcoded fallback so the table never says
# "PMID 12345" without context).
PAPER_TITLE_FALLBACK = {
    "33326746": "High-Throughput Discovery and Characterization of Human Transcriptional Effectors.",
    "39487265": "Development of compact transcriptional effectors using high-throughput measurements in diverse contexts.",
    "37020022": "Large-scale mapping and mutagenesis of human transcriptional effector domains.",
    "35016035": "Identification and functional characterization of transcriptional activators in human cells.",
    "36545802": "Systematic identification and characterization of repressive domains in Drosophila transcription factors.",
    "30006452": "A high-throughput method to identify trans-activation domains within transcription factor sequences.",
}

# Many source names encode the domain type directly — used when Domain type col is empty
SOURCE_TYPE_HINT = {
    "Tycko2020_NucRep": "RD",
    "Tycko2020_NucAct": "AD",
    "Tycko2020_TilingRep": "RD",
    "Tycko2024_RD": "RD",
    "Tycko2024_AD": "AD",
    "DelRosso_RD": "RD",
    "DelRosso_AD": "AD",
    "Alerasool2022_AD": "AD",
    "Klaus2023_RD": "RD",
    "Arnold2018_AD": "AD",
}

# Map v1 domain types to our DomainType union
DOMAIN_TYPE_MAP = {
    "AD": "AD",
    "RD": "RD",
    "Bif": "Bifunctional",
    "bif": "Bifunctional",
    "Bifunctional": "Bifunctional",
    "bifunctional": "Bifunctional",
    "BIFUNCTIONAL": "Bifunctional",
    "activation": "AD",
    "repression": "RD",
    "Pfam": "DNA",  # for DBDs from TF_dbd_data_
}

# Normalise the family label so v1 variants collapse into one canonical name
FAMILY_NORMALIZE = {
    "p53": "P53", "P53": "P53",
    "HMG/Sox": "SOX / HMG-box", "HMG_box": "SOX / HMG-box", "HMG box": "SOX / HMG-box",
    "Ets": "ETS", "Ets;": "ETS", "ETS family": "ETS",
    "bHLH": "bHLH", "HLH": "bHLH", "Basic helix-loop-helix": "bHLH",
    "bZIP": "bZIP", "Basic leucine zipper": "bZIP", "bZIP_1": "bZIP", "bZIP_2": "bZIP",
    "C2H2 ZF": "C2H2 ZF", "zf-C2H2": "C2H2 ZF", "C2H2 zinc finger": "C2H2 ZF",
    "Homeobox": "Homeodomain", "Homeo": "Homeodomain",
    "Nuclear receptor": "Nuclear receptor", "zf-C4": "Nuclear receptor",
    "Forkhead": "Forkhead", "Fork head": "Forkhead",
    "Unknown": "Unclassified", "Unclassified": "Unclassified",
    "Pfam: UNKNOWN": "Unclassified", "UNKNOWN": "Unclassified",
    "nan": "Unclassified",
    "AT hook": "AT-hook", "AT_hook": "AT-hook",
    "E2F": "E2F", "E2F_TDP": "E2F",
    "Homeodomain; POU": "Homeodomain · POU",
    "Myb/SANT": "Myb/SANT",
    "AP-2": "AP-2",
    "THAP finger": "THAP",
    "CENPB": "CENPB",
}


def normalize_family(name: str) -> str:
    if not name or str(name).lower() in ("nan", "none", "unclassified"):
        return "Unclassified"
    n = str(name).strip()
    # Strip leading "Pfam: " prefix when the rest is a recognised key
    if n.startswith("Pfam: "):
        bare = n[6:]
        if bare in FAMILY_NORMALIZE:
            return FAMILY_NORMALIZE[bare]
        return n
    return FAMILY_NORMALIZE.get(n, n)


# Pfam DBD name -> human-readable TF family
PFAM_FAMILY_MAP = {
    "zf-C2H2": "C2H2 ZF",
    "zf-C4": "Nuclear receptor",
    "Homeobox": "Homeodomain",
    "HLH": "bHLH",
    "bZIP_1": "bZIP",
    "bZIP_2": "bZIP",
    "HMG_box": "SOX / HMG-box",
    "Forkhead": "Forkhead",
    "Ets": "ETS",
    "E2F_TDP": "E2F",
    "P53": "P53",
    "RHD_DNA_bind": "Rel/NFkB",
    "STAT_bind": "STAT",
    "T-box": "T-box",
    "Pou": "POU",
    "Pax": "PAX",
    "GATA": "GATA",
    "zf-MIZ": "zf-MIZ",
    "zf-CXXC": "CxxC ZF",
    "zf-BED": "BED ZF",
    "AT_hook": "AT-hook",
    "KRAB": "KRAB ZF",
    "Runt": "Runt",
    "MH1": "MH1",
    "TF_AP-2": "AP-2",
    "Myb_DNA-binding": "Myb",
    "HSF_DNA-bind": "HSF",
    "SAND": "SAND",
    "SRF-TF": "MADS-box",
    "TBP": "TBP",
    "TF_Otx": "OTX",
    "Tub": "Tubby",
    "Zn_clus": "Zn cluster",
    "GCM": "GCM",
    "CG-1": "CG-1",
    "Nito_FN3_assoc": "Nito",
}

# ---------------------------------------------------------------------- helpers
def _family_from_pfam(dbd_rows: pd.DataFrame) -> str:
    """Pick a family label from the Pfam DBDs annotated for this protein."""
    if dbd_rows is None or dbd_rows.empty:
        return "Unclassified"
    # Most common Pfam wins
    counts = dbd_rows["Pfam"].dropna().value_counts()
    if counts.empty:
        return "Unclassified"
    top = counts.index[0]
    mapped = PFAM_FAMILY_MAP.get(top)
    return mapped if mapped else f"Pfam: {top}"


# ---------------------------------------------------------------------- PAE
def write_pae(uniprot: str, length: int, max_cells: int = 400) -> Optional[int]:
    """Read the local AF PAE matrix, downsample to <=max_cells x max_cells, write JSON.

    Source format: list[{predicted_aligned_error: number[][], ...}]
    Output: {length: int, matrix: number[][]} where matrix is integer-rounded.
    """
    src = os.path.join(SRC_PAE, f"AF-{uniprot}-F1-pae.json")
    if not os.path.exists(src):
        return None
    with open(src) as f:
        raw = json.load(f)
    entry = raw[0] if isinstance(raw, list) else raw
    mat = entry.get("predicted_aligned_error")
    if not mat:
        return None

    n = len(mat)
    if n != length:
        # Trust the matrix; native plot labels use length, mapping handled in viewer
        length = n

    if n > max_cells:
        stride = n / max_cells
        sample = [int(i * stride) for i in range(max_cells)]
    else:
        sample = list(range(n))

    out = [[round(mat[i][j], 1) for j in sample] for i in sample]
    dest = os.path.join(DEST_PAE, f"{uniprot}.json")
    with open(dest, "w") as f:
        json.dump({"length": length, "matrix": out}, f)
    return os.path.getsize(dest)


# ---------------------------------------------------------------------- conservation
def write_conservation(uniprot: str) -> Optional[int]:
    """Trim ConSurf JSON to {pos, grade} per residue and write."""
    src = os.path.join(SRC_CONS, f"{uniprot}.json")
    if not os.path.exists(src):
        return None
    try:
        with open(src) as f:
            raw = json.load(f)
    except Exception:
        return None
    out = [
        {"pos": int(r.get("pos", i + 1)), "grade": int(r.get("color", 0)), "aa": r.get("aa", "X")}
        for i, r in enumerate(raw)
        if r.get("color") is not None
    ]
    dest = os.path.join(DEST_CONS, f"{uniprot}.json")
    with open(dest, "w") as f:
        json.dump(out, f)
    return os.path.getsize(dest)


# ---------------------------------------------------------------------- main builder
def build(args):
    # Start the per-TF JSON dir clean so a shrinking TF set (e.g. after a
    # source-table swap) can't leave stale JSONs behind — the index step
    # (which globs this dir) would otherwise pick them up.
    if os.path.isdir(DEST_MOCK):
        shutil.rmtree(DEST_MOCK)
    os.makedirs(DEST_MOCK, exist_ok=True)
    os.makedirs(DEST_PDB, exist_ok=True)
    os.makedirs(DEST_PAE, exist_ok=True)
    os.makedirs(DEST_CONS, exist_ok=True)

    print(f"[1/6] Loading canonical IDs + Lambert2018 family map...")
    # Prefer the new CSV with the `Ensembl canonical` column when available;
    # otherwise fall back to the legacy xlsx. The new file is the source of
    # truth for which ENSP is canonical (column literal: "Canonical").
    if os.path.exists(SRC_IDS_CSV):
        ids_df = pd.read_csv(SRC_IDS_CSV)
        print(f"  source: {SRC_IDS_CSV}  rows={len(ids_df)}")
        if "Ensembl canonical" in ids_df.columns:
            print(
                f"  Ensembl-canonical-tagged rows: "
                f"{(ids_df['Ensembl canonical'] == 'Canonical').sum()}"
            )
    else:
        ids_df = read_excel_cached(SRC_IDS)
        print(f"  source: {SRC_IDS}  rows={len(ids_df)}  (legacy fallback)")
    # Lambert 2018 (DatabaseExtract_v_1.01.txt) — the gold standard human TF
    # roster. We use its 'DBD' column verbatim as the family label, keyed by
    # HGNC symbol. Falls through to Pfam/v1 family for genes not in Lambert.
    try:
        lambert = pd.read_csv(SRC_LAMBERT, sep="\t")
        lambert_yes = lambert[lambert["Is TF?"] == "Yes"]
        lambert_family: dict[str, str] = {
            str(r["HGNC symbol"]).strip().upper(): str(r["DBD"]).strip()
            for _, r in lambert_yes.iterrows()
            if pd.notna(r["HGNC symbol"]) and pd.notna(r["DBD"])
        }
        # Pull the human-readable protein name from the EntrezGene Description
        # column ("tumor protein p53 [Source:HGNC Symbol;Acc:HGNC:11998]" → strip
        # the trailing source tag). This drives the full-name shown next to the
        # symbol in the heatmap rows.
        def _clean_descr(s: str) -> str:
            s = str(s).split("[Source:")[0].strip()
            # Lambert descriptions are lowercase ("tumor protein p53"). Capitalize
            # the first character so the heatmap renders "Tumor protein p53".
            return s[:1].upper() + s[1:] if s else s
        lambert_name: dict[str, str] = {
            str(r["HGNC symbol"]).strip().upper(): _clean_descr(r["EntrezGene Description"])
            for _, r in lambert_yes.iterrows()
            if pd.notna(r["HGNC symbol"]) and pd.notna(r["EntrezGene Description"])
        }
        # Build per-family universe sizes by normalising each Lambert DBD label
        # through the same FAMILY_NORMALIZE map we use for our TF families. This
        # is the canonical denominator for "curated / total" — using a hardcoded
        # hint table caused E2F to show 9/8 when Lambert actually has 11 E2F TFs
        # (E2F1–E2F8 + TFDP1/2/3).
        lambert_universe: dict[str, int] = collections.Counter()
        for _, r in lambert_yes.iterrows():
            if pd.isna(r["DBD"]):
                continue
            fam = normalize_family(str(r["DBD"]).strip())
            lambert_universe[fam] += 1
        print(
            f"  Lambert 2018 TF families indexed: {len(lambert_family)} "
            f"({len(lambert_name)} with names, {len(lambert_universe)} family buckets)"
        )
    except Exception as e:
        lambert_family = {}
        lambert_name = {}
        lambert_universe = {}
        print(f"  (warning: Lambert load failed: {e})")
    # Canonicalise the column shapes the downstream code reads.
    if "UniProt_ID" in ids_df.columns:
        ids_df["UniProt_Clean"] = (
            ids_df["UniProt_ID"].astype(str).str.split("-").str[0].str.strip()
        )
    else:
        ids_df["UniProt_Clean"] = ""
    ids_df["Translation_Clean"] = (
        ids_df["Translation ID"].astype(str).str.split(".").str[0].str.strip()
    )
    ids_df["Gene_Upper"] = (
        ids_df["Gene name (HGNC)"].astype(str).str.strip().str.upper()
    )

    # Canonical-row picker: group rows by Gene Name (HGNC) — NOT by UniProt,
    # which used to lump alt isoforms together and pick arbitrarily within the
    # group. Inside each gene's rows, pick the canonical row via a tiered
    # preference chain, with the new CSV's `Ensembl canonical = "Canonical"`
    # flag at the top:
    #   1. Row tagged `Ensembl canonical == "Canonical"`.
    #   2. APPRIS PRINCIPAL:1.
    #   3. Any APPRIS PRINCIPAL tier.
    #   4. Longest protein in the group.
    #   5. First row (last resort).
    canon_records: dict[str, pd.Series] = {}  # gene_upper -> chosen row
    canon_by_ensp: dict[str, pd.Series] = {}  # ENSP -> chosen row (for direct lookups)
    has_canon_col = "Ensembl canonical" in ids_df.columns
    for gene, group in ids_df.dropna(subset=["Translation_Clean"]).groupby("Gene_Upper"):
        if not gene or gene == "NAN":
            continue
        # Drop rows that don't have a usable ENSP.
        usable = group[group["Translation_Clean"].str.startswith("ENSP", na=False)]
        if usable.empty:
            continue
        pick = None
        if has_canon_col:
            ens_canon = usable[usable["Ensembl canonical"].astype(str) == "Canonical"]
            if len(ens_canon):
                pick = ens_canon.iloc[0]
        if pick is None:
            p1 = usable[
                usable["APPRIS Annotation"].astype(str).str.contains("PRINCIPAL:1", na=False)
            ]
            if len(p1):
                pick = p1.iloc[0]
        if pick is None:
            pany = usable[
                usable["APPRIS Annotation"].astype(str).str.contains("PRINCIPAL", na=False)
            ]
            if len(pany):
                pick = pany.iloc[0]
        if pick is None:
            with_len = usable.dropna(subset=["Protein length"]) if "Protein length" in usable.columns else usable
            if len(with_len):
                pick = with_len.sort_values("Protein length", ascending=False).iloc[0]
        if pick is None:
            pick = usable.iloc[0]
        canon_records[gene] = pick
        canon_by_ensp[str(pick["Translation_Clean"])] = pick
    print(f"  canonical records by Gene Name: {len(canon_records)}")

    # ---- Reviewed (Swiss-Prot) DISPLAY accession per gene (Route 2) --------
    # The canonical picker keys the page + every asset (AlphaFold model, PAE,
    # conservation) and the domain coordinate frame on the Ensembl-canonical
    # isoform's accession, which for ~30 TFs is a TrEMBL (A0A...) entry. We do
    # NOT move that key (swapping to a reviewed isoform of a different length
    # would misplace curated domains). Instead we resolve a reviewed accession
    # used ONLY for the public UniProt cross-link / identity, leaving structure
    # and coordinates untouched.
    _SWISSPROT_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$")  # high-confidence reviewed
    # Reviewed accessions for the handful of genes whose reviewed entry is not
    # present anywhere in TF_completeIDs. VERIFY against UniProt before release.
    # All five verified against the UniProt REST API (reviewed, human), 2026-08-02.
    REVIEWED_OVERRIDE = {
        "POU6F1": "Q14863", "AHRR": "A9YTQ3", "POU3F4": "P49335",
        "TPRX1": "Q8N7U7", "ZBTB34": "Q8NCN2",
    }
    reviewed_acc_by_gene: dict[str, str] = {}
    for gene, grp in ids_df.groupby("Gene_Upper"):
        # Among the gene's reviewed (Swiss-Prot) accessions, prefer the one on
        # the LONGEST protein, so a short secondary fragment isn't chosen over
        # the real reviewed entry (e.g. GLI2 Q1PSW9 82 aa vs P10070 1586 aa).
        cand: dict[str, int] = {}
        for _, r in grp.iterrows():
            a = str(r.get("UniProt_Clean") or "").split("-")[0].strip().upper()
            if a and _SWISSPROT_RE.match(a):
                try:
                    plen = int(float(r.get("Protein length")))
                except (TypeError, ValueError):
                    plen = 0
                cand[a] = max(cand.get(a, 0), plen)
        if cand:
            reviewed_acc_by_gene[gene] = max(cand, key=lambda a: cand[a])
    for g, acc in REVIEWED_OVERRIDE.items():
        reviewed_acc_by_gene.setdefault(g.upper(), acc)

    def _reviewed_accession(gene: str, canonical_up: str) -> str:
        """Reviewed accession for display; falls back to the canonical one."""
        up = str(canonical_up).split("-")[0].strip().upper()
        if _SWISSPROT_RE.match(up):
            return up  # canonical is already a reviewed Swiss-Prot accession
        return reviewed_acc_by_gene.get(gene, up)

    print(f"[2/6] Loading Unified_Valid_Domains.xlsx + Hoja 1 paper titles...")
    unified = read_excel_cached(SRC_UNIFIED)
    # Hoja 1 has PMID + Title for ~3137 unique papers — the literature
    # backbone that HumanizedTable was curated from. Use it to enrich the
    # per-TF papers list with real titles.
    raw_xlsx = os.environ.get(
        "RAW_XLSX", os.path.expanduser("~/Desktop/tfregdb/dev/UpdatedTable.xlsx")
    )
    pmid_titles: dict[str, str] = {}
    try:
        # We cache by xlsx mtime, just like the other tables.
        hoja1 = pd.read_excel(raw_xlsx, sheet_name="Hoja 1")
        hoja1["PMID_Clean"] = (
            hoja1["PMID"].astype(str).str.split(".").str[0].str.strip()
        )
        for _, r in hoja1.dropna(subset=["PMID", "Title"]).iterrows():
            pmid = str(r["PMID_Clean"])
            if pmid and pmid not in pmid_titles:
                pmid_titles[pmid] = str(r["Title"]).strip()
        print(f"  Hoja 1 paper titles indexed: {len(pmid_titles)}")
    except Exception as e:
        print(f"  (warning: Hoja 1 title load failed: {e})")
    # Filter to rows with a clean human UniProt mapping
    unified = unified.dropna(subset=["Final_Human_UNIPROT_ID"]).copy()
    unified["UP_Clean"] = (
        unified["Final_Human_UNIPROT_ID"].astype(str).str.split("-").str[0].str.strip()
    )
    unified["Gene_Upper"] = (
        unified["Original_Gene"]
        .astype(str)
        .str.strip()
        .str.upper()
        .map(normalize_bed_gene)
    )
    # Re-bin every row under the canonical HGNC gene that owns its
    # Final_Human_UNIPROT_ID. We trust the BLAST projection's UniProt over the
    # curator's free-form Original_Gene, because:
    #   • aliases (ZN274 / BMAL1 / NRF1) collapse to canonical (ZNF274 / ARNTL
    #     / NFE2L1) automatically;
    #   • curator typos that look superficially valid (Original_Gene=TEF but
    #     UniProt=P28347→TEAD1, Original_Gene=TWIST1 but UniProt→TWIST2) get
    #     moved to the right TF page instead of polluting the wrong one.
    # Only fall back to Original_Gene when the UniProt isn't in any canonical
    # record (rare — usually means the UniProt is a non-human accession that
    # slipped through the BLAST filter).
    up_to_gene: dict[str, str] = {}
    for g, row in canon_records.items():
        up = str(row.get("UniProt_Clean") or "").strip().upper()
        if up and up not in up_to_gene:
            up_to_gene[up] = g
    canon_gene_set = set(canon_records.keys())

    def _canon_seq(g: str) -> str:
        r = canon_records.get(g)
        if r is None:
            return ""
        s = str(r.get("Sequence aa") or "")
        return s.upper() if s and s != "nan" else ""

    def _dom_seq(row: pd.Series) -> str:
        for c in ("Final_Human_Domain_Sequence", "Domain_Sequence"):
            v = row.get(c)
            if pd.notna(v) and str(v).strip() and str(v) != "nan":
                return re.sub(r"[^A-Za-z]", "", str(v)).upper()
        return ""

    def _resolve_gene(row: pd.Series) -> str:
        # Re-bin a curated row to the canonical gene that owns its BLAST
        # Final_Human_UNIPROT_ID — but SEQUENCE-AWARE, not blindly. The BLAST
        # projection sometimes lands a curated domain on a paralog/readthrough
        # (e.g. CLOCK->BMAL1, RORC->RORA, ATF7->ATF7-NPFF), which used to steal
        # the whole TF's rows and drop it from the site.
        up = str(row["UP_Clean"] or "").strip().upper()
        og = row["Gene_Upper"]
        owner = up_to_gene.get(up)
        if owner is None:
            return og                       # UniProt owned by no canon gene -> keep curator label
        if owner == og:
            return owner
        # Rebin would move to a DIFFERENT canonical gene. Trust the BLAST owner
        # only when Original_Gene is not itself a canonical TF (real alias/typo
        # fix, e.g. ZN274 -> ZNF274). When Original_Gene IS canonical, the
        # curator named a real TF; keep it UNLESS the curated domain sequence
        # proves the label wrong (present in the owner's protein but absent from
        # the curator's — e.g. mislabeled EBF3 tiles that are really EBF2 seq).
        if og not in canon_gene_set:
            return owner
        dseq = _dom_seq(row)
        if dseq:
            if dseq in _canon_seq(owner) and dseq not in _canon_seq(og):
                return owner                # curator label is demonstrably wrong
        return og                           # trust the curator's canonical gene
    n_orig_canon = (unified["Gene_Upper"].isin(canon_gene_set)).sum()
    new_genes = unified.apply(_resolve_gene, axis=1)
    n_moved = (new_genes != unified["Gene_Upper"]).sum()
    unified["Gene_Upper"] = new_genes
    n_after = (unified["Gene_Upper"].isin(canon_gene_set)).sum()
    if n_moved:
        print(
            f"  Re-binned {n_moved} unified rows to the canonical-gene owner"
            f" of their Final_Human_UNIPROT_ID"
            f" (Original_Gene canonical-hit rate: {n_orig_canon} → {n_after})"
        )
    # Coalesce PMID across the two source columns — `Reference (PMID)` (TFRegDB1)
    # and `PUBMED_ID` (HumanizedTable). pd.fillna handles real NaN; the .isin
    # filter cleans up stringified-NaN that bled through prior type coercion.
    unified["PMID_Clean"] = (
        unified["Reference (PMID)"]
        .fillna(unified["PUBMED_ID"])
        .astype(str)
        .str.split(".")
        .str[0]
        .str.strip()
    )
    unified.loc[
        unified["PMID_Clean"].str.lower().isin(["nan", "none", ""]), "PMID_Clean"
    ] = ""
    # Coalesce domain type across columns + source-database fallback
    def _coalesce_type(row):
        for col in ("Domain type", "DomainType"):
            v = row.get(col)
            if pd.notna(v) and str(v).strip() and str(v) != "nan":
                return str(v).strip()
        src = str(row.get("Source_Database", ""))
        return SOURCE_TYPE_HINT.get(src, "")
    unified["DomainType_Clean"] = unified.apply(_coalesce_type, axis=1)
    print(f"  domain reports: {len(unified)}")
    print(f"  unique UniProts: {unified['UP_Clean'].nunique()}, unique genes: {unified['Gene_Upper'].nunique()}")

    print(f"[3/6] Loading PTM_aggregated.xlsx... (this takes a moment)")
    ptms_df = read_excel_cached(SRC_PTMS)
    ptms_df["UP_Clean"] = ptms_df["UniProt"].astype(str).str.split("-").str[0].str.strip()
    ptms_df["Translation_Clean"] = (
        ptms_df["Translation ID"].astype(str).str.split(".").str[0].str.strip()
    )
    ptms_df["PMID"] = ptms_df["PMID"].astype(str).str.split(";").str[0].str.strip()
    # Index by canonical Translation_ID for fast filter
    ptms_by_trans: dict[str, pd.DataFrame] = dict(tuple(ptms_df.groupby("Translation_Clean")))
    print(f"  PTM rows: {len(ptms_df)} -> per-Translation buckets: {len(ptms_by_trans)}")

    print(f"[4/6] Loading TF_dbd_data_.xlsx...")
    dbd_df = read_excel_cached(SRC_DBD)
    dbd_df = dbd_df.dropna(subset=["Translation ID", "Pfam", "From", "To"]).copy()
    dbd_df["Translation_Clean"] = (
        dbd_df["Translation ID"].astype(str).str.split(".").str[0].str.strip()
    )
    dbd_by_trans: dict[str, pd.DataFrame] = dict(tuple(dbd_df.groupby("Translation_Clean")))
    print(f"  DBD rows: {len(dbd_df)}")

    # -------- Per-isoform effector domains (AD / RD / Bifunctional) --------
    # The new csv has one row per (ENSP × effector lineage) with the iso's own
    # From/To plus SW-alignment quality for projected calls. Resolve type via
    # (gene, source_db, source_row) lookup back into unified, where the curated
    # type lives in `Domain type` (TFRegDB1) or `DomainType` (HumanizedTable),
    # with a SOURCE_TYPE_HINT fallback for the experimental MPRA sources whose
    # database name already encodes _AD_ / _RD_.
    print(f"[4b/6] Loading per-isoform effectors from new csv...")
    effector_by_ensp: dict[str, list[dict]] = {}
    if os.path.exists(SRC_EFFECTOR_CSV):
        eff_df = pd.read_csv(SRC_EFFECTOR_CSV)
        eff_df["ensp_clean"] = (
            eff_df["Translation ID"].astype(str).str.split(".").str[0].str.strip()
        )
        eff_df["gene_upper"] = eff_df["Gene name (HGNC)"].astype(str).str.strip().str.upper()
        # (gene, source_db_prefix, row_idx) -> type via coalesced unified columns.
        # `Source_Database` in unified can be "DelRosso_AD"; we keep the prefix
        # before the first underscore since the new csv's Source string uses
        # "<GENE>_<DB_PREFIX>_..._RowN" (DB prefix is everything between the
        # gene token and the trailing "_Row<N>").
        unified_type_map: dict[tuple[str, str, int], str] = {}
        for _, ur in unified.iterrows():
            g = str(ur.get("Gene_Upper") or "").strip()
            sd = str(ur.get("Source_Database") or "").strip()
            ri = ur.get("Original_Row_Index")
            if not g or g == "NAN" or not sd or pd.isna(ri):
                continue
            t_clean = str(ur.get("DomainType_Clean") or "").strip()
            if t_clean:
                mapped = DOMAIN_TYPE_MAP.get(t_clean, "Other")
                if mapped == "Other":
                    tl = t_clean.lower()
                    if "bifunc" in tl: mapped = "Bifunctional"
                    elif "activ" in tl: mapped = "AD"
                    elif "repress" in tl: mapped = "RD"
            else:
                mapped = SOURCE_TYPE_HINT.get(sd, "")
            if not mapped or mapped == "Other":
                continue
            try:
                key = (g, sd.split("_")[0], int(ri))
            except (ValueError, TypeError):
                continue
            unified_type_map.setdefault(key, mapped)

        _src_re = re.compile(r"^\s*([^_]+)_(.+?)_Row(\d+)\s*$")

        def _resolve_eff_type(gene: str, source_str: str) -> str | None:
            # Fast path: explicit _AD_/_RD_ in the source database token.
            s = str(source_str)
            if "_AD_" in s or "_NucAct" in s or "TilingAct" in s:
                return "AD"
            if "_RD_" in s or "_NucRep" in s or "TilingRep" in s:
                return "RD"
            # Cross-reference into unified by (gene, db_prefix, row).
            for piece in s.split("|"):
                m = _src_re.match(piece.strip())
                if not m:
                    continue
                _g, db_full, row = m.group(1).upper(), m.group(2), int(m.group(3))
                db_prefix = db_full.split("_")[0]
                if (gene, db_prefix, row) in unified_type_map:
                    return unified_type_map[(gene, db_prefix, row)]
            return None

        eff_df["origin_kind"] = (
            eff_df["Domain_origin"].astype(str).str.split("_").str[0]
        )
        eff_df["lineage"] = (
            eff_df["Domain_origin"].astype(str).str.extract(r"_([A-Z])$")[0]
        )

        # Coord-based unified type map for the rebuilt BLAST pipeline:
        # the effector csv's Source strings carry legacy row indices into the
        # old unified xlsx (e.g. RARA_TFRegDB1_Row596), but the rebuilt unified
        # uses different row indices, so the (gene, db_prefix, row) lookup
        # misses. The effector csv's `From`/`To` coords still match unified
        # rows exactly, so we key by (gene_upper, from, to) -> mapped type.
        unified_coord_type_map: dict[tuple[str, int, int], str] = {}
        _coord_re = re.compile(r"\s*(\d+)\s*-\s*(\d+)\s*")
        for _, ur in unified.iterrows():
            g = str(ur.get("Gene_Upper") or "").strip()
            if not g or g == "NAN":
                continue
            t_clean = str(ur.get("DomainType_Clean") or "").strip()
            if not t_clean:
                continue
            mapped = DOMAIN_TYPE_MAP.get(t_clean, "Other")
            if mapped == "Other":
                tl = t_clean.lower()
                if "bifunc" in tl: mapped = "Bifunctional"
                elif "activ" in tl: mapped = "AD"
                elif "repress" in tl: mapped = "RD"
            if not mapped or mapped == "Other":
                continue
            for coord_col in ("DomainCoordinates", "Coordinates"):
                c = ur.get(coord_col)
                if pd.isna(c):
                    continue
                raw = str(c).replace("–", "-").replace("—", "-").lstrip(":")
                m = _coord_re.match(raw)
                if m:
                    s_, e_ = int(m.group(1)), int(m.group(2))
                    unified_coord_type_map.setdefault((g, s_, e_), mapped)
                    break

        # Step 1: per (gene, lineage) resolve type once, preferring an ORIGINAL row.
        gene_lin_type: dict[tuple[str, str], str] = {}
        eff_clean = eff_df.dropna(subset=["lineage", "Source", "gene_upper", "ensp_clean"])
        for (gene, lin), grp in eff_clean.groupby(["gene_upper", "lineage"]):
            anchor = grp[grp["origin_kind"] == "ORIGINAL"]
            anchor = anchor.iloc[0] if len(anchor) else grp.iloc[0]
            # Coord-based lookup first (works for the rebuilt unified); fall
            # back to the legacy source-string parsing if coords don't match.
            t = None
            try:
                a_from = int(float(anchor["From"]))
                a_to = int(float(anchor["To"]))
                t = unified_coord_type_map.get((str(gene), a_from, a_to))
            except (TypeError, ValueError):
                pass
            if not t:
                t = _resolve_eff_type(str(gene), str(anchor["Source"]))
            if t:
                gene_lin_type[(str(gene), str(lin))] = t

        # Step 2: walk every row, build the per-ENSP effector list.
        unresolved = 0
        for _, r in eff_clean.iterrows():
            gene = str(r["gene_upper"])
            lin = str(r["lineage"])
            t = gene_lin_type.get((gene, lin))
            if not t:
                # Last-ditch fallback: try the row's own Source.
                t = _resolve_eff_type(gene, str(r["Source"]))
            if not t:
                unresolved += 1
                continue
            try:
                s_ = int(float(r["From"]))
                e_ = int(float(r["To"]))
            except (TypeError, ValueError):
                continue
            if e_ < s_:
                continue
            origin = str(r["origin_kind"]).lower()
            entry: dict = {
                "type": t,
                "name": t,
                "start": s_,
                "end": e_,
                "origin": origin,  # "original" | "derivada" | "alieneamiento"
                "lineage": lin,
                "source": str(r["Source"]),
            }
            if origin == "alieneamiento":
                ai = r.get("alignment_identity")
                ac = r.get("alignment_coverage")
                # A sequence-aligned transfer is accepted as present only at
                # full (100%) identity and full (100%) coverage of the domain;
                # otherwise the isoform is treated as not carrying it.
                if not (
                    pd.notna(ai)
                    and pd.notna(ac)
                    and float(ai) >= 99.999
                    and float(ac) >= 99.999
                ):
                    continue
                entry["identity"] = round(float(ai), 1)
                entry["coverage"] = round(float(ac), 1)
            effector_by_ensp.setdefault(str(r["ensp_clean"]), []).append(entry)
        print(
            f"  effector csv rows: {len(eff_df)} (per-ENSP buckets: {len(effector_by_ensp)}, "
            f"lineages typed: {len(gene_lin_type)}, unresolved: {unresolved})"
        )
    else:
        print(f"  (no per-isoform effector csv at {SRC_EFFECTOR_CSV}; skipping)")

    # The gene -> canonical map is already built upstream as `canon_records`
    # (with the Ensembl-canonical priority chain). Reuse it; only require that
    # the chosen row also has a Sequence aa entry — without one we can't emit
    # the per-TF JSON, so skip those genes.
    gene_to_canon: dict[str, pd.Series] = {
        g: row
        for g, row in canon_records.items()
        if "Sequence aa" in row.index and pd.notna(row["Sequence aa"])
    }
    print(f"  gene -> canonical map (with sequence): {len(gene_to_canon)} genes")

    # -------- Determine TF set: group reports by gene
    genes_with_reports = sorted(unified["Gene_Upper"].dropna().unique())
    # Drop genes with no canonical record (unusable)
    eligible_genes = [g for g in genes_with_reports if g in gene_to_canon]
    print(f"  genes with reports: {len(genes_with_reports)}, with canonical record: {len(eligible_genes)}")

    if args.symbols:
        wanted = {s.strip().upper() for s in args.symbols.split(",")}
        eligible_genes = [g for g in eligible_genes if g in wanted]
        print(f"  filtered to {len(eligible_genes)} TFs from --symbols")
    if args.limit:
        eligible_genes = eligible_genes[: args.limit]
        print(f"  limited to first {len(eligible_genes)} TFs")

    print(f"[5/6] Building per-TF JSON for {len(eligible_genes)} TFs...")

    written: list[str] = []
    skipped_no_pdb: list[str] = []
    skipped_no_canon: list[str] = []
    pae_bytes_total = 0
    cons_bytes_total = 0
    pdb_bytes_total = 0
    SOURCE_SET = set()
    # Rows in Unified_Valid_Domains whose Final_Human_DomainCoords end > the
    # canonical protein length. We clip them so the build still runs, but log
    # every drop so the curator can fix the xlsx (commonly: gene-name vs
    # UniProt mismatch — e.g. an "Original_Gene = TEF" row whose
    # Final_Human_UNIPROT_ID actually points to TEAD1's P28347).
    OUT_OF_BOUNDS_LOG: list[dict] = []

    # Pre-group unified by gene for fast lookup
    unified_by_gene = dict(tuple(unified.groupby("Gene_Upper")))

    for idx, gene in enumerate(eligible_genes, 1):
        canon = gene_to_canon.get(gene)
        if canon is None:
            skipped_no_canon.append(gene)
            continue
        symbol = gene
        uniprot = str(canon["UniProt_Clean"]).strip()

        # The canonical sequence MUST come from TF_completeIDs (Ensembl
        # canonical) — not from AlphaFold's PDB, which tracks UniProt's
        # canonical isoform and can disagree. E.g. CIC: Ensembl canonical is
        # ENSP00000505728 = 2517 aa, but UniProt canonical Q96RK0-1 is 1608 aa
        # so AF-Q96RK0.pdb is only 1608 aa. The curated coords (2410-2489 etc.)
        # are in Ensembl-canonical frame and only fit the 2517-aa version.
        # Only fall back to PDB if TF_completeIDs has no sequence at all.
        seq_from_ids = str(canon.get("Sequence aa", "") or "").strip()
        # A valid canonical sequence is all standard amino-acid letters. Guard
        # against upstream placeholder / junk cells: e.g. NOBOX's Ensembl-
        # canonical row carries the literal string "Not found" (9 chars), which
        # would otherwise be accepted as a 9-residue protein and drop every
        # curated domain. Also reject a sequence far shorter than the row's own
        # `Protein length` (a truncated extraction). In either case fall through
        # to the full-length AlphaFold sequence for the (reviewed) accession.
        _AA = set("ACDEFGHIKLMNPQRSTVWYXUBZ")
        try:
            row_len = int(float(canon.get("Protein length")))
        except (TypeError, ValueError):
            row_len = None
        seq_ok = (
            bool(seq_from_ids)
            and seq_from_ids.lower() != "nan"
            and set(seq_from_ids.upper()) <= _AA
            and (row_len is None or len(seq_from_ids) >= 0.5 * row_len)
        )
        if seq_ok:
            sequence = seq_from_ids
            length = len(sequence)
        else:
            if seq_from_ids and seq_from_ids.lower() != "nan":
                print(
                    f"  [seq-guard] {symbol}: rejected bad 'Sequence aa' "
                    f"({seq_from_ids[:20]!r}, len {len(seq_from_ids)}"
                    f"{f', row Protein length {row_len}' if row_len else ''}); "
                    f"using AlphaFold sequence for {uniprot}"
                )
            length, sequence = parse_pdb(uniprot)
            if length is None:
                skipped_no_pdb.append(symbol)
                continue

        ensg = str(canon.get("Ensembl Gene ID", "") or "")
        ensp = str(canon.get("Translation_Clean", "") or "")

        # ---- Domain reports for this gene (across all UniProt variants)
        tf_rows = unified_by_gene.get(gene, pd.DataFrame())
        reports = []
        used_pmids: set[str] = set()
        for _, r in tf_rows.iterrows():
            coords = str(r.get("Final_Human_DomainCoords", "") or "").strip()
            if "-" not in coords:
                # A handful of curated effector rows (HOXB2, FOXO1, FOXO4, SOX4,
                # SOX12, HOXA2, E2F3, E2F5, PLAGL2, ZBED6) never got
                # Final_Human_DomainCoords populated even though the BLAST step
                # aligned the domain to the human canonical (BLAST_DomainCoordinates
                # + Final_Identity are set). Fall back to the BLAST coords so the
                # domain isn't silently dropped from the protein.
                coords = str(r.get("BLAST_DomainCoordinates", "") or "").strip()
            if "-" not in coords:
                continue
            try:
                a, b = coords.split("-")
                start, end = int(a), int(b)
            except ValueError:
                continue
            if start < 1:
                start = 1
            if end > length:
                # The row's coords don't fit the canonical. Either:
                #   • boundary noise (≤20 aa over — curator off-by-a-few or a
                #     near-canonical isoform). We clip and keep.
                #   • frame mismatch (>20 aa over — curator referenced a much
                #     longer isoform like XBP1s 376aa vs canonical XBP1u 261aa,
                #     or CIC long isoform 2517aa vs canonical 1608aa). We drop
                #     so the page doesn't show a misleadingly-truncated chip;
                #     the audit log carries the original coords for review.
                overflow = end - length
                drop = overflow > 20
                OUT_OF_BOUNDS_LOG.append({
                    "symbol": symbol,
                    "length": length,
                    "uniprot": uniprot,
                    "coords": coords,
                    "end_over_by": overflow,
                    "action": "dropped" if drop else "clipped",
                    "source": str(r.get("Source_Database") or ""),
                    "pmid": str(r.get("PMID_Clean") or ""),
                    "type": str(r.get("DomainType_Clean") or ""),
                    "row_id": str(r.get("Effector domain ID") or ""),
                })
                if drop:
                    continue
                end = length
            if end < start:
                continue
            dtype = str(r.get("DomainType_Clean", "")).strip()
            mapped = DOMAIN_TYPE_MAP.get(dtype, "Other")
            if mapped == "Other" and dtype:
                # Some literature wrote e.g. "Repression" / "Activation" verbatim
                dl = dtype.lower()
                if "bifunc" in dl: mapped = "Bifunctional"
                elif "activ" in dl: mapped = "AD"
                elif "repress" in dl: mapped = "RD"
            source = str(r.get("Source_Database", "")).strip()
            # Rename the legacy curation source name in flight so per-TF JSONs
            # ship the user-facing "TFRegDB2" instead of the internal
            # data-engine name "HumanizedTable".
            if source == "HumanizedTable":
                source = "TFRegDB2"
            pmid = str(r.get("PMID_Clean", "")).strip()
            if pmid in ("", "nan", "None"):
                pmid = SOURCE_PAPER_PMID.get(source, "")
            elif pmid in BAD_PMIDS and source in SOURCE_PAPER_PMID:
                # A handful of curated rows carry the same wrong flagship PMIDs
                # this file used to hardcode (they resolve to unrelated papers —
                # see BAD_PMIDS). Trust the dataset's flagship PMID over the
                # sheet for those rows only.
                pmid = SOURCE_PAPER_PMID[source]
            assay = str(r.get("Assay") or "").strip()
            if assay.lower() in ("nan", "none", ""): assay = ""
            activity = str(r.get("Activity (H, M or L)") or "").strip()
            if activity.lower() in ("nan", "none", ""): activity = ""
            confidence = str(r.get("Confidence (H, M or L)") or "").strip()
            if confidence.lower() in ("nan", "none", ""): confidence = ""
            n_or_s = str(r.get("N or S?") or "").strip()
            if n_or_s.lower() in ("nan", "none", ""): n_or_s = ""
            identity = r.get("Final_Identity")
            note_parts = []
            for col in ("Notes", "Note"):
                v = r.get(col)
                if pd.notna(v) and str(v).strip() not in ("nan", ""):
                    note_parts.append(str(v).strip())
            method = SOURCE_METHOD.get(source, source)
            SOURCE_SET.add(method)
            # Choose a human-readable name for the bar label
            if mapped == "DNA":
                domain_label = "DBD"
            elif mapped == "AD":
                domain_label = "AD"
            elif mapped == "RD":
                domain_label = "RD"
            elif mapped == "Bifunctional":
                domain_label = "Bif"
            else:
                domain_label = dtype or "Region"
            reports.append({
                "name": domain_label,
                "type": mapped,
                "start": start,
                "end": end,
                "pmid": pmid,
                "source": source,
                "method": method,
                "assay": assay or None,
                "activity": activity or None,
                "confidence": confidence or None,
                "necessaryOrSufficient": n_or_s or None,
                "identity": float(identity) if pd.notna(identity) else None,
                "note": " · ".join(note_parts) or None,
            })
            if pmid:
                used_pmids.add(pmid)

        # ---- Add Pfam DBDs (from TF_dbd_data) — only if Translation matches canonical
        dbd_rows = dbd_by_trans.get(ensp, pd.DataFrame())
        for _, r in dbd_rows.iterrows():
            try:
                s, e = int(r["From"]), int(r["To"])
            except (TypeError, ValueError):
                continue
            if s < 1 or e > length or e < s:
                continue
            nm = str(r.get("Pfam", "DBD")).strip()
            # Skip unresolved Pfam scans: "UNKNOWN" is not a real family and the
            # source lists it spanning ~the whole protein (bogus full-length DBD,
            # e.g. TET2). Real DBDs carry a proper Pfam family name (zf-C2H2 …).
            if nm.upper().replace("PFAM:", "").strip() in ("UNKNOWN", "", "NAN", "NONE"):
                continue
            reports.append({
                "name": nm,
                "type": "DNA",
                "start": s,
                "end": e,
                "pmid": "",
                "source": "CIS-BP 2.0",
                "method": "CIS-BP 2.0",
                "assay": None,
                "activity": None,
                "confidence": None,
                "necessaryOrSufficient": None,
                "identity": None,
                "note": "CIS-BP DBD",
            })
            SOURCE_SET.add("CIS-BP 2.0")

        # ---- Canonical "domains" (collapsed) for the 3D type-color view.
        # Group reports by (name, type), then merge ONLY overlapping intervals
        # within each group. Non-overlapping reports with the same name (e.g.
        # ADNP has four separate literature RD calls at positions 2-81,
        # 392-521, 592-671, 612-691) stay as four distinct domain entries.
        # Previously they collapsed to a single max-span entry, which hid most
        # domains in the heatmap stripe.
        by_key: dict[tuple, list[dict]] = {}
        for rep in reports:
            key = (rep["name"], rep["type"])
            by_key.setdefault(key, []).append(rep)
        domain_entries: list[dict] = []
        for (name, dtype), reps in by_key.items():
            reps_sorted = sorted(reps, key=lambda r: (r["start"], r["end"]))
            merged: list[tuple[int, int]] = []
            for r in reps_sorted:
                if merged and r["start"] <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], r["end"]))
                else:
                    merged.append((r["start"], r["end"]))
            # Keep Bifunctional as a first-class domain type (was previously
            # collapsed to AD for legacy rendering). All colour maps now
            # support it explicitly.
            display_type = dtype
            for s, e in merged:
                domain_entries.append({
                    "name": name,
                    "type": display_type,
                    "start": s,
                    "end": e,
                })
        domains = sorted(domain_entries, key=lambda d: d["start"])

        # ---- PTMs (canonical translation only, clipped to length)
        ptm_rows = ptms_by_trans.get(ensp)
        ptms = []
        if ptm_rows is not None and not ptm_rows.empty:
            seen: set[tuple] = set()
            for _, r in ptm_rows.iterrows():
                try:
                    pos = int(r["Real Position"])
                except (TypeError, ValueError):
                    continue
                if pos < 1 or pos > length:
                    continue
                ptype_long = str(r.get("PTM", ""))
                ptype = PTM_TYPE_MAP.get(ptype_long, "Other")
                if ptype == "Other":
                    continue  # drop sumo/glyco/etc. for now (extend later)
                key = (pos, ptype)
                if key in seen:
                    continue
                seen.add(key)
                aa = sequence[pos - 1] if 1 <= pos <= len(sequence) else "X"
                pmid = str(r.get("PMID", "")).split(".")[0].strip()
                ptms.append({"pos": pos, "type": ptype, "residue": f"{aa}{pos}",
                             **({"pmid": pmid} if pmid and pmid != "nan" else {})})
            ptms.sort(key=lambda p: (p["pos"], p["type"]))

        # ---- Isoforms (canonical first, then alternates from TF_completeIDs)
        # The source xlsx occasionally has two rows that share the same
        # Isoform name string (e.g. "ARID3B-201" at length 561 and 560). The
        # downstream UI keys "this is the canonical row" on the *name* string,
        # so a duplicate would light up both rows as canonical. Disambiguate
        # by appending the length to the second occurrence and forbid more
        # than one canonical flag per TF.
        # Emit every translated ENSP for this gene. Dedupe by Translation ID
        # (the ENSP), not by length — TP53 alone has six distinct ENSPs that
        # share the canonical 393 aa length, and the length-keyed dedupe was
        # collapsing all of them into a single row. The canonical ENSP is also
        # excluded so it doesn't appear twice in the stack.
        gene_rows = ids_df[ids_df["Gene_Upper"] == gene]
        canonical_ensp = (
            str(canon.get("Translation_Clean") or canon.get("Translation ID") or "")
            .split(".")[0]
            .strip()
        )
        # Per-isoform Pfam DBDs come from TF_dbd_data_ which is keyed by
        # Translation ID and gives From/To in each isoform's own residue frame.
        def _pfam_dbds_for(ensp: str, iso_len: int) -> list[dict]:
            rows = dbd_by_trans.get(ensp, pd.DataFrame())
            out: list[dict] = []
            for _, r in rows.iterrows():
                try:
                    s_raw = int(r["From"])
                    e_raw = int(r["To"])
                except (TypeError, ValueError):
                    continue
                # The source sometimes reports From=0 for isoforms whose DBD
                # starts at the N-terminus after truncation — clamp instead of
                # rejecting (the canonical convention treats From as 0-indexed
                # inclusive, To as 1-indexed inclusive, which matches the
                # existing canonical chip).
                # Drop DBDs that don't fit this isoform's length (>2 aa overshoot
                # = wrong-length sequence version); clamp only tiny off-by-ones.
                if s_raw > iso_len or e_raw > iso_len + 2:
                    continue
                s = max(1, s_raw if s_raw > 0 else 1)
                e = min(iso_len, e_raw)
                if e < s:
                    continue
                nm = str(r.get("Pfam", "DBD")).strip()
                # Skip unresolved "UNKNOWN" Pfam scans (bogus full-length DBD).
                if nm.upper().replace("PFAM:", "").strip() in ("UNKNOWN", "", "NAN", "NONE"):
                    continue
                out.append({
                    "name": nm,
                    "start": s,
                    "end": e,
                })
            return out

        def _sorted_effectors_for(ensp: str, iso_len: int) -> list[dict]:
            rows = effector_by_ensp.get(ensp, [])
            kept: list[dict] = []
            for r in rows:
                # Drop domains that don't fit THIS isoform's length (a >2 aa
                # overshoot means the projection used a different-length sequence
                # version, e.g. IKZF2 ENSP00000395203 — 452 aa in the source vs
                # 409 aa here). Clamping would silently fabricate a truncated
                # domain; dropping is correct. ≤2 aa overshoot is tolerated as
                # off-by-one and clamped.
                if r["start"] > iso_len or r["end"] > iso_len + 2:
                    continue
                s_, e_ = max(1, r["start"]), min(iso_len, r["end"])
                if e_ < s_:
                    continue
                kept.append({**r, "start": s_, "end": e_})
            kept.sort(key=lambda x: (
                {"original": 0, "derivada": 1, "alieneamiento": 2}.get(x["origin"], 3),
                x["start"],
                x["end"],
            ))
            return kept

        aa_isoforms = []
        seen_ensps: set[str] = {canonical_ensp} if canonical_ensp else set()
        seen_names: set[str] = set()
        # Transcript metadata accumulator (one entry per ENST). Cross-links
        # iso.name ↔ ENST so per-transcript variant data (e.g. COSMIC
        # HGVS.p per ENST) can attach to the right iso in the UI.
        aa_transcripts: list[dict] = []
        canonical_name = str(canon.get("Isoform name") or symbol)
        canonical_enst = str(canon.get("Transcript ID") or "").split(".")[0].strip()
        aa_isoforms.append({
            "name": canonical_name,
            "length": length,
            "canonical": True,
            "exons": [{"start": 1, "end": length}],
            "appris": str(canon.get("APPRIS Annotation") or "") or None,
            "proteinId": canonical_ensp or None,
            "pfamDBDs": _pfam_dbds_for(canonical_ensp, length) if canonical_ensp else [],
            "effectorDomains": _sorted_effectors_for(canonical_ensp, length) if canonical_ensp else [],
            "sequence": sequence,
            "uniprotIso": uniprot,
        })
        if canonical_enst:
            # length(nt) = 3 × aa + 3 stop codon. Real exon model lands in
            # genomicStructures from prot2exon — this Transcript[] entry
            # is for ENST cross-linking only.
            aa_transcripts.append({
                "name": canonical_name,
                "ensemblId": canonical_enst,
                "length": length * 3 + 3,
                "exons": [{"start": 1, "end": length * 3 + 3}],
                "proteinIsoform": canonical_name,
                "canonical": True,
            })
        seen_names.add(canonical_name)
        # Sort alternates: longest first so the stack mirrors the protein-size
        # gradient the canonical anchors at the top.
        gene_rows_sorted = gene_rows.dropna(subset=["Protein length"]).sort_values(
            "Protein length", ascending=False
        )
        for _, r in gene_rows_sorted.iterrows():
            iso_ensp = str(r.get("Translation_Clean") or "").split(".")[0].strip()
            if not iso_ensp or iso_ensp in seen_ensps:
                continue
            iso_len = int(r.get("Protein length"))
            # Skip 0/1-aa stubs that prot2exon can't usefully draw, but allow
            # genuine short fragments (e.g. TP53's 35-aa ENSP).
            if iso_len < 20:
                continue
            raw_name = str(r.get("Isoform name") or f"isoform-{iso_ensp}")
            # Collide-safe name: append (Naa) for duplicates so the UI never
            # ends up with two isoforms sharing a single key.
            iso_name = raw_name
            if iso_name in seen_names:
                iso_name = f"{raw_name} ({iso_len} aa)"
            # Exon range is a placeholder unless prot2exon overrides it later.
            # Anchor short isoforms at the C-terminus (offset = length - iso_len)
            # for the canonical case and full-range for longer ones.
            if iso_len < length:
                offset = length - iso_len
                exons = [{"start": offset + 1, "end": length}]
            else:
                exons = [{"start": 1, "end": iso_len}]
            # Iso sequence from TF_completeIDs — drives the sequence panel
            # that pops up under the AA tab when a non-canonical iso is
            # selected (FASTA preview + diff vs canonical + download).
            iso_seq_raw = str(r.get("Sequence aa") or "").strip()
            iso_seq = iso_seq_raw if iso_seq_raw and iso_seq_raw.lower() != "nan" else None
            iso_uniprot_full = str(r.get("UniProt_ID") or "").strip()
            iso_uniprot = (
                iso_uniprot_full if iso_uniprot_full and iso_uniprot_full.lower() != "nan" else None
            )
            entry = {
                "name": iso_name,
                "length": iso_len,
                "exons": exons,
                "appris": str(r.get("APPRIS Annotation") or "") or None,
                "proteinId": iso_ensp,
                "pfamDBDs": _pfam_dbds_for(iso_ensp, iso_len),
                "effectorDomains": _sorted_effectors_for(iso_ensp, iso_len),
            }
            if iso_seq:
                entry["sequence"] = iso_seq
            if iso_uniprot:
                entry["uniprotIso"] = iso_uniprot
            aa_isoforms.append(entry)
            seen_ensps.add(iso_ensp)
            seen_names.add(iso_name)
            iso_enst = str(r.get("Transcript ID") or "").split(".")[0].strip()
            if iso_enst:
                aa_transcripts.append({
                    "name": iso_name,
                    "ensemblId": iso_enst,
                    "length": iso_len * 3 + 3,
                    "exons": [{"start": 1, "end": iso_len * 3 + 3}],
                    "proteinIsoform": iso_name,
                    "canonical": False,
                })

        # ---- Copy PDB
        pdb_src = os.path.join(SRC_PDB, f"AF-{uniprot}-F1.pdb")
        pdb_dest = os.path.join(DEST_PDB, f"AF-{uniprot}-F1.pdb")
        if os.path.exists(pdb_src):
            if not os.path.exists(pdb_dest):
                shutil.copy(pdb_src, pdb_dest)
            pdb_bytes_total += os.path.getsize(pdb_dest)

        # ---- PAE
        if not args.skip_pae:
            pae_bytes = write_pae(uniprot, length)
            if pae_bytes:
                pae_bytes_total += pae_bytes

        # ---- Conservation
        cons_bytes = write_conservation(uniprot)
        if cons_bytes:
            cons_bytes_total += cons_bytes

        # ---- Papers list — enrich every PMID with whatever title we have:
        # Hoja 1 covers most curated/HumanizedTable rows; PAPER_TITLE_FALLBACK
        # carries the HT-MPRA flagship paper titles.
        papers = []
        for p in sorted(used_pmids):
            if not p or p == "nan":
                continue
            entry: dict[str, str] = {"pmid": p}
            if p in pmid_titles:
                entry["title"] = pmid_titles[p]
            elif p in PAPER_TITLE_FALLBACK:
                entry["title"] = PAPER_TITLE_FALLBACK[p]
            papers.append(entry)

        out = {
            "symbol": symbol,
            "uniprot": uniprot,
            # Reviewed Swiss-Prot accession for the public UniProt link / page
            # identity (Route 2). `uniprot` above stays the asset + coordinate
            # key; this only changes what the cross-reference strip shows/links.
            "reviewedAccession": _reviewed_accession(gene, uniprot),
            "ensemblGene": ensg if ensg and ensg != "nan" else None,
            "ensemblProtein": ensp if ensp and ensp != "nan" else None,
            # Prefer the Lambert 2018 EntrezGene Description ("tumor protein p53")
            # over the bare HGNC symbol so the heatmap rows can show a real
            # human-readable protein name.
            "name": lambert_name.get(gene) or str(canon.get("Gene name (HGNC)") or symbol),
            # Family preference: Lambert 2018 > v1 'TF Family' col > Pfam fallback.
            "family": normalize_family(
                lambert_family.get(gene)
                or (
                    str(tf_rows["TF Family"].dropna().iloc[0])
                    if not tf_rows.empty and tf_rows["TF Family"].notna().any()
                    else _family_from_pfam(dbd_by_trans.get(ensp, pd.DataFrame()))
                )
            ),
            "length": length,
            "sequence": sequence,
            "domains": domains,
            "ptms": ptms,
            "mutations": [],
            "isoforms": aa_isoforms,
            "transcripts": aa_transcripts,
            "papers": papers,
            "domainReports": reports,
            "hasConservation": cons_bytes is not None,
            "hasPDB": os.path.exists(pdb_dest),
            "hasPAE": not args.skip_pae and os.path.exists(os.path.join(DEST_PAE, f"{uniprot}.json")),
        }

        with open(os.path.join(DEST_MOCK, f"{symbol}.json"), "w") as f:
            json.dump(out, f, separators=(",", ":"))
        written.append(symbol)

        if idx % 100 == 0 or idx == len(eligible_genes):
            print(f"  [{idx}/{len(eligible_genes)}] {symbol:<10} {uniprot:<8} {length:>4}aa "
                  f"reports={len(reports):<3} ptms={len(ptms):<4} cons={cons_bytes is not None}")

    print()
    print(f"Written: {len(written)} TFs")
    print(f"  skipped (no canonical record): {len(skipped_no_canon)}")
    print(f"  skipped (no PDB and no IDs seq): {len(skipped_no_pdb)}")
    print(f"  PDB bytes copied:  {pdb_bytes_total/1024/1024:.0f} MB")
    print(f"  PAE bytes written: {pae_bytes_total/1024/1024:.0f} MB")
    print(f"  Cons bytes written:{cons_bytes_total/1024/1024:.1f} MB")
    print(f"  Method sources seen: {sorted(SOURCE_SET)}")

    # Emit a TSV of out-of-bounds domain rows for curator review. Each row is
    # tagged "clipped" (≤20 aa overflow, kept on the page) or "dropped"
    # (>20 aa overflow, excluded from the page).
    if OUT_OF_BOUNDS_LOG:
        audit_path = os.path.join(ROOT, "scripts", "out_of_bounds_audit.tsv")
        with open(audit_path, "w") as f:
            f.write("symbol\tlength\tuniprot\tcoords\tend_over_by\taction\tsource\tpmid\ttype\trow_id\n")
            for r in OUT_OF_BOUNDS_LOG:
                f.write("\t".join(str(r[k]) for k in [
                    "symbol", "length", "uniprot", "coords",
                    "end_over_by", "action", "source", "pmid", "type", "row_id",
                ]) + "\n")
        n_drop = sum(1 for r in OUT_OF_BOUNDS_LOG if r["action"] == "dropped")
        n_clip = len(OUT_OF_BOUNDS_LOG) - n_drop
        print(
            f"\n[validate] {len(OUT_OF_BOUNDS_LOG)} out-of-bounds domain rows "
            f"({n_drop} dropped, {n_clip} clipped). Audit written to:\n  {audit_path}\n"
            f"  Dropped rows usually mean curator referenced a longer isoform"
            f" than canonical (e.g. XBP1s 376aa vs XBP1u 261aa)."
        )

    print(f"\n[6/6] Regenerating public/mock/index.json...")
    regen_index(lambert_universe)

    # Optional final step: fold prot2exon genomic mapping into each per-TF JSON.
    # Skipped silently if the prot2exon output directory hasn't been produced
    # yet (run `python3 scripts/wire_genomic.py` separately, or rerun this
    # script after producing prot2exon results).
    if args.skip_isoform_wire:
        print(f"\n[bonus] Skipping prot2exon isoform wiring (--skip-isoform-wire).")
    else:
        try:
            # scripts/ isn't on sys.path when running as `python3 scripts/foo.py`,
            # so make sibling imports work.
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            import wire_isoform_genomics  # type: ignore
            if os.path.exists(wire_isoform_genomics.PROT2EXON_BIN) and \
               os.path.exists(wire_isoform_genomics.INDEX):
                print(f"\n[bonus] Mapping per-isoform genomics via prot2exon...")
                wire_isoform_genomics.main()
            else:
                print(
                    f"\n[bonus] Skipping prot2exon isoform wiring (binary or index "
                    f"missing): {wire_isoform_genomics.PROT2EXON_BIN} / "
                    f"{wire_isoform_genomics.INDEX}"
                )
        except Exception as e:
            print(f"  (warning: prot2exon isoform wiring step failed: {e})")


# ---------------------------------------------------------------------- index
# Fallback denominators used only when Lambert 2018 doesn't bucket a family
# under the same normalised name. Lambert-derived counts always take priority
# (see lambert_universe in main()).
FAMILY_SIZE_HINT = {
    "C2H2 ZF": 747, "C2H2 zinc finger": 747,
    "Homeo": 196, "Homeodomain": 196,
    "bHLH": 108, "Basic helix-loop-helix": 108,
    "Forkhead": 49,
    "Nuclear receptor": 48,
    "bZIP": 54, "Basic leucine zipper": 54,
    "ETS": 27, "Ets": 27,
    "SOX": 30, "HMG box": 30, "HMG_box": 30,
    "E2F": 11, "E2F_TDP": 11,
    "P53": 3,
    "AT-hook": 16, "AT hook": 16,
}


def regen_index(lambert_universe: dict[str, int] | None = None):
    by_family = collections.defaultdict(list)
    method_tfs = collections.defaultdict(set)
    domain_type_tfs = collections.defaultdict(set)
    # Collect unique PMIDs across every per-TF JSON so the index's
    # `papers` total counts distinct references, not the sum-across-TFs
    # double-counting fix the audit flagged.
    unique_pmids: set[str] = set()

    for fn in sorted(os.listdir(DEST_MOCK)):
        if not fn.endswith(".json"):
            continue
        with open(os.path.join(DEST_MOCK, fn)) as f:
            tf = json.load(f)
        sym = tf["symbol"]
        for p in tf.get("papers", []) or []:
            pid = str(p.get("pmid") or "").strip()
            if pid and pid.lower() not in ("nan", "none"):
                unique_pmids.add(pid)
        fam = normalize_family(tf.get("family", "Unclassified") or "Unclassified")
        n_ad = sum(1 for r in tf.get("domainReports", []) if r["type"] == "AD")
        n_rd = sum(1 for r in tf.get("domainReports", []) if r["type"] == "RD")
        n_dna = sum(1 for r in tf.get("domainReports", []) if r["type"] == "DNA")
        n_bif = sum(1 for r in tf.get("domainReports", []) if r["type"] == "Bifunctional")
        n_ptms = len(tf.get("ptms", []))
        n_papers = len(tf.get("papers", []))
        # "bifunctional" counts genuine Bifunctional-typed domain calls (like
        # ad/rd count AD/RD calls) — NOT the old "has both an AD and an RD"
        # heuristic, which conflated functionally-bifunctional TFs with TFs that
        # actually have a bifunctional domain and left the Browse "Bifunctional"
        # tab listing TFs whose rows never rendered (no Bifunctional-typed row).
        bifunc = n_bif
        # Per-TF Pfam DBD names, deduped + ordered by first appearance.
        # Surfaces e.g. "P53" / "Homeobox" / "zf-C2H2" on the browse table so
        # users can scan the DBD family at a glance.
        pfam_names: list[str] = []
        seen_pf: set[str] = set()
        for rep in tf.get("domainReports", []):
            if rep.get("source") != "CIS-BP 2.0":
                continue
            nm = str(rep.get("name") or "").strip()
            if nm and nm not in seen_pf:
                seen_pf.add(nm)
                pfam_names.append(nm)
        by_family[fam].append({
            "symbol": sym,
            "name": tf.get("name") or sym,
            "uniprot": tf["uniprot"],
            # Ensembl IDs included so the homepage SearchBar can match by
            # ENSG/ENSP — pulled from the per-TF JSON (already populated
            # from TF_completeIDs). Falls back to the legacy ensemblId
            # field for older builds.
            "ensemblGene": tf.get("ensemblGene") or tf.get("ensemblId") or None,
            "ensemblProtein": tf.get("ensemblProtein") or None,
            "length": tf["length"],
            "counts": {
                "ad": n_ad, "rd": n_rd, "dna": n_dna,
                "ptms": n_ptms, "papers": n_papers, "bifunctional": bifunc,
            },
            "domains": [
                {"start": d["start"], "end": d["end"], "type": d["type"]}
                for d in tf.get("domains", [])
            ],
            "pfamDBDs": pfam_names,
        })
        if n_ad: domain_type_tfs["Activation domain"].add(sym)
        if n_rd: domain_type_tfs["Repression domain"].add(sym)
        if n_bif: domain_type_tfs["Bifunctional"].add(sym)
        if n_dna: domain_type_tfs["DNA-binding domain"].add(sym)
        for rep in tf.get("domainReports", []):
            m = rep.get("method") or rep.get("source")
            if m and m != "CIS-BP 2.0":
                method_tfs[m].add(sym)

    families_out = []
    lu = lambert_universe or {}
    for fam, tfs in sorted(by_family.items(), key=lambda kv: -len(kv[1])):
        tfs_sorted = sorted(tfs, key=lambda t: t["symbol"])
        agg = {k: sum(t["counts"][k] for t in tfs_sorted) for k in ["ad", "rd", "dna", "ptms", "papers", "bifunctional"]}
        # Pick the universe denominator with this priority:
        #   1. Lambert 2018 count for the same normalised family name.
        #   2. FAMILY_SIZE_HINT (curated fallback for families Lambert lumps
        #      differently, e.g. multi-DBD combinations).
        #   3. 2× the curated count (last-ditch heuristic).
        # Whatever we pick must be >= the curated count, so we never report a
        # nonsense fraction like 9 / 8.
        curated_n = len(tfs_sorted)
        universe = lu.get(fam) or FAMILY_SIZE_HINT.get(fam) or max(10, curated_n * 2)
        if universe < curated_n:
            universe = curated_n
        families_out.append({
            "name": fam,
            "curated": curated_n,
            "totalInFamily": universe,
            "counts": {**agg, "tfs": curated_n},
            "tfs": tfs_sorted,
        })

    domain_types_out = [
        {"name": name, "tfs": sorted(tfs), "count": len(tfs), "subtitle": ""}
        for name, tfs in sorted(domain_type_tfs.items())
    ]
    methods_out = [
        {"name": name, "tfs": sorted(tfs), "count": len(tfs),
         "subtitle": f"{len(tfs)} TFs"}
        for name, tfs in sorted(method_tfs.items(), key=lambda kv: -len(kv[1]))
    ]
    totals = {
        "tfs": sum(f["counts"]["tfs"] for f in families_out),
        "families": len(families_out),
        # Unique PMIDs across the whole catalog (not sum-across-TFs, which
        # double-counted papers that cite multiple TFs).
        "papers": len(unique_pmids),
        # Sum-across-TFs is still useful for "total paper-TF curation
        # events" — keep it under a less misleading name.
        "paperRefs": sum(f["counts"]["papers"] for f in families_out),
        "ad": sum(f["counts"]["ad"] for f in families_out),
        "rd": sum(f["counts"]["rd"] for f in families_out),
        "ptms": sum(f["counts"]["ptms"] for f in families_out),
    }
    with open(INDEX_OUT, "w") as f:
        json.dump({
            "totals": totals,
            "families": families_out,
            "domainTypes": domain_types_out,
            "methods": methods_out,
        }, f, separators=(",", ":"))
    print(f"  index.json: {totals}")


# ---------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="build only first N TFs")
    ap.add_argument("--symbols", type=str, default=None, help="comma-separated HGNC symbols")
    ap.add_argument("--skip-pae", action="store_true", help="skip PAE write (saves ~530 MB + ~30 s)")
    ap.add_argument("--skip-isoform-wire", action="store_true",
                    help="skip the prot2exon isoform-genomic mapping step (saves ~60 s)")
    ap.add_argument("--fast", action="store_true",
                    help="dev shortcut: implies --skip-pae --skip-isoform-wire")
    args = ap.parse_args()
    if args.fast:
        args.skip_pae = True
        args.skip_isoform_wire = True

    # Sanity: source folder must exist
    for p in [SRC_UNIFIED, SRC_PTMS, SRC_IDS, SRC_DBD, SRC_PDB, SRC_PAE, SRC_CONS]:
        if not os.path.exists(p):
            print(f"ERROR: missing source: {p}", file=sys.stderr)
            sys.exit(1)

    build(args)


if __name__ == "__main__":
    main()
