"""Run DSSP against each bundled AlphaFold PDB and emit a per-residue
secondary-structure string the TF page can ingest at build time.

Replaces `mockSecondaryStructure` for the 1,265 TFs that have an
AlphaFold PDB on disk. The 44 TFs without a PDB fall back to the mock
(StructureViewer's "SecStr" colour mode handles the empty-array case).

DSSP outputs 8 states (H, G, I, E, B, T, S, ' '). We collapse to the 3
buckets the UI uses:
    H / G / I → H  (helix family — α, 3₁₀, π)
    E / B     → E  (strand family — extended, β-bridge)
    T / S / ' ' → C  (turn / bend / coil)

Writes `public/mock/secondary/{uniprot}.json` of shape `{"ss": "HHCCC..."}`.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP


PUBLIC = Path("public")
PDB_DIR = PUBLIC / "pdb"
OUT_DIR = PUBLIC / "mock" / "secondary"
# mkdssp from the conda env (bioconda's dssp 4.x); Bio.PDB.DSSP accepts
# either an absolute path or a bare command name.
MKDSSP_BIN = "/home/goguxor/miniconda3/envs/mpratree_app_env/bin/mkdssp"

# DSSP 8-state → 3-state collapse.
SS_MAP = {
    "H": "H", "G": "H", "I": "H",
    "E": "E", "B": "E",
    "T": "C", "S": "C", " ": "C", "-": "C", "P": "C",
}


def extract_one(pdb_path: Path, parser: PDBParser) -> str | None:
    """Returns a per-residue SS string ('HHHCCCEEE…'), or None on failure."""
    try:
        struct = parser.get_structure(pdb_path.stem, str(pdb_path))
    except Exception:
        return None
    model = next(iter(struct), None)
    if model is None:
        return None
    try:
        dssp = DSSP(model, str(pdb_path), dssp=MKDSSP_BIN)
    except Exception:
        return None
    # DSSP iterates residues in PDB order. Each key is (chain_id, residue id).
    # AlphaFold monomers ship a single chain A with residues 1..N. Pull the
    # 3-state SS letter in order.
    out: list[str] = []
    last_resi = 0
    for key in dssp.keys():
        try:
            resi = key[1][1]
        except Exception:
            continue
        if resi <= last_resi:
            continue  # guard against alt-loc duplicates
        last_resi = resi
        info = dssp[key]
        # info shape: (idx, aa, ss, accessibility, phi, psi, NH→O_1_relidx,
        #              NH→O_1_energy, O→NH_1_relidx, O→NH_1_energy, ...)
        ss = info[2]
        out.append(SS_MAP.get(ss, "C"))
    return "".join(out) if out else None


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    parser = PDBParser(QUIET=True)
    pdbs = sorted(PDB_DIR.glob("AF-*-F1.pdb"))
    print(f"[extract_secondary] running DSSP across {len(pdbs)} AlphaFold PDBs")
    n_ok = n_skip = 0
    for i, pdb in enumerate(pdbs, 1):
        uniprot = pdb.stem.removeprefix("AF-").removesuffix("-F1")
        ss = extract_one(pdb, parser)
        if ss is None:
            n_skip += 1
            continue
        out = OUT_DIR / f"{uniprot}.json"
        out.write_text(json.dumps({"ss": ss}, separators=(",", ":")))
        n_ok += 1
        if i % 100 == 0:
            print(f"  [{i}/{len(pdbs)}] {n_ok} ok · {n_skip} skipped")
    print(f"  wrote {n_ok:,} secondary-structure files · skipped {n_skip} · → {OUT_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
