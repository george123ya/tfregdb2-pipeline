"""Extract real per-residue pLDDT from the bundled AlphaFold PDBs.

AlphaFold encodes pLDDT in the B-factor column of every ATOM record (cols
61-66 in PDB fixed-width). Per-residue pLDDT = the B-factor of the
residue's CA atom (every atom in a residue shares the same B-factor in
AlphaFold output, but CA is the canonical pick).

Writes `public/mock/plddt/{uniprot}.json` shaped `[float, …]` of length =
residue count. [symbol].astro reads this at build time when present and
falls back to mockPlddt() otherwise.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


PUBLIC = Path("public")
PDB_DIR = PUBLIC / "pdb"
OUT_DIR = PUBLIC / "mock" / "plddt"


def extract_one(pdb_path: Path) -> list[float] | None:
    """Return ordered list of per-residue pLDDT (B-factor of each CA)."""
    plddt: list[float] = []
    last_resi = None
    with pdb_path.open() as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                resi = int(line[22:26])
                b = float(line[60:66])
            except ValueError:
                continue
            # PDB is residue-sequential; guard against duplicate CAs from
            # alt-loc records (col 17) which AlphaFold doesn't emit but
            # being defensive is cheap.
            if last_resi is not None and resi == last_resi:
                continue
            plddt.append(round(b, 2))
            last_resi = resi
    return plddt if plddt else None


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdbs = sorted(PDB_DIR.glob("AF-*-F1.pdb"))
    print(f"[extract_plddt] found {len(pdbs)} AlphaFold PDBs")
    n_ok = n_skip = 0
    for pdb in pdbs:
        # Filename pattern: AF-{UNIPROT}-F1.pdb
        uniprot = pdb.stem.removeprefix("AF-").removesuffix("-F1")
        plddt = extract_one(pdb)
        if plddt is None:
            n_skip += 1
            continue
        out = OUT_DIR / f"{uniprot}.json"
        out.write_text(json.dumps(plddt, separators=(",", ":")))
        n_ok += 1
    print(f"  wrote {n_ok:,} pLDDT files · skipped {n_skip} · → {OUT_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
