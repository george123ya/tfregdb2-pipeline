"""Convert returned per-isoform ConSurf grades into per-isoform JSON.

Run locally after the SCC array job (consurf_array_isoforms.sh) returns its
`conservation_scores/{base}_scores.txt` files (the flattened
no_model_consurf_grades.txt). Using consurf_isoforms/manifest.json:

  - parse each scores file (the COLOR column is the 1-9 grade; '*' = low
    confidence, stripped) into the site's {pos, grade, aa} shape
  - write public/conservation/{isoformName}.json for EVERY isoform sharing
    that sequence (dedup expansion)
  - for canonical-equal isoforms, copy the existing canonical
    public/conservation/{uniprot}.json to {isoformName}.json

Mirrors tfregdb/data_pipeline/parse_consurf.py, but keys by isoform name and
emits `grade` (not `color`) to match ConservationViewer.tsx.

Usage:
    protein2genomic/.venv/bin/python scripts/parse_consurf_isoforms.py \
        [scores_dir]   # default: consurf_isoforms/conservation_scores
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MANIFEST = REPO / "consurf_isoforms" / "manifest.json"
OUT_DIR = REPO / "public" / "conservation"


def parse_grades(path: Path) -> list[dict]:
    residues: list[dict] = []
    data = False
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("POS") and "SCORE" in s:
            data = True
            continue
        if not data:
            continue
        cols = s.split()
        if not cols[0].isdigit():
            continue
        try:
            residues.append({
                "pos": int(cols[0]),
                "grade": int(cols[3].replace("*", "")),  # COLOR column, 1-9
                "aa": cols[1],
            })
        except (ValueError, IndexError):
            continue
    return residues


def main() -> None:
    scores_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        REPO / "consurf_isoforms" / "conservation_scores"
    manifest = json.loads(MANIFEST.read_text())
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    written = missing_scores = 0
    for base, names in manifest["runs"].items():
        sf = scores_dir / f"{base}_scores.txt"
        if not sf.exists():
            missing_scores += 1
            continue
        residues = parse_grades(sf)
        if not residues:
            continue
        payload = json.dumps(residues, separators=(",", ":"))
        for name in names:
            (OUT_DIR / f"{name}.json").write_text(payload)
            written += 1

    copied = missing_canon = 0
    for name, uni in manifest["canonical_equal"].items():
        src = OUT_DIR / f"{uni}.json"
        if not src.exists():
            missing_canon += 1
            continue
        shutil.copyfile(src, OUT_DIR / f"{name}.json")
        copied += 1

    print(f"[consurf-out] isoform files written from runs: {written}")
    print(f"[consurf-out] canonical-equal copied: {copied}")
    if missing_scores:
        print(f"[consurf-out] WARNING: {missing_scores} runs missing a scores "
              f"file (job incomplete/failed for those)")
    if missing_canon:
        print(f"[consurf-out] WARNING: {missing_canon} canonical-equal isoforms "
              f"had no canonical JSON to copy")


if __name__ == "__main__":
    main()
