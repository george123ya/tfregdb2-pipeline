"""Merge the intronic gnomAD genome variants (fetch_gnomad_introns.py output)
into the per-TF variant JSONs the genomic view reads.

Each new row becomes a record in the same shape the existing gnomAD coding
variants use (source=gnomAD, pos=0 → no protein coordinate, so it only appears
in the Genomic view, never the protein lollipop). Rows whose (genomicPos, ref,
alt) already exist in the file (i.e. coding gnomAD already pulled by the exome
pipeline) are skipped, so we never double-count — the net addition is the
intronic / deep-non-coding population variation.

Idempotent: re-running drops previously-merged intronic rows first (id prefix
"gnomad-intron:") and re-adds, so the file converges regardless of run count.

Run from repo root:
    python3 scripts/variants/merge_gnomad_introns.py
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

TSV = Path("data/variants/gnomad_introns.tsv")
TF_DIR = Path("public/mock/variants")
ID_PREFIX = "gnomad-intron:"


def main() -> int:
    if not TSV.exists():
        sys.exit(f"missing {TSV} — run fetch_gnomad_introns.py first")

    by_sym: dict[str, list[dict]] = defaultdict(list)
    with TSV.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            by_sym[row["symbol"]].append(row)
    print(f"[merge] {sum(len(v) for v in by_sym.values()):,} intronic rows "
          f"across {len(by_sym)} genes", flush=True)

    n_added = n_dup = n_files = 0
    before = after = 0
    for sym, rows in by_sym.items():
        fp = TF_DIR / f"{sym}.json"
        if not fp.exists():
            continue
        before += fp.stat().st_size
        d = json.loads(fp.read_text())
        variants = d.get("variants", [])
        # Drop any prior intronic merge, then rebuild the seen-set from the
        # remaining (coding) variants so we dedup against real coding calls.
        variants = [v for v in variants if not str(v.get("id", "")).startswith(ID_PREFIX)]
        seen = {(v.get("genomicPos"), v.get("ref"), v.get("alt")) for v in variants}
        for row in rows:
            pos = int(row["pos"])
            ref, alt = row["ref"], row["alt"]
            if (pos, ref, alt) in seen:
                n_dup += 1
                continue
            af = float(row["af"]) if row["af"] else None
            ac = int(row["ac"]) if row["ac"] else None
            an = int(row["an"]) if row["an"] else None
            pop = float(row["popmax_af"]) if row["popmax_af"] else None
            variants.append({
                "id": f"{ID_PREFIX}{row['chrom']}:{pos}{ref}>{alt}",
                "pos": 0,
                "type": "Other",
                "count": max(1, ac or 1),
                "source": "gnomAD",
                "chrom": row["chrom"],
                "genomicPos": pos,
                "ref": ref,
                "alt": alt,
                "aaRef": None,
                "aaAlt": None,
                "consequence": "intronic",
                "rsid": None,
                "clinvar": None,
                "gnomad": {"af": af, "ac": ac, "an": an, "popmaxAf": pop},
                "cosmicId": None,
                "isoforms": None,
            })
            n_added += 1
        d["variants"] = variants
        fp.write_text(json.dumps(d, separators=(",", ":")))
        after += fp.stat().st_size
        n_files += 1

    print(f"[merge] +{n_added:,} intronic gnomAD variants into {n_files} files "
          f"({n_dup:,} coding dups skipped)", flush=True)
    print(f"[merge] variants size {before/1e6:.0f} MB → {after/1e6:.0f} MB "
          f"(+{(after-before)/1e6:.0f} MB for merged genes)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
