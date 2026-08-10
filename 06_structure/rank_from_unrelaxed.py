#!/usr/bin/env python3
"""Produce ranked_*.pdb + ranking_debug.json from already-predicted unrelaxed models.

Why this exists: the >2700 aa giants finish inference fine, then die in the AMBER
relax step with

    MemoryError: Unable to allocate 34.2 GiB for an array with shape
    (3953, 3953, 14, 14, 3)   in all_atom.between_residue_clash_loss

That allocation is the violation check inside amber_minimize, and it scales as
O(N^2 * 14 * 14 * 3), so it is hopeless for a 4 kaa chain. The models themselves
are already on disk. AlphaFold under --models_to_relax=none ranks the UNRELAXED
models and writes ranked_*.pdb from them; this script does exactly that, so the
structures become usable without burning another ~30 GPU-hours.

Ranking metric: for the monomer / monomer_ptm presets AlphaFold's
`ranking_confidence` is mean pLDDT, and it writes per-residue pLDDT into the
B-factor column of the unrelaxed PDB. So the mean CA B-factor reproduces the
ranking exactly, without opening the multi-GB result_model_*.pkl files.

NOTE: ranked_0.pdb written here is UNRELAXED, unlike the other 414 targets which
were relaxed with -r best. A dropped `RANKED_FROM_UNRELAXED.txt` marks each one.
"""

import argparse
import json
import os
import re
import shutil
import sys

MODEL_RE = re.compile(r"^unrelaxed_(model_\d+.*)\.pdb$")


def mean_ca_bfactor(pdb_path):
    """Mean B-factor over CA atoms == mean pLDDT == AF2's ranking_confidence."""
    total, n = 0.0, 0
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            # columns are fixed-width per the PDB spec: atom name 13-16, B 61-66
            if line[12:16].strip() != "CA":
                continue
            total += float(line[60:66])
            n += 1
    if n == 0:
        raise ValueError(f"no CA atoms found in {pdb_path}")
    return total / n


def rank_one(target_dir, force=False):
    name = os.path.basename(target_dir.rstrip("/"))
    ranked0 = os.path.join(target_dir, "ranked_0.pdb")
    if os.path.exists(ranked0) and not force:
        return f"{name}: ranked_0.pdb already exists, skipped"

    models = {}
    for fn in sorted(os.listdir(target_dir)):
        m = MODEL_RE.match(fn)
        if m:
            models[m.group(1)] = os.path.join(target_dir, fn)
    if not models:
        return f"{name}: NO unrelaxed models, nothing to rank"

    scores = {k: mean_ca_bfactor(v) for k, v in models.items()}
    order = sorted(scores, key=scores.get, reverse=True)

    for i, model in enumerate(order):
        shutil.copyfile(models[model], os.path.join(target_dir, f"ranked_{i}.pdb"))

    with open(os.path.join(target_dir, "ranking_debug.json"), "w") as fh:
        json.dump({"plddts": scores, "order": order}, fh, indent=4)

    with open(os.path.join(target_dir, "RANKED_FROM_UNRELAXED.txt"), "w") as fh:
        fh.write(
            "ranked_*.pdb here were ranked from the UNRELAXED models by mean CA\n"
            "B-factor (= mean pLDDT = AlphaFold's ranking_confidence), because the\n"
            "AMBER relax step OOMs on chains this long. Equivalent to running\n"
            "AlphaFold with --models_to_relax=none. No relaxation was applied.\n"
            f"models ranked: {len(order)}\n"
        )

    top = order[0]
    return (
        f"{name}: {len(order)} models ranked, best={top} "
        f"(mean pLDDT {scores[top]:.2f}) -> ranked_0.pdb"
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="af2_output",
                    help="AlphaFold output root holding one dir per target")
    ap.add_argument("--ids", nargs="+", required=True,
                    help="target IDs to rank (dir names under --out)")
    ap.add_argument("--force", action="store_true",
                    help="rewrite ranked_*.pdb even if ranked_0.pdb exists")
    args = ap.parse_args()

    if not os.path.isdir(args.out):
        sys.exit(f"error: --out '{args.out}' is not a directory")

    rc = 0
    for tid in args.ids:
        d = os.path.join(args.out, tid)
        if not os.path.isdir(d):
            print(f"{tid}: NO output directory under {args.out}")
            rc = 1
            continue
        try:
            print(rank_one(d, force=args.force))
        except Exception as exc:                      # keep going over the rest
            print(f"{tid}: FAILED - {exc}")
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
