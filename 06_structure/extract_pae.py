#!/usr/bin/env python
"""Extract the PAE matrix of the top-ranked model into the site's JSON format.

AF2 keeps predicted_aligned_error only inside result_model_*.pkl, which for a
4 kaa chain is ~8 GB (the L x L x 64 distogram / aligned_confidence_probs
tensors are the bulk). We load the one pkl that ranked_0 came from, keep the
L x L PAE, and drop the rest.

Output matches what the site already serves in public/pae/:
    {"length": N, "matrix": [[int, ...], ...]}
PAE is rounded to the nearest integer angstrom, same as the AFDB downloads.

Also stages ranked_0.pdb as AF-{ID}-F1.pdb so the layout mirrors public/pdb/.
"""

import argparse
import json
import os
import pickle
import shutil
import sys

import numpy as np

PAE_KEYS = ("predicted_aligned_error", "predicted_aligned_error_ptm")


def top_model(out_dir, tid):
    """Name of the model that became ranked_0, per ranking_debug.json."""
    path = os.path.join(out_dir, tid, "ranking_debug.json")
    with open(path) as fh:
        dbg = json.load(fh)
    order = dbg.get("order")
    if not order:
        raise ValueError("%s has no 'order' key" % path)
    return order[0]


def load_pae(pkl_path):
    with open(pkl_path, "rb") as fh:
        result = pickle.load(fh)
    for key in PAE_KEYS:
        if key in result:
            return np.asarray(result[key])
    raise KeyError(
        "no PAE in %s (keys: %s)" % (pkl_path, ", ".join(sorted(result)))
    )


def write_pae_json(pae, dest):
    """Stream the matrix out row by row; json.dumps on 16M ints is wasteful."""
    mat = np.rint(pae).astype(np.int16)
    n = mat.shape[0]
    tmp = dest + ".part"
    with open(tmp, "w") as fh:
        fh.write('{"length":%d,"matrix":[' % n)
        for i in range(n):
            if i:
                fh.write(",")
            fh.write("[" + ",".join(map(str, mat[i].tolist())) + "]")
        fh.write("]}")
    os.replace(tmp, dest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="af2_output directory")
    ap.add_argument("--dest", required=True, help="staging dir (gets pae/ and pdb/)")
    ap.add_argument("--ids", nargs="+", required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    pae_dir = os.path.join(args.dest, "pae")
    pdb_dir = os.path.join(args.dest, "pdb")
    os.makedirs(pae_dir, exist_ok=True)
    os.makedirs(pdb_dir, exist_ok=True)

    failures = 0
    for tid in args.ids:
        src_pdb = os.path.join(args.out, tid, "ranked_0.pdb")
        dst_pdb = os.path.join(pdb_dir, "AF-%s-F1.pdb" % tid)
        dst_pae = os.path.join(pae_dir, "%s.json" % tid)

        if not os.path.exists(src_pdb):
            print("SKIP %s: no ranked_0.pdb" % tid, flush=True)
            failures += 1
            continue

        if os.path.exists(dst_pae) and os.path.exists(dst_pdb) and not args.force:
            print("HAVE %s" % tid, flush=True)
            continue

        try:
            model = top_model(args.out, tid)
            pkl = os.path.join(args.out, tid, "result_%s.pkl" % model)
            pae = load_pae(pkl)
            write_pae_json(pae, dst_pae)
            shutil.copyfile(src_pdb, dst_pdb)
            print(
                "OK %s model=%s L=%d pae_max=%.1f"
                % (tid, model, pae.shape[0], float(pae.max())),
                flush=True,
            )
        except Exception as exc:  # noqa: BLE001 - want the id alongside the error
            print("FAIL %s: %s: %s" % (tid, type(exc).__name__, exc), flush=True)
            failures += 1

    if failures:
        print("failures: %d" % failures, flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
