"""Extract finches per-pair 2D interaction maps as COMPACT BINARY (SCC, himem).

Full-JSON per-pair matrices would be ~17 GB (over the R2 10 GB free tier), so we
quantise each epsilon matrix to int8 and write a tiny binary file per pair:

    scores2d/{TF}/{coords}__{partner}.bin
      bytes 0-1 : uint16 LE  rows  (partner positions)
      bytes 2-3 : uint16 LE  cols  (domain window positions)
      bytes 4.. : int8 x rows*cols, row-major
                  value q in [-127,127] = round(clamp(eps,-3,3)/3 * 127)
                  q == -128  => no data (NaN)

The site decodes this (eps = q/127*3) and computes the marginal profiles
client-side. ~1 byte/cell instead of ~6 → the whole set is ~3 GB. Matrix
orientation (from s5): rows = partner surface, cols = TF-IDR window;
eps<0 attractive (purple), eps>0 repulsive (green). Copy scores2d/ to the site
under public/interactions/scores2d/.
"""

from __future__ import annotations

import gzip
import pickle
import struct
import sys
from pathlib import Path

import numpy as np

PKL = sys.argv[1] if len(sys.argv) > 1 else "all_interactions.pkl"
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "scores2d")


def _open(p: str):
    return gzip.open(p, "rb") if p.endswith(".gz") else open(p, "rb")


def main() -> None:
    print(f"[intermaps] loading {PKL} (needs ~30-48 GB RAM)...", flush=True)
    with _open(PKL) as fh:
        data = pickle.load(fh)
    print(f"[intermaps] {len(data)} IDR keys -> binary per-pair maps in {OUT}/", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)

    n_idr = n_pair = 0
    for idr_name, partners in data.items():
        tf = idr_name.rsplit("_", 1)[0]
        coords = idr_name.rsplit("_", 1)[1]
        d = OUT / tf
        d.mkdir(exist_ok=True)
        for puni, mat in partners.items():
            M = np.asarray(mat, dtype=float)
            if M.ndim != 2:
                continue
            rows, cols = M.shape
            q = np.full((rows, cols), -128, dtype=np.int8)
            m = ~np.isnan(M)
            q[m] = np.round(np.clip(M[m], -3.0, 3.0) / 3.0 * 127.0).astype(np.int8)
            with open(d / f"{coords}__{puni}.bin", "wb") as f:
                f.write(struct.pack("<HH", rows, cols))
                f.write(q.tobytes())
            n_pair += 1
        n_idr += 1
        if n_idr % 100 == 0:
            print(f"[intermaps] {n_idr} domains / {n_pair} pairs ...", flush=True)

    print(f"[intermaps] done: {n_pair} binary per-pair maps ({n_idr} domains)", flush=True)


if __name__ == "__main__":
    main()
