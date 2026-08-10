"""Run PADDLE-noSS (Sanborn et al. 2021, eLife) across every curated TF to
produce a per-residue predicted activation-domain (AD) track.

PADDLE is a 10-CNN ensemble that scores 53-aa tiles for acidic-AD activation
(yeast MPRA-trained, transfers to human). We use the **noSS** variant — it
needs only sequence (no PSIPRED / IUPRED2), so no BLAST DB or disorder
pre-compute is required, and it is "nearly as accurate as PADDLE" (the paper).

Speed trick (same as extract_adpred.py): rather than the library's per-tile
calls, we stack every 53-aa tile of a protein into one (N,53,20) batch and run
each model's serving signature once. Loading via `tf.saved_model.load` +
`signatures['serving_default']` bypasses Keras entirely, so the TF-2.2-trained
SavedModels run unmodified under the env's TF 2.20.

Per-tile Z-scores (length L-52) are averaged across the 10 splits, the model
offset (2.4845) is added, and a 9-aa sliding mean is applied — matching
`PADDLE.predict_protein`. Each residue r takes the score of the tile centred
on it (tile start = clamp(r-26, 0, L-53)), giving a length-L per-residue track.
A predicted Z-score ≳ 4 marks a strong acidic AD (paper threshold).

Output: public/mock/paddle/{uniprot}.json — packed float array length L,
2-decimal Z-scores. [symbol].astro renders it as a track alongside ADpred.

Run from repo root:
    PADDLE_DIR=/home/goguxor/Desktop/PADDLE \
      /home/goguxor/miniconda3/envs/adpred/bin/python scripts/extract_paddle.py
Test one: PADDLE_ONLY=TP53 ... scripts/extract_paddle.py
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")  # CPU — tiny CNN, dodges driver issues

import numpy as np  # noqa: E402
import tensorflow as tf  # noqa: E402

PADDLE_DIR = Path(os.environ.get("PADDLE_DIR", "/home/goguxor/Desktop/PADDLE"))
TFS_DIR = Path("public/mock/tfs")
OUT_DIR = Path("public/mock/paddle")
AAS = list("ACDEFGHIKLMNPQRSTVWY")
AA_IDX = {a: i for i, a in enumerate(AAS)}
TILE = 53
OFFSET = 2.4845441766212284  # PADDLE model offset (Z-score shift)
SMOOTH = 9


def load_signatures() -> list:
    sigs = []
    for i in range(10):
        mp = PADDLE_DIR / "models" / f"PADDLE_noSS_split{i}.model"
        if not mp.exists():
            sys.exit(f"missing model: {mp}")
        m = tf.saved_model.load(str(mp))
        sigs.append(m.signatures["serving_default"])
    return sigs


def onehot_tiles(seq: str) -> np.ndarray:
    """(N,53,20) one-hot for every 53-aa tile; N = len(seq)-52."""
    n = len(seq) - TILE + 1
    out = np.zeros((n, TILE, len(AAS)), dtype=np.float32)
    for i in range(n):
        for k, aa in enumerate(seq[i:i + TILE]):
            j = AA_IDX.get(aa)
            if j is not None:
                out[i, k, j] = 1.0
    return out


def predict_protein(seq: str, sigs: list) -> np.ndarray:
    """Per-residue smoothed Z-score, length == len(seq). Empty if seq < 53."""
    L = len(seq)
    if L < TILE:
        return np.array([], dtype=np.float32)
    batch = tf.constant(onehot_tiles(seq))
    cols = []
    for sig in sigs:
        # Each split's SavedModel has its own input tensor name
        # (conv1d_634_input, conv1d_643_input, …) — read it per-signature.
        in_key = list(sig.structured_input_signature[1].keys())[0]
        out = sig(**{in_key: batch})
        cols.append(list(out.values())[0].numpy().reshape(-1))
    tile_z = np.mean(np.stack(cols, axis=0), axis=0) + OFFSET  # length L-52

    if SMOOTH > 1:
        half = (SMOOTH - 1) / 2
        tile_z = np.array([
            tile_z[max(0, int(n - half)): int(n + half) + 1].mean()
            for n in range(len(tile_z))
        ])

    # Map each residue to the tile centred on it.
    n_tiles = len(tile_z)
    res = np.empty(L, dtype=np.float32)
    for r in range(L):
        res[r] = tile_z[min(max(0, r - TILE // 2), n_tiles - 1)]
    return res


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    only = (os.environ.get("PADDLE_ONLY") or "").strip()
    tf_files = sorted(TFS_DIR.glob("*.json"))
    if only:
        wanted = {s.strip().upper() for s in only.split(",")}
        tf_files = [f for f in tf_files if f.stem.upper() in wanted]

    print(f"[paddle] loading 10 noSS models from {PADDLE_DIR}/models …", flush=True)
    sigs = load_signatures()
    print(f"[paddle] running across {len(tf_files)} TFs"
          + (f" (filter: {only})" if only else ""), flush=True)

    n_ok = n_short = n_skip = n_resume = 0
    t0 = time.time()
    for i, fp in enumerate(tf_files, 1):
        d = json.loads(fp.read_text())
        seq = (d.get("sequence") or "").upper()
        upro = (d.get("uniprot") or "").strip()
        if not seq or not upro:
            n_skip += 1
            continue
        out_path = OUT_DIR / f"{upro}.json"
        if out_path.exists():
            n_resume += 1
            continue
        scores = predict_protein(seq, sigs)
        if scores.size == 0:
            n_short += 1
            continue
        out_path.write_text(
            json.dumps([round(float(s), 2) for s in scores], separators=(",", ":"))
        )
        n_ok += 1
        if i % 50 == 0:
            dt = time.time() - t0
            rate = i / dt
            eta = (len(tf_files) - i) / rate
            print(f"  [{i}/{len(tf_files)}] {n_ok} ok · {dt:.0f}s · "
                  f"{rate:.1f} TFs/s · ETA {eta/60:.1f} min", flush=True)
    dt = time.time() - t0
    print(f"[paddle] canonical: wrote {n_ok:,} files in {dt/60:.1f} min · "
          f"resumed:{n_resume} short(<53aa):{n_short} skip:{n_skip} → {OUT_DIR}",
          flush=True)

    # ── Per-isoform pass: PADDLE in each isoform's own frame, keyed by iso
    # NAME (e.g. RARA-208.json) — matches the conservation / adpred files so the
    # track reframes when an isoform is selected. noSS needs only the sequence,
    # so (unlike ADpred) there is no secondary-structure dependency.
    n_iok = n_ishort = n_ires = 0
    t_iso = time.time()
    for i, fp in enumerate(tf_files, 1):
        d = json.loads(fp.read_text())
        for iso in d.get("isoforms") or []:
            name = iso.get("name")
            iseq = (iso.get("sequence") or "").upper()
            if not name or not iseq:
                continue
            out_path = OUT_DIR / f"{name}.json"
            if out_path.exists():
                n_ires += 1
                continue
            scores = predict_protein(iseq, sigs)
            if scores.size == 0:
                n_ishort += 1
                continue
            out_path.write_text(
                json.dumps([round(float(s), 2) for s in scores], separators=(",", ":"))
            )
            n_iok += 1
        if i % 100 == 0:
            el = time.time() - t_iso
            print(f"  [iso {i}/{len(tf_files)}] {n_iok} written · {el/60:.1f} min",
                  flush=True)
    print(f"[paddle] isoforms: wrote {n_iok:,} files in {(time.time()-t_iso)/60:.1f} min · "
          f"resumed:{n_ires} short(<53aa):{n_ishort} → {OUT_DIR}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
