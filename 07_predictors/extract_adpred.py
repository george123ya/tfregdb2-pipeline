"""Run ADpred (Erijman 2020) across every curated TF, replacing the
9aaTAD heuristic in `biophys.ts` with the real CNN trained on yeast
MPRA activation-domain data.

Why this is faster than ADpred's own `predict()`:
  The library calls `ADPred.predict()` once per 30-residue window. Each
  call goes through Keras's full call overhead (model graph traversal,
  device pinning, etc.) ~22 ms each. On a 393-aa protein that's 363
  windows × 22 ms = 8 s pure overhead, with ~25 s total wall time.

  We stack all windows for a protein into a single (N, 30, 23, 1)
  tensor and submit one `predict()` call. Keras then vectorises the
  CNN forward pass across the batch and the per-protein cost drops to
  the actual compute time (~1-3 s).

DSSP secondary structure (extract_secondary.py) is reused as ADpred's
ss input — saves a PSIPRED webserver call per TF and uses the same
high-quality predictions we already paid for. C → '-' relabel matches
ADpred's training vocabulary {E, H, -}.

Output: `public/mock/adpred/{uniprot}.json` as a packed float array
of length N (per-residue P(AD)). [symbol].astro reads this and falls
back to the in-app `predictADs` heuristic when missing.
"""

from __future__ import annotations

import json
import os
import sys
import time
import warnings
from pathlib import Path

# Quiet the TF stack — we don't have GPU and the dlerror noise drowns the log.
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
warnings.filterwarnings("ignore")

import numpy as np
from adpred.ADpred import ADPred, calculate_psipred  # noqa: E402  (after env vars)
from adpred.lib.utils import aa as AA_LEX, ss as SS_LEX, make_ohe  # noqa: E402


TFS_DIR = Path("public/mock/tfs")
SS_DIR = Path("public/mock/secondary")
OUT_DIR = Path("public/mock/adpred")
WINDOW = 30
PAD = 15  # ADpred pads the protein with 15 'G' / '-' at each end


def predict_batched(seq: str, struct: str) -> np.ndarray:
    """Returns per-residue P(AD), length == len(seq).

    `struct` must be the same length as seq, using ADpred's vocab {E,H,-}.
    """
    if len(seq) != len(struct):
        raise ValueError(f"seq/struct length mismatch: {len(seq)} vs {len(struct)}")
    if any(c not in AA_LEX for c in seq):
        # ADpred's lexicon is the standard 20. Selenocysteine (U) and
        # pyrrolysine (O) are rare in TFs but break make_ohe. Replace
        # any out-of-vocab character with 'X' → drop the protein.
        return np.array([], dtype=np.float32)

    padded_seq = ("G" * PAD) + seq + ("G" * PAD)
    padded_ss = ("-" * PAD) + struct + ("-" * PAD)
    ohe = make_ohe(padded_seq, padded_ss)  # shape (L+30, 23)
    n_windows = len(seq)
    # Build (N, 30, 23, 1) by sliding the window.
    batch = np.stack([ohe[i:i + WINDOW] for i in range(n_windows)], axis=0)
    batch = batch.reshape(n_windows, WINDOW, 23, 1).astype(np.float32)
    preds = ADPred.predict(batch, verbose=0)
    return preds.reshape(-1).astype(np.float32)


def run_one(seq: str, ss_raw: str, out_path: Path) -> str:
    """Run ADpred for one protein, write out_path. Returns a status string."""
    ss = ss_raw.replace("C", "-")
    if len(ss) != len(seq):
        ss = ss[: len(seq)] if len(ss) > len(seq) else ss + "-" * (len(seq) - len(ss))
    try:
        scores = predict_batched(seq, ss)
    except Exception as e:  # noqa: BLE001
        return f"err:{e}"
    if scores.size == 0:
        return "skip"
    out_path.write_text(json.dumps([round(float(s), 3) for s in scores], separators=(",", ":")))
    return "ok"


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Optional symbol filter for testing: ADPRED_ONLY=RARA python ... extract_adpred.py
    only = (os.environ.get("ADPRED_ONLY") or "").strip()
    tf_files = sorted(TFS_DIR.glob("*.json"))
    if only:
        wanted = {s.strip().upper() for s in only.split(",")}
        tf_files = [f for f in tf_files if f.stem.upper() in wanted]
    print(f"[extract_adpred] running ADpred across {len(tf_files)} TFs"
          + (f" (filter: {only})" if only else ""))

    n_ok = n_skip_seq = n_skip_ss = n_resume = 0
    t_start = time.time()
    for i, fp in enumerate(tf_files, 1):
        tf = json.loads(fp.read_text())
        seq = (tf.get("sequence") or "").upper()
        upro = (tf.get("uniprot") or "").strip()
        if not seq or not upro:
            n_skip_seq += 1
            continue
        out_path = OUT_DIR / f"{upro}.json"
        if out_path.exists():
            # Resume mode: if the disconnect killed a prior run, pick up
            # where it left off. Existing files are trusted; rm them
            # manually to force a re-run.
            n_resume += 1
            continue
        ss_path = SS_DIR / f"{upro}.json"
        if not ss_path.exists():
            n_skip_ss += 1
            continue
        ss_raw = json.loads(ss_path.read_text()).get("ss") or ""
        # DSSP C → ADpred '-'. Length should already match the protein
        # (extract_secondary.py clamps for off-by-one cases).
        ss = ss_raw.replace("C", "-")
        if len(ss) != len(seq):
            # Conservative clamp/pad — mismatched lengths would crash
            # make_ohe. Pad with '-' (coil) when ss is short.
            if len(ss) > len(seq):
                ss = ss[:len(seq)]
            else:
                ss = ss + "-" * (len(seq) - len(ss))
        try:
            scores = predict_batched(seq, ss)
        except Exception as e:
            print(f"  WARN {upro}: {e}")
            n_skip_seq += 1
            continue
        if scores.size == 0:
            n_skip_seq += 1
            continue
        # 3-decimal precision keeps the JSON ~3 KB per TF (vs 8 KB at full
        # float). Plenty for a 0..1 score.
        (OUT_DIR / f"{upro}.json").write_text(
            json.dumps([round(float(s), 3) for s in scores], separators=(",", ":"))
        )
        n_ok += 1
        if i % 50 == 0:
            dt = time.time() - t_start
            rate = i / dt
            eta = (len(tf_files) - i) / rate
            print(f"  [{i}/{len(tf_files)}] {n_ok} ok · {dt:.0f}s elapsed · "
                  f"{rate:.1f} TFs/s · ETA {eta/60:.1f} min")
    dt = time.time() - t_start
    print(f"  wrote {n_ok:,} ADpred files in {dt/60:.1f} min · "
          f"resumed:{n_resume} skipped seq:{n_skip_seq} ss:{n_skip_ss} · → {OUT_DIR}")

    # ── Per-isoform pass: ADpred in each isoform's own frame, keyed by iso
    # NAME (e.g. RARA-208.json), matching the conservation / conway files.
    # Needs the isoform's secondary structure (public/mock/secondary/{acc}.json);
    # isoforms that share the canonical accession reuse the canonical SS.
    n_iso_ok = n_iso_ss = n_iso_resume = 0
    t_iso = time.time()
    for fp in tf_files:
        tf = json.loads(fp.read_text())
        for iso in tf.get("isoforms") or []:
            name = iso.get("name")
            iseq = (iso.get("sequence") or "").upper()
            acc = (iso.get("uniprotIso") or "").strip()
            if not name or not iseq:
                continue
            out_path = OUT_DIR / f"{name}.json"
            if out_path.exists():
                n_iso_resume += 1
                continue
            ss_path = SS_DIR / f"{acc}.json" if acc else None
            if not ss_path or not ss_path.exists():
                n_iso_ss += 1
                continue
            ss_raw = json.loads(ss_path.read_text()).get("ss") or ""
            status = run_one(iseq, ss_raw, out_path)
            if status == "ok":
                n_iso_ok += 1
    print(f"  wrote {n_iso_ok:,} isoform ADpred files in {(time.time()-t_iso)/60:.1f} min · "
          f"resumed:{n_iso_resume} no-ss:{n_iso_ss}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
