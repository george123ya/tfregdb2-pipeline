"""Stage 2: run blastp queries against the local human-proteome database.

Wraps the NCBI BLAST+ `blastp` binary. Same parameter set as the original
notebook (cell 11) — evalue 1000 with strict seg/comp filters off so short
effector domains aren't masked away. The output is tabular outfmt 6 with
the 12 standard columns; downstream projection only needs a subset but the
full table is kept for debugging.

Requires `blastp` to be on $PATH (apt: `ncbi-blast+`, brew: `blast`).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

BLAST_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore"
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", type=Path, required=True)
    ap.add_argument("--db", type=Path, required=True,
                    help="Path stem of the BLAST DB (no .phr extension).")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--evalue", type=float, default=1000)
    args = ap.parse_args()

    if not shutil.which("blastp"):
        print("ERROR: blastp not on $PATH — install ncbi-blast+ first.",
              file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "blastp",
        "-query", str(args.query),
        "-db", str(args.db),
        "-out", str(args.out),
        "-evalue", str(args.evalue),
        "-seg", "no",
        "-comp_based_stats", "F",
        "-max_target_seqs", "1",
        "-num_threads", str(args.threads),
        "-outfmt", f"6 {BLAST_FIELDS}",
    ]
    print(f"$ {' '.join(cmd)}")
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        print(f"blastp failed with exit code {proc.returncode}", file=sys.stderr)
        return proc.returncode

    # Quick sanity report
    with args.out.open() as fh:
        n_hits = sum(1 for _ in fh)
    print(f"\n{n_hits} hit rows → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
