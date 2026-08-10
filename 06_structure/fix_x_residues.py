#!/usr/bin/env python3
"""Resolve ambiguous X (UNK) residues in the fold-target FASTAs so AlphaFold's
Amber relaxation works (UNK has no atoms → relax pre-check crashes).

These X's come from non-curated TrEMBL/Ensembl isoform translations:
  • A leading X (position 1) is the unresolved START residue → set to M
    (translation initiates at Met; verified: every X-leading isoform's canonical
    starts with M).
  • A rare internal X → A (neutral placeholder; only 1 case, A0A5K1VW50@428).

Edits fold_targets/fastas/*.fasta IN PLACE (length preserved, so domain / PAE
coords are unaffected), annotates the header, and writes a provenance TSV.
NOTE: the sequence changes, so any already-computed features.pkl / MSAs for the
old X sequences are stale — regenerate features before re-folding.

Run:  python3 scripts/fix_x_residues.py
"""
import glob, os

FASTA_DIR = "fold_targets/fastas"
LOG = "fold_targets/x_substitutions.tsv"


def main():
    rows = []
    for fp in sorted(glob.glob(os.path.join(FASTA_DIR, "*.fasta"))):
        lines = open(fp).read().splitlines()
        header = lines[0]
        seq = "".join(l.strip() for l in lines[1:] if l.strip())
        if "X" not in seq:
            continue
        acc = os.path.basename(fp)[:-6]
        subs = []
        chars = list(seq)
        for i, c in enumerate(chars):
            if c != "X":
                continue
            repl = "M" if i == 0 else "A"
            chars[i] = repl
            subs.append(f"X{i + 1}{repl}")
        newseq = "".join(chars)
        # rewrite fasta (60-col), annotate header
        h = header if "Xsub=" in header else f"{header} Xsub={','.join(subs)}"
        with open(fp, "w") as fh:
            fh.write(h + "\n")
            for j in range(0, len(newseq), 60):
                fh.write(newseq[j : j + 60] + "\n")
        rows.append((acc, len(newseq), ";".join(subs)))

    with open(LOG, "w") as fh:
        fh.write("accession\tlength\tsubstitutions\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")
    lead = sum(1 for _, _, s in rows if s.startswith("X1M"))
    print(f"fixed {len(rows)} FASTAs ({lead} leading X1->M, {len(rows) - lead} with an internal X->A)")
    print(f"provenance → {LOG}")
    # sanity: confirm no X remains
    left = sum(
        1
        for fp in glob.glob(os.path.join(FASTA_DIR, "*.fasta"))
        if "X" in "".join(l.strip() for l in open(fp) if not l.startswith(">"))
    )
    print(f"FASTAs still containing X: {left}")


if __name__ == "__main__":
    main()
