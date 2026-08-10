# 06 - Structure

Predicted structures and per-residue confidence for every canonical protein and
alternative isoform.

Most models come from the AlphaFold Database. The 475 sequences it lacks (15
canonical proteins + 460 isoforms) were folded in house with AlphaFold 2.3.1 on
the cluster (`af2_features.qsub` -> `af2_inference.qsub`), and their PAE matrices
extracted with `extract_pae.qsub` / `extract_pae.py`. See `../scc/README.md`.

Local steps:

- `download_isoform_pdbs.py`, `download_isoform_pae.py` - fetch models and PAE
  from the AlphaFold Database.
- `extract_plddt.py` - per-residue pLDDT from PDB B-factors.
- `extract_secondary.py`, `extract_secondary_dssp.py` - secondary structure with
  DSSP (`mkdssp 4`).
- `fix_x_residues.py` - substitute unknown (X) residues before folding.
- `rank_from_unrelaxed.py` - rank models by mean pLDDT when the AMBER relax runs
  out of memory on very long chains.
