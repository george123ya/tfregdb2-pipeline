# 11 - Interactions

IDR-mediated interaction maps from finches. `s1`-`s5` collect STRING partners,
resolve IDs, and fetch and filter PDBs. The heavy epsilon-matrix extraction runs
on the cluster (`extract_intermaps.qsub`); `convert_interactions.py` turns the
summary into per-TF interaction JSON.
