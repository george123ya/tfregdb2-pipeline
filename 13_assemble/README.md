# 13 - Assemble

`build_tf_mocks.py` unifies the canonical record, effector domains, DBDs, PTMs and
structure assets into one record per TF, plus the family index. The scripts in
`overlays/` apply the manual re-curation (`recuration_apply.py`,
`remove_dbdonly.py`, `relabel_cisbp3.py`, reconciled by `rebuild_index.py`).
`build_api_index.py`, `compute_constraint_stats.py` and `build_downloads.py`
produce the query index, statistics and download tables.
