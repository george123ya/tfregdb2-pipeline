# 12 - Variants

ClinVar, COSMIC v103 and gnomAD v4.1 variants. Each source is fetched and parsed
to a TSV, loaded into a SQLite database (`build_variant_db.py`), then sliced by
each canonical CDS and projected to protein position (`project_variants_to_tfs.py`).
Intronic variants are added for the genomic view.
