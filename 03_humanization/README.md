# tfregdb-blast

Reproducible pipeline that takes the **raw curated effector-domain tables**
(`HumanizedTable`, `TFRegDB1`, and the new effector csv) and projects every
domain coordinate onto the canonical human Ensembl protein via BLASTp.

The output is the `Unified_Valid_Domains.xlsx` schema downstream consumers
need — specifically the `Final_Human_UNIPROT_ID`, `Final_Human_DomainCoords`,
and `Final_Human_Domain_Sequence` columns that `build_tf_mocks.py` reads.

This replaces the exploratory `curate_table.ipynb` notebook with a scripted
flow you can run end-to-end and version.

## What it does, in order

```
input tables  ──►  01_prepare_queries.py  ──►  queries.fasta
                                                    │
                       human_proteome.fasta  ──►  02_build_db.py  ──►  custom_db/
                                                    │
                                              blastp (subprocess)  ──►  blast_results.tsv
                                                    │
                                          03_project_canonical.py
                                                    │
                                                    ▼
                                  Unified_Valid_Domains.xlsx
```

Each step is idempotent and reads / writes plain files; intermediate state
lives in `data/output/` so reruns skip already-done work.

## Inputs

Drop the following files into `data/input/`:

| File | Source |
|---|---|
| `TF_completeIDs_with_ensembl_canonical.csv` | Lambert 2018 + Ensembl APPRIS picker output (the "canonical record" lookup table the rest of the pipeline keys on) |
| `HumanizedTable.xlsx` | Cleaned `UpdatedTable.xlsx`, "Hoja 1" sheet, with curated AD/RD calls + PMIDs |
| `TFRegDB1.xlsx` | Same UpdatedTable, "TFRegDB1" sheet |
| `TF_effector_data_with_derivadas_and_alignment_ideal.csv` | The per-isoform effector csv with SW-alignment quality |
| `human_proteome.fasta` | UniProt Swiss-Prot human reference (download from UniProt) |

## Running

```bash
# One-time prep (build the BLAST DB from the human proteome FASTA)
make db

# Build queries (FASTA of every effector domain sequence, headered by row ID)
make queries

# Run BLAST
make blast

# Project hits → canonical, write Unified_Valid_Domains.xlsx
make project

# All of the above
make all
```

## Outputs

Everything lands in `data/output/`:

- `queries.fasta` — one record per source row, header `>{source}_row{N}|{gene}|{uniprot}`
- `blast_results.tsv` — raw `blastp -outfmt 6` table
- `Unified_Valid_Domains.xlsx` — final coordinate-projected curation table
- `blast_failures.tsv` — rows that didn't BLAST cleanly (<95% identity, no hit, etc.)

## Why the rewrite

The notebook flow worked but: (1) the BLAST step was buried inside a `subprocess.run` cell that's hard to rerun on partial data, (2) the canonical-projection logic was duplicated for `HumanizedTable` and `TFRegDB1`, (3) intermediate state (the `df_final` variables) lived only in kernel memory and was easy to lose. Splitting into three idempotent scripts + a Makefile fixes all three and lets you run any single stage with `make <stage>`.

## License

MIT (or pick your own at `git init` time).
