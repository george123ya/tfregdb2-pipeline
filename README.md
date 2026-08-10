<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/logo-dark.svg">
    <img alt="TFRegDB v2.0" src="assets/logo-light.svg" width="300">
  </picture>
</p>

# TFRegDB v2.0 pipeline

The curation and annotation pipeline behind TFRegDB v2.0, a database of human
transcription factor effector domains (activation, repression, bifunctional)
with structural, evolutionary, sequence, and variant annotation.

This repository holds the code that turns raw curated tables and public datasets
into the per-TF records served by the website. The website itself is a separate
repository (github.com/george123ya/tfregdb2).

## Pipeline

Each numbered folder is one stage. Stages 06-11 are independent annotation
branches that all feed the per-TF records assembled in stage 12.

```
raw tables: curated domains, TFRegDB1, HT screens, Lambert, CIS-BP
      |
      v
 03 canonical mapping         BLAST curated domains to the canonical human protein
      |
      v
 04 isoforms + 05 projection  project domains across all isoforms; add CIS-BP DBDs
      |
      |  annotation branches (independent, run in parallel):
      +-- 06 structure         AlphaFold, DSSP, pLDDT, PAE         [HPC]
      +-- 07 predictors        ADpred, PADDLE
      +-- 08 features          short linear motifs (ELM)
      +-- 09 genomic           fastCDS exon-intron structure
      +-- 10 conservation      phyloP, ConSurf                     [HPC]
      +-- 11 variants          ClinVar, COSMIC, gnomAD
      |
      v
 12 assemble                  one record per TF, curation overlays, build aggregates
```

| Stage | What it does | Where |
| --- | --- | --- |
| 01 curation | PubMed query + LLM triage + manual review (placeholder) | local |
| 02 ht_screens | Parse the four high-throughput screens | local |
| 03 canonical mapping | BLAST curated domains onto the canonical human protein | local |
| 04 isoforms | Isoform repertoire from APPRIS / Ensembl / UniProt | local |
| 05 projection | Project domains onto all isoforms; add CIS-BP DBDs | local |
| 06 structure | AlphaFold models, pLDDT, PAE, secondary structure | HPC + local |
| 07 predictors | ADpred and PADDLE activation-domain tracks | local |
| 08 features | Short linear motifs (ELM) | local |
| 09 genomic | Exon-intron structure with fastCDS | local |
| 10 conservation | phyloP tracks and per-isoform ConSurf | HPC + local |
| 11 variants | ClinVar / COSMIC / gnomAD, projected to domains | local |
| 12 assemble | Build per-TF records, apply curation, build aggregates | local |

The two HPC stages (structure and ConSurf) ran on an SGE cluster; their
job scripts and environment are described in `scc/README.md`.

## What you need to start

Only the mapping stage (03) needs external input files; later stages read what
earlier ones produce. To run from stage 03 onward you need:

- **Curated tables** (distributed with the release): `UpdatedTable.xlsx` (or
  `HumanizedTable.xlsx` + `TFRegDB1.xlsx`) with the effector-domain reports and
  the high-throughput study sheets.
- **`TF_completeIDs_with_ensembl_canonical.csv`** — the gene to canonical
  Ensembl protein lookup.
- **`human_proteome.fasta`** — UniProt Swiss-Prot human, for the BLAST database.
- **`TF_dbd_data_.xlsx`** — CIS-BP DNA-binding domains.

The public reference datasets each stage pulls (APPRIS, Ensembl, CIS-BP,
AlphaFold, ELM, phyloP, ClinVar/COSMIC/gnomAD) and where to get them are listed
in `data/DATA_SOURCES.md`.

## Running it

1. Create the environment and install external tools: see `env/`.
2. Obtain the raw inputs listed in `data/DATA_SOURCES.md`.
3. Set paths and thresholds in `config/config.yaml`.
4. Work through the stages in order; `run_all.sh` documents the sequence. Several
   stages are long-running or run on the cluster, so the pipeline is meant to be
   run stage by stage rather than end to end in one command.

## Updating the database

When the curated table changes (for example when the remaining rows are curated),
point `RECUR_XLSX` and the source tables at the new files and run:

    bash update.sh

It re-maps the domains, rebuilds the per-TF records, re-applies the curation
overlay, rebuilds the aggregates, and prints the publish step. The annotation
branches (06-11) only need re-running if new TFs or isoforms were added.

## Notes

- Stage 01 (automated literature curation) is a placeholder; curation is manual
  for now, with PMIDs carried through the tables.
- Biophysical property tracks (hydropathy, net charge, low complexity, phase
  separation) are computed in the browser at view time and live in the website
  repository, not here.
- External tools and their versions are listed in `env/EXTERNAL_TOOLS.md`.
