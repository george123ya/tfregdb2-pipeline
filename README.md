# TFRegDB 2.0 pipeline

The curation and annotation pipeline behind TFRegDB 2.0, a database of human
transcription factor effector domains (activation, repression, bifunctional)
with structural, evolutionary, sequence, and variant annotation.

This repository holds the code that turns raw curated tables and public datasets
into the per-TF records served by the website. The website itself is a separate
repository (github.com/george123ya/tfregdb2).

## Pipeline

Each numbered folder is one stage. Stages 06-12 are independent branches that all
feed the per-TF records assembled in stage 13.

```
raw tables (curated domains, TFRegDB1, HT screens, Lambert, CIS-BP)
      |
      v
 03 humanization (BLAST to canonical human protein) -> Unified_Valid_Domains
      |
      v
 04/05 isoform dataset + domain projection across isoforms
      |
      v
 13 assemble per-TF records -> curation overlays -> aggregates
      ^   (branches below run in parallel and feed the per-TF records)
      |-- 06 structure        (AlphaFold, DSSP, pLDDT, PAE)      [HPC]
      |-- 07 predictors       (ADpred, PADDLE)
      |-- 08 features         (SLiMs / ELM)
      |-- 09 genomic          (fastCDS exon-intron structure)
      |-- 10 conservation     (phyloP, ConSurf)                 [HPC]
      |-- 11 interactions     (finches IDR)                     [HPC]
      |-- 12 variants         (ClinVar, COSMIC, gnomAD)
```

| Stage | What it does | Where |
| --- | --- | --- |
| 01 curation | PubMed query + LLM triage + manual review (placeholder) | local |
| 02 ht_screens | Parse the four high-throughput screens | local |
| 03 humanization | BLAST curated domains onto the canonical human protein | local |
| 04 isoforms | Isoform repertoire from APPRIS / Ensembl / UniProt | local |
| 05 projection | Project domains onto all isoforms; add CIS-BP DBDs | local |
| 06 structure | AlphaFold models, pLDDT, PAE, secondary structure | HPC + local |
| 07 predictors | ADpred and PADDLE activation-domain tracks | local |
| 08 features | Short linear motifs (ELM) | local |
| 09 genomic | Exon-intron structure with fastCDS | local |
| 10 conservation | phyloP tracks and per-isoform ConSurf | HPC + local |
| 11 interactions | finches IDR interaction maps | HPC + local |
| 12 variants | ClinVar / COSMIC / gnomAD, projected to domains | local |
| 13 assemble | Build per-TF records, apply curation, build aggregates | local |

The three HPC stages (structure, ConSurf, finches) ran on an SGE cluster; their
job scripts and environment are described in `scc/README.md`.

## Running it

1. Create the environment and install external tools: see `env/`.
2. Obtain the raw inputs listed in `data/DATA_SOURCES.md`.
3. Set paths and thresholds in `config/config.yaml`.
4. Work through the stages in order; `run_all.sh` documents the sequence. Several
   stages are long-running or run on the cluster, so the pipeline is meant to be
   run stage by stage rather than end to end in one command.

## Notes

- Stage 01 (automated literature curation) is a placeholder; curation is manual
  for now, with PMIDs carried through the tables.
- Biophysical property tracks (hydropathy, net charge, low complexity, phase
  separation) are computed in the browser at view time and live in the website
  repository, not here.
- External tools and their versions are listed in `env/EXTERNAL_TOOLS.md`.
