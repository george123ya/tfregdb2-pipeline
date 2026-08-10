# External tools and data versions

Versions used to build the current release. Python packages are in
`requirements.txt`. Several scripts still hardcode local paths and conda
environment names (for example the DSSP and PADDLE binaries); parameterize
those through `config/config.yaml` before running elsewhere.

## Tools

| Tool | Version | Stage | Notes |
| --- | --- | --- | --- |
| Biopython | 1.87 | 01, 05, 06 | Entrez queries, pairwise alignment, DSSP wrapper |
| NCBI BLAST+ (blastp, makeblastdb) | unpinned | 03 | run permissive: `-evalue 1000 -seg no -comp_based_stats F -max_target_seqs 1` |
| GPT-5-mini | gpt-5-mini-2025-08-07 | 01 | article triage and extraction (not yet scripted here) |
| metapub | 0.6.4 | 01 | full-text retrieval |
| AlphaFold | 2.3.1 (CUDA 11.8) | 06 | in-house folding on HPC; most models from the AlphaFold Database |
| DSSP (mkdssp) | 4.x | 06 | secondary structure |
| ADpred | Erijman 2020 | 07 | CNN activation-domain predictor (Keras/TensorFlow) |
| PADDLE-noSS | Sanborn 2021 | 07 | 10-model ensemble, TensorFlow 2.20 SavedModels |
| ELM database | release 1.4 | 08 | short linear motif classes |
| fastCDS | 2.2.0 | 09 | exon-intron structure (C++ tool) |
| ConSurf (standalone) | HMMER algorithm | 10 | per-isoform conservation, on HPC |

## Reference datasets

| Dataset | Version | Used for |
| --- | --- | --- |
| Lambert et al. TF catalogue | DatabaseExtract v1.01 (2018) | the human TF universe and family assignment |
| APPRIS | 2025_07.v50 (GRCh38) | principal isoform / canonical selection |
| Ensembl REST API | release 115 | isoform sequences and identifiers |
| UniProt | Swiss-Prot (TrEMBL fallback) | canonical proteins, cross-species sequences |
| CIS-BP | 3.0 | DNA-binding domain coordinates and Pfam families |
| GENCODE | v49 | fastCDS genomic index |
| phyloP | cactus241way + hg38.phyloP470way (UCSC) | nucleotide conservation |
| UniRef90 | current | ConSurf homolog search |
| ClinVar | variant_summary | pathogenic / benign variants |
| COSMIC | v103 | somatic variants |
| gnomAD | v4.1 (exomes + genomes) | population variants |
