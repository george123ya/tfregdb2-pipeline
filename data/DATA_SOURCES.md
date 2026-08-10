# Data sources

Raw inputs the pipeline consumes. They are not stored in this repository; obtain
them as noted and point `config/config.yaml` at their location. The curated
tables are distributed separately.

## Curated tables (distributed with the release)

| File | Contents |
| --- | --- |
| `UpdatedTable.xlsx` / `HumanizedTable.xlsx` | Curated effector-domain reports: gene, species, domain coordinates, type, assay, activity, PMID. Includes the high-throughput study sheets. |
| `TFRegDB1.xlsx` | The 2021 compendium (TF name, UniProt ID, sequence, PMID). |
| `TF_completeIDs_with_ensembl_canonical.csv` | Canonical record lookup: gene, Ensembl gene/transcript/translation IDs, protein length, APPRIS annotation, canonical flag, sequence, UniProt ID. |
| `TF_dbd_data_.xlsx` | DNA-binding-domain coordinates and Pfam families (from CIS-BP). |
| `PTM_aggregated.xlsx` | Per-residue PTM sites, keyed by translation ID. |
| `human_proteome.fasta` | UniProt Swiss-Prot human reference; source for the BLAST database. |

## Public datasets (download)

| Source | Where | Used in |
| --- | --- | --- |
| Lambert et al. TF catalogue (DatabaseExtract v1.01) | published supplement | TF universe, families |
| High-throughput screens: Tycko 2020 (33326746), Tycko 2024 (39487265), Del Rosso 2023 (37020022), Alerasool 2022 (35016035) | each study's supplement | stage 02 |
| APPRIS 2025_07.v50, Ensembl REST 115 | apprisws.bioinfo.cnio.es, rest.ensembl.org | stage 04 |
| CIS-BP 3.0 | cisbp.ccbr.utoronto.ca | stage 05 (DBDs) |
| AlphaFold Database | alphafold.ebi.ac.uk | stage 06 |
| ELM classes (release 1.4) | elm.eu.org | stage 08 |
| GENCODE v49 GTF | gencodegenes.org | stage 09 |
| phyloP cactus241way + hg38.phyloP470way | UCSC | stage 10 |
| UniRef90 | UniProt | stage 10 (ConSurf) |
| ClinVar variant_summary, COSMIC v103, gnomAD v4.1 | NCBI, COSMIC, gnomAD | stage 12 |
