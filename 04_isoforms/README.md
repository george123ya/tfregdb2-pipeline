# 04 - Isoform dataset

The isoform repertoire per gene from APPRIS and the Ensembl REST API, linked to
UniProt (reviewed Swiss-Prot first, TrEMBL fallback). The canonical protein is
chosen by Ensembl canonical, then APPRIS principal, then longest. The build of
`TF_completeIDs_with_ensembl_canonical.csv` is external; see `data/DATA_SOURCES.md`.
