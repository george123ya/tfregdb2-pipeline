#!/usr/bin/env bash
# Stage order for the pipeline. Most stages are long-running, need external data,
# or run on the cluster, so run them one at a time and check the output rather
# than executing this top to bottom. This file is the map, not a one-click build.
set -euo pipefail

echo "01 curation      - manual for now (see 01_curation/README.md)"

echo "02 ht_screens    - parse the four high-throughput screens"
#   python 02_ht_screens/assemble_from_master.py

echo "03 humanization  - BLAST curated domains onto canonical human proteins"
#   cd 03_humanization && make            # enrich -> fetch-animals -> db -> queries -> blast -> project

echo "04/05 projection - isoform dataset + project domains + CIS-BP DBDs"
#   python 05_projection/map_canonical_to_isoforms.py

echo "13 assemble      - build per-TF records from the projected tables"
#   python 13_assemble/build_tf_records.py

echo "-- annotation branches (independent; feed the per-TF records) --"
echo "06 structure     - HPC folding + local pLDDT/PAE/secondary  (see scc/README.md)"
echo "07 predictors    - python 07_predictors/extract_adpred.py ; extract_paddle.py"
echo "08 features      - python 08_features/extract_slims_elm.py"
echo "09 genomic       - fastCDS via 09_genomic/wire_isoform_genomics.py"
echo "10 conservation  - local phyloP + HPC ConSurf  (see scc/README.md)"
echo "11 interactions  - HPC finches extract + local convert_interactions.py"
echo "12 variants      - 12_variants: fetch -> build_variant_db.py -> project_variants_to_tfs.py"

echo "13 assemble      - apply curation overlays, then build aggregates"
#   python 13_assemble/overlays/recuration_apply.py
#   python 13_assemble/build_api_index.py
#   python 13_assemble/compute_constraint_stats.py
#   python 13_assemble/build_downloads.py
