#!/usr/bin/env bash
# Rebuild the database from an updated curated table (e.g. when Sheyla finishes
# the remaining rows) and republish. Point RECUR_XLSX / the source tables at the
# new files first (see config/config.yaml).
set -euo pipefail
cd "$(dirname "$0")"

echo "== 1. map curated domains to canonical human proteins =="
( cd 03_canonical_mapping && make )                 # -> Unified_Valid_Domains

echo "== 2. project across isoforms + rebuild per-TF records =="
python3 05_projection/map_canonical_to_isoforms.py
python3 12_assemble/build_tf_records.py

echo "== 3. apply curation overlay + reconcile counts =="
python3 12_assemble/overlays/recuration_apply.py
python3 12_assemble/overlays/remove_dbdonly.py
python3 12_assemble/overlays/relabel_cisbp3.py
python3 12_assemble/overlays/rebuild_index.py --write

echo "== 4. rebuild aggregates =="
python3 12_assemble/build_api_index.py
python3 12_assemble/compute_constraint_stats.py
python3 12_assemble/build_downloads.py

# Re-run the annotation branches (06-11) only if new TFs or isoforms were added;
# a curation-only update of existing TFs does not need them.

echo "== 5. publish: R2_BUCKET=tfregdb2-data bash 12_assemble/upload_r2.sh, then push the site repo =="
