#!/usr/bin/env bash
# Upload TFRegDB2 bulk data to Cloudflare R2 under a versioned prefix, and build
# the build-input tarball that the CI build (scripts/build_pages.sh) pulls.
#
# Requires an rclone remote named "r2" (S3-compatible, Cloudflare provider):
#   rclone config create r2 s3 provider=Cloudflare \
#     access_key_id=<R2_ACCESS_KEY_ID> secret_access_key=<R2_SECRET> \
#     endpoint=https://<ACCOUNT_ID>.r2.cloudflarestorage.com
#
# Then:  R2_BUCKET=tfregdb2-data bash scripts/upload_r2.sh
# PUBLIC_DATA_BASE for the site = https://<public-bucket-url>/<VER>  (e.g. .../v1)
set -euo pipefail
cd "$(dirname "$0")/.."

REMOTE=${R2_REMOTE:-r2}
BUCKET=${R2_BUCKET:?set R2_BUCKET (e.g. tfregdb2-data)}
VER=${R2_VERSION:-v1}
DEST="$REMOTE:$BUCKET/$VER"
FLAGS="--transfers=48 --checkers=48 --fast-list --s3-no-check-bucket"

# Refresh the precomputed Statistics aggregate (ED-vs-DBD constraint + Figure-1
# stats) from the full local variant / domain / pLDDT data — the CF build env
# only has the tarball, so this must be baked in here, not regenerated there.
echo "== computing stats_constraint.json =="
python3 scripts/compute_constraint_stats.py || echo "  (skipped — variant data not present)"
echo "== building API index (domains.json + D1 in-domain variants) =="
python3 scripts/build_api_index.py || echo "  (skipped — variant data not present)"
echo "== building bulk download files =="
python3 scripts/build_downloads.py || echo "  (skipped — source data not present)"

# 1) Build-input tarball: the dirs the build reads from disk (mock/{index, tfs,
#    slims, secondary, plddt, adpred, paddle} + stats_constraint.json).
echo "== building build_data.tar.gz =="
tar czf /tmp/build_data.tar.gz \
  public/mock/index.json public/mock/tfs public/mock/slims \
  public/mock/secondary public/mock/plddt public/mock/adpred public/mock/paddle \
  public/mock/stats_constraint.json
echo "   $(du -h /tmp/build_data.tar.gz | cut -f1)"

# 2) Runtime data the client fetches via dataUrl() (PUBLIC_DATA_BASE + path).
# mock/adpred + mock/slims are also bundled in the build tarball (canonical
# tracks read them via fs at build), but the per-isoform files (keyed by iso
# NAME, e.g. RARA-208.json) are fetched at RUNTIME by the biophysics tracks
# when an isoform is selected — so they must live as individual objects on R2.
# mock/api + mock/tfs + mock/index.json + mock/stats_constraint.json are fetched
# at RUNTIME by the dynamic API (functions/api/v1) — mock/tfs/index also feed the
# build (tarball), but the API needs them as individual R2 objects too.
# mock/conway_genomic: per-gene genomic phyloP .bin + _index.json, fetched by
#   the conservation Genomic view. mock/conway (per-isoform CDS) is still used
#   by the protein phyloP overlay in SequenceText. mock/paddle: per-isoform
#   PADDLE (none yet — canonical is baked into the build tarball; listed so the
#   _index/canonical objects exist on R2 for future per-iso files).
RUNTIME=(conservation pae pdb mock/conway mock/conway_genomic mock/plddt mock/secondary mock/structures mock/variants mock/adpred mock/paddle mock/ptms mock/slims mock/api mock/tfs downloads)
for d in "${RUNTIME[@]}"; do
  echo "== sync public/$d -> $DEST/$d =="
  rclone sync "public/$d" "$DEST/$d" $FLAGS
done

echo "== upload build_data.tar.gz -> $DEST/ =="
rclone copy /tmp/build_data.tar.gz "$DEST/" $FLAGS

# Single-file runtime assets the API reads (the /stats aggregate + the TF
# directory index for /tfs).
echo "== upload stats_constraint.json + index.json -> $DEST/mock/ =="
rclone copy public/mock/stats_constraint.json "$DEST/mock/" $FLAGS
rclone copy public/mock/index.json "$DEST/mock/" $FLAGS

echo "Done. Set the site's PUBLIC_DATA_BASE to:  <bucket public URL>/$VER"
