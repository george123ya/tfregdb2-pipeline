#!/usr/bin/env bash
# Production build for Cloudflare Pages.
#
# Builds the site, then removes the data dirs that are served from Cloudflare
# R2 (via PUBLIC_DATA_BASE) so the Pages deploy stays under the 20,000-file cap
# and well under size limits. Build-only data (mock/tfs, slims, adpred, index)
# is already inlined into the generated HTML at build time, so it's dropped too.
#
# PUBLIC_DATA_BASE must be set at BUILD time (Astro inlines it into the client
# islands). On Cloudflare Pages set it as an env var; locally:
#   PUBLIC_DATA_BASE=https://<bucket>.r2.dev bash scripts/build_pages.sh
#
# The same data dirs must be uploaded to R2 (scripts/upload_r2.sh, added once
# the bucket exists). Local `astro dev` needs none of this — with
# PUBLIC_DATA_BASE unset, dataUrl() stays relative and serves from public/.
set -euo pipefail
cd "$(dirname "$0")/.."

# The R2 public base is a non-secret constant (also the default baked into the
# API Functions). Default it here so the build never depends on a dashboard env
# var being wired correctly — a CLI/dashboard-provided value still overrides it.
# Astro inlines PUBLIC_DATA_BASE into client islands at build, so export it.
: "${PUBLIC_DATA_BASE:=https://pub-67d32012fe4c4df48658f879a21c57b8.r2.dev/v1}"
export PUBLIC_DATA_BASE
echo "== PUBLIC_DATA_BASE=${PUBLIC_DATA_BASE} =="

# CI builds clone a SOURCE-ONLY repo (no data). The build inlines per-TF data at
# build time, so pull the build-input tarball from R2 first. Local builds (data
# already in public/) skip this. Built + uploaded by scripts/upload_r2.sh.
if [ ! -d public/mock/tfs ]; then
  : "${PUBLIC_DATA_BASE:?build data missing and PUBLIC_DATA_BASE unset}"
  echo "== fetching build_data.tar.gz from ${PUBLIC_DATA_BASE} =="
  curl -fsSL "${PUBLIC_DATA_BASE%/}/build_data.tar.gz" -o /tmp/build_data.tar.gz
  tar xzf /tmp/build_data.tar.gz   # extracts public/mock/{tfs, slims, secondary, ...}
  echo "   extracted: $(find public/mock/tfs -type f | wc -l) TF files"
fi

npm run build

echo "== stripping off-host data from dist/ =="
# Served from R2 at runtime:
rm -rf dist/pae dist/pdb dist/conservation
rm -rf dist/mock/conway dist/mock/plddt dist/mock/secondary dist/mock/structures dist/mock/variants
# Build-only (already inlined into HTML — never fetched at runtime):
rm -rf dist/mock/tfs dist/mock/slims dist/mock/adpred dist/mock/index dist/mock/index.json
# stray data archives that shouldn't ship to Pages:
rm -f dist/mock/*.tar.gz dist/mock/*.zip

echo "== Pages deploy footprint =="
printf "files: %s\n" "$(find dist -type f | wc -l)"
du -sh dist
echo "top-level:"; ls dist
