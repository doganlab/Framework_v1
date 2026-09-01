#!/usr/bin/env bash
set -euo pipefail

# Framework_v1 - run Ensembl VEP on current mGWAS-derived variants (GRCh38)

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="${PROJECT_ROOT}/data/intermediate/variant_annotation/vep_coordinates_GRCh38.txt"
OUTPUT="${PROJECT_ROOT}/data/intermediate/variant_annotation/VEP_annotation_GRCh38.tsv"

VEP_BIN="${VEP_BIN:-vep}"
VEP_CACHE_DIR="${VEP_CACHE_DIR:-/home/dplab3/utilities/vep_cache}"
VEP_VERSION="${VEP_VERSION:-116}"
CACHE_DIR="${VEP_CACHE_DIR}/homo_sapiens/${VEP_VERSION}_GRCh38"

if [[ ! -s "${INPUT}" ]]; then
  echo "ERROR: VEP input is missing or empty: ${INPUT}" >&2
  exit 1
fi

if ! command -v "${VEP_BIN}" >/dev/null 2>&1; then
  echo "ERROR: VEP executable not found: ${VEP_BIN}" >&2
  exit 1
fi

if [[ ! -d "${CACHE_DIR}" ]]; then
  echo "ERROR: matching GRCh38 VEP cache not found:" >&2
  echo "  ${CACHE_DIR}" >&2
  echo "Install/check the homo_sapiens ${VEP_VERSION}_GRCh38 cache before continuing." >&2
  exit 1
fi

echo "VEP executable: $(command -v "${VEP_BIN}")"
"${VEP_BIN}" --version || true
echo "Input:  ${INPUT}"
echo "Output: ${OUTPUT}"
echo "Cache:  ${CACHE_DIR}"

"${VEP_BIN}" \
  --assembly GRCh38 \
  --species homo_sapiens \
  --cache \
  --offline \
  --dir_cache "${VEP_CACHE_DIR}" \
  --format ensembl \
  --input_file "${INPUT}" \
  --output_file "${OUTPUT}" \
  --tab \
  --symbol \
  --biotype \
  --canonical \
  --variant_class \
  --nearest symbol \
  --distance 5000 \
  --fork 8 \
  --force_overwrite

echo
echo "VEP completed."
echo "Output: ${OUTPUT}"
echo "Summary HTML should be: ${OUTPUT}_summary.html"
