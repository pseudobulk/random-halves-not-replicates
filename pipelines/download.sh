#!/usr/bin/env bash
# Fetch the public data listed by code/make_manifests.py into DATA_ROOT (config.env) and check sizes and checksums.
# Re-running skips files that are already complete. About 1.04 TB in 3,401 files: CD4 single cells 0.62 TB,
# Tahoe-100M shards 0.34 TB, Replogle 0.09 TB.
#   bash pipelines/download.sh             build the manifest, then download
#   bash pipelines/download.sh <manifest>  download an existing manifest
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CFG="${RHNR_CONFIG:-$REPO/config.env}"
[ -f "$CFG" ] || { echo "missing $CFG: copy config.example.env to config.env and edit it"; exit 1; }
_env="$(export -p)"; set -a; . "$CFG"; set +a; eval "$_env"   # variables already exported take precedence
M="${1:-$DATA_ROOT/manifest.tsv}"
mkdir -p "$DATA_ROOT"
[ -f "$M" ] || "${PYTHON:-python3}" "$REPO/code/make_manifests.py" "$M"
bad=0
while IFS=$'\t' read -r url dst size ctype csum; do
  f="$DATA_ROOT/$dst"; mkdir -p "$(dirname "$f")"
  if [ -f "$f" ] && [ "$(stat -c %s "$f")" = "$size" ]; then continue; fi
  echo "get $dst"
  curl -L --fail --retry 5 -C - -o "$f" "$url" || { echo "FAILED $dst"; bad=1; continue; }
  case "$ctype" in
    sha256) [ "$(sha256sum "$f" | cut -d' ' -f1)" = "$csum" ] || { echo "CHECKSUM $dst"; bad=1; } ;;
    md5)    [ "$(md5sum "$f" | cut -d' ' -f1)" = "$csum" ] || { echo "CHECKSUM $dst"; bad=1; } ;;
    *)      [ "$(stat -c %s "$f")" = "$size" ] || { echo "SIZE $dst"; bad=1; } ;;
  esac
done < "$M"
echo "Not in the manifest, fetch by hand (README, Data):"
echo "  $DATA_ROOT/nadig/jurkat.h5ad, $DATA_ROOT/nadig/hepg2.h5ad   (GEO GSE264667)"
echo "  $DATA_ROOT/hpa/rna_celline.tsv                               (Human Protein Atlas)"
echo "  $DATA_ROOT/vcc2025/adata_Training.h5ad                       (VCC 2025 H1 training data, Arc Virtual Cell Atlas)"
exit $bad
