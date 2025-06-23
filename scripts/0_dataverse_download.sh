#!/bin/bash
#SBATCH --job-name=dataverse_pull
#SBATCH --output=../log/dataverse_pull_%j.out
#SBATCH --error=../log/dataverse_pull_%j.err
#SBATCH --time=01-00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=sapphire


set -euo pipefail

# ── 1. Move to repo root ─────────────────────────────────────────────────────
cd "$(dirname "$0")/.."

# ── 2. Config ────────────────────────────────────────────────────────────────
DOI="doi:10.7910/DVN/NW0ZN1"
ZIP_NAME="dataset_bundle.zip"
TMP_DIR="tmp_extracted"

# ── 3. Ensure required folders exist ─────────────────────────────────────────
mkdir -p data/2C data/4C ref output/seurat log

# ── 4. Download from Dataverse ───────────────────────────────────────────────
echo "📥 Downloading dataset from Dataverse..."
curl -L "https://dataverse.harvard.edu/api/access/dataset/:persistentId/?persistentId=${DOI}" -o "$ZIP_NAME"

# ── 5. Extract ───────────────────────────────────────────────────────────────
echo "📦 Extracting archive..."
UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE unzip "$ZIP_NAME" -d "$TMP_DIR"
rm "$ZIP_NAME"

# ── 6. Move files ────────────────────────────────────────────────────────────
# Move data folders
if [ -d "${TMP_DIR}/data/2C" ]; then mv "${TMP_DIR}/data/2C" data/; fi
if [ -d "${TMP_DIR}/data/4C" ]; then mv "${TMP_DIR}/data/4C" data/; fi

# Move ref
if [ -d "${TMP_DIR}/ref" ]; then mv "${TMP_DIR}/ref/"* ref/; fi

# Move seurat
if [ -d "${TMP_DIR}/seurat" ]; then mv "${TMP_DIR}/seurat/"* output/seurat/; fi

# ── 7. Clean up ──────────────────────────────────────────────────────────────
rm -rf "$TMP_DIR"

echo "✅ Dataverse dataset downloaded and unpacked into repo folders."