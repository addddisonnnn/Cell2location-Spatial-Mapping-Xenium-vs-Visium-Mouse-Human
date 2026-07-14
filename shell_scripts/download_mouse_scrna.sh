#!/usr/bin/env bash
# =============================================================================
# download_mouse_scrna.sh
# =============================================================================
# Downloads the mouse brain snRNA-seq reference data (E-MTAB-11115) from the
# EBI BioStudies FTP server and saves it to data/mousescRNAseq/.
#
# Source:
#   European Bioinformatics Institute BioStudies
#   Accession: E-MTAB-11115
#   Paper: Kleshchevnikov et al. 2022, Nature Biotechnology
#
# Total size: ~2.8GB
#
# Files USED by the pipeline:
#   5705STDY8058280_filtered_feature_bc_matrix.h5  (x6 samples)
#   cell_annotation.csv
#
# Usage:
#   bash shell_scripts/download_mouse_scrna.sh
#   bash shell_scripts/download_mouse_scrna.sh --data-dir /custom/path/data
# =============================================================================

set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/../data"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir) DATA_DIR="$2"; shift 2 ;;
        -h|--help)
            sed -n 's/^# \?//p' "$0" | head -30
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

DATA_DIR="$(realpath "${DATA_DIR}")"
OUT_DIR="${DATA_DIR}/mousescRNAseq"
FTP_BASE="ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/115/E-MTAB-11115/Files"

# ── Setup ─────────────────────────────────────────────────────────────────────
echo "============================================================"
echo " DOWNLOADING MOUSE BRAIN snRNA-seq (E-MTAB-11115)"
echo "============================================================"
echo " Output directory : ${OUT_DIR}"
echo " Source           : ${FTP_BASE}/"
echo " Estimated size   : ~2.8GB"
echo ""

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# ── Download individual files ─────────────────────────────────────────────────
# Downloading individual files rather than recursive wget gives better
# progress reporting and allows skipping files already present.

SAMPLES=(
    "5705STDY8058280"
    "5705STDY8058281"
    "5705STDY8058282"
    "5705STDY8058283"
    "5705STDY8058284"
    "5705STDY8058285"
)

FILE_TYPES=(
    "filtered_feature_bc_matrix.h5"
    "metrics_summary.csv"
    "molecule_info.h5"
    "raw_feature_bc_matrix.h5"
    "web_summary.html"
)

echo "[1/3] Downloading sample files..."
for SAMPLE in "${SAMPLES[@]}"; do
    for FTYPE in "${FILE_TYPES[@]}"; do
        FNAME="${SAMPLE}_${FTYPE}"
        if [[ -f "${FNAME}" && -s "${FNAME}" ]]; then
            echo "  Skipping (exists): ${FNAME}"
        else
            echo "  Downloading: ${FNAME}"
            wget --quiet --show-progress --continue \
                "${FTP_BASE}/${FNAME}" -O "${FNAME}" || \
            echo "  WARNING: Failed to download ${FNAME}"
        fi
    done
done

echo ""
echo "[2/3] Downloading annotation and metadata files..."

META_FILES=(
    "cell_annotation.csv"
    "E-MTAB-11115.idf.txt"
    "E-MTAB-11115.sdrf.txt"
)

for MFILE in "${META_FILES[@]}"; do
    if [[ -f "${MFILE}" && -s "${MFILE}" ]]; then
        echo "  Skipping (exists): ${MFILE}"
    else
        echo "  Downloading: ${MFILE}"
        wget --quiet --show-progress --continue \
            "${FTP_BASE}/${MFILE}" -O "${MFILE}" || \
        echo "  WARNING: Failed to download ${MFILE}"
    fi
done

# ── Verify required files ─────────────────────────────────────────────────────
echo ""
echo "[3/3] Verifying required pipeline files..."

REQUIRED=(
    "5705STDY8058280_filtered_feature_bc_matrix.h5"
    "5705STDY8058281_filtered_feature_bc_matrix.h5"
    "5705STDY8058282_filtered_feature_bc_matrix.h5"
    "5705STDY8058283_filtered_feature_bc_matrix.h5"
    "5705STDY8058284_filtered_feature_bc_matrix.h5"
    "5705STDY8058285_filtered_feature_bc_matrix.h5"
    "cell_annotation.csv"
)

MISSING=0
for f in "${REQUIRED[@]}"; do
    if [[ -f "$f" && -s "$f" ]]; then
        SIZE=$(du -h "$f" | cut -f1)
        printf "  OK      %-65s %s\n" "$f" "(${SIZE})"
    else
        printf "  MISSING %-65s\n" "$f"
        MISSING=$((MISSING + 1))
    fi
done

echo ""
echo "Total directory size: $(du -sh "${OUT_DIR}" | cut -f1)"

echo ""
echo "============================================================"
if [[ $MISSING -eq 0 ]]; then
    echo " Download complete. All required files present."
    echo ""
    echo " Next step:"
    echo "   bash shell_scripts/download_mouse_visium.sh"
else
    echo " WARNING: ${MISSING} required file(s) missing."
    echo " Check your network connection and re-run this script."
    exit 1
fi
echo "============================================================"