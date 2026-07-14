#!/usr/bin/env bash
# =============================================================================
# download_mouse_visium.sh
# =============================================================================
# Downloads the mouse brain Visium spatial transcriptomics data (E-MTAB-11114)
# from the EBI BioStudies FTP server and extracts the spatial archives.
#
# Source:
#   European Bioinformatics Institute BioStudies
#   Accession: E-MTAB-11114
#   Paper: Kleshchevnikov et al. 2022, Nature Biotechnology
#
# Total size: ~5.5GB (includes .cloupe files not used by pipeline)
#
# Sections downloaded (5 sections across 2 mice):
#   ST8059048  — used for reproduction (spatial mapping)
#   ST8059049  — available for additional analysis
#   ST8059050  — available for additional analysis
#   ST8059051  — available for additional analysis
#   ST8059052  — available for additional analysis
#
# Files USED by the pipeline:
#   ST8059048_filtered_feature_bc_matrix.h5
#   ST8059048_spatial/spatial/  (extracted from ST8059048_spatial.tar.gz)
#       tissue_hires_image.png
#       tissue_lowres_image.png
#       tissue_positions_list.csv
#       scalefactors_json.json
#       aligned_fiducials.jpg
#       detected_tissue_image.jpg
#
# IMPORTANT: spatial path has extra nesting: ST8059048_spatial/spatial/
# NOT ST8059048_spatial/ — script 04 points to the inner directory.
#
# Usage:
#   bash shell_scripts/download_mouse_visium.sh
#   bash shell_scripts/download_mouse_visium.sh --data-dir /custom/path/data
#   bash shell_scripts/download_mouse_visium.sh --section ST8059048  # one section only
# =============================================================================

set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/../data"
SECTION_FILTER=""   # empty = download all sections

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir) DATA_DIR="$2"; shift 2 ;;
        --section)  SECTION_FILTER="$2"; shift 2 ;;
        -h|--help)
            sed -n 's/^# \?//p' "$0" | head -40
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

DATA_DIR="$(realpath "${DATA_DIR}")"
OUT_DIR="${DATA_DIR}/mouseVisium"
FTP_BASE="ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/114/E-MTAB-11114/Files"

ALL_SECTIONS=(
    "ST8059048"
    "ST8059049"
    "ST8059050"
    "ST8059051"
    "ST8059052"
)

if [[ -n "$SECTION_FILTER" ]]; then
    SECTIONS=("$SECTION_FILTER")
    echo "Downloading single section: ${SECTION_FILTER}"
else
    SECTIONS=("${ALL_SECTIONS[@]}")
fi

# ── Setup ─────────────────────────────────────────────────────────────────────
echo "============================================================"
echo " DOWNLOADING MOUSE BRAIN VISIUM (E-MTAB-11114)"
echo "============================================================"
echo " Output directory : ${OUT_DIR}"
echo " Source           : ${FTP_BASE}/"
echo " Sections         : ${SECTIONS[*]}"
echo " Estimated size   : ~5.5GB (all sections) / ~1.1GB (ST8059048 only)"
echo ""

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# ── Download metadata files ───────────────────────────────────────────────────
echo "[1/3] Downloading metadata files..."

META_FILES=(
    "E-MTAB-11114.idf.txt"
    "E-MTAB-11114.sdrf.txt"
)

for MFILE in "${META_FILES[@]}"; do
    if [[ -f "${MFILE}" && -s "${MFILE}" ]]; then
        echo "  Skipping (exists): ${MFILE}"
    else
        echo "  Downloading: ${MFILE}"
        wget --quiet --show-progress --continue \
            "${FTP_BASE}/${MFILE}" -O "${MFILE}"
    fi
done

# ── Download section files ────────────────────────────────────────────────────
echo ""
echo "[2/3] Downloading section files..."

# Per-section file suffixes available on EBI
FILE_SUFFIXES=(
    "cloupe.cloupe"
    "filtered_feature_bc_matrix.h5"
    "metrics_summary.csv"
    "molecule_info.h5"
    "raw_feature_bc_matrix.h5"
    "spatial.tar.gz"
    "web_summary.html"
)

for SECTION in "${SECTIONS[@]}"; do
    echo ""
    echo "  Section: ${SECTION}"
    for SUFFIX in "${FILE_SUFFIXES[@]}"; do
        FNAME="${SECTION}_${SUFFIX}"

        # Skip large files not needed by the pipeline
        if [[ "${SUFFIX}" == "cloupe.cloupe" ]]; then
            echo "    Skipping (not needed): ${FNAME}  (~750MB)"
            continue
        fi

        if [[ -f "${FNAME}" && -s "${FNAME}" ]]; then
            SIZE=$(du -h "${FNAME}" | cut -f1)
            echo "    Skipping (exists): ${FNAME}  (${SIZE})"
        else
            echo "    Downloading: ${FNAME}"
            wget --quiet --show-progress --continue \
                "${FTP_BASE}/${FNAME}" -O "${FNAME}" || \
            echo "    WARNING: Failed to download ${FNAME}"
        fi
    done
done

# ── Extract all spatial archives ──────────────────────────────────────────────
echo ""
echo "[3/3] Extracting spatial archives..."
echo "      (Creates {section}_spatial/spatial/ for each section)"

for SECTION in "${SECTIONS[@]}"; do
    TAR="${SECTION}_spatial.tar.gz"
    SPATIAL_DIR="${SECTION}_spatial"

    if [[ ! -f "${TAR}" ]]; then
        echo "  MISSING: ${TAR} — skipping extraction"
        continue
    fi

    if [[ -d "${SPATIAL_DIR}/spatial" ]]; then
        echo "  Already extracted: ${SPATIAL_DIR}/spatial/"
        continue
    fi

    echo "  Extracting: ${TAR} -> ${SPATIAL_DIR}/"
    mkdir -p "${SPATIAL_DIR}"
    tar -xzf "${TAR}" -C "${SPATIAL_DIR}"

    # Verify expected files exist in the spatial subfolder
    SPATIAL_INNER="${SPATIAL_DIR}/spatial"
    if [[ -d "${SPATIAL_INNER}" ]]; then
        echo "    OK: ${SPATIAL_INNER}/"
        ls "${SPATIAL_INNER}/" | sed 's/^/      /'
    else
        # Some versions extract without the inner spatial/ subfolder
        # Move files into the expected location
        echo "    Restructuring to expected path: ${SPATIAL_INNER}/"
        mkdir -p "${SPATIAL_INNER}"
        mv "${SPATIAL_DIR}"/*.png "${SPATIAL_INNER}/" 2>/dev/null || true
        mv "${SPATIAL_DIR}"/*.jpg "${SPATIAL_INNER}/" 2>/dev/null || true
        mv "${SPATIAL_DIR}"/*.json "${SPATIAL_INNER}/" 2>/dev/null || true
        mv "${SPATIAL_DIR}"/*.csv "${SPATIAL_INNER}/" 2>/dev/null || true
        echo "    Files moved to: ${SPATIAL_INNER}/"
    fi
done

# ── Verify required files ─────────────────────────────────────────────────────
echo ""
echo "Verifying required pipeline files for ST8059048 (primary section)..."

REQUIRED=(
    "ST8059048_filtered_feature_bc_matrix.h5"
    "ST8059048_spatial/spatial/tissue_hires_image.png"
    "ST8059048_spatial/spatial/tissue_lowres_image.png"
    "ST8059048_spatial/spatial/tissue_positions_list.csv"
    "ST8059048_spatial/spatial/scalefactors_json.json"
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
    echo " Download and extraction complete."
    echo ""
    echo " Next step:"
    echo "   bash shell_scripts/download_gbm_visium.sh"
else
    echo " WARNING: ${MISSING} required file(s) missing."
    echo " Check your network connection and re-run."
    exit 1
fi
echo "============================================================"