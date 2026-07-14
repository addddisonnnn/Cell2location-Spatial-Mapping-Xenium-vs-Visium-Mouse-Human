#!/usr/bin/env bash
# =============================================================================
# download_gbm_data.sh
# =============================================================================
# Downloads two datasets needed for the human glioblastoma extension:
#
# 1. GBmap scRNA-seq reference (NCBI GEO GSE211376, Ruiz-Moreno et al. 2022)
#       data/GBmap_reference/GBmap_raw_counts.tsv.gz       (~300MB)
#       data/GBmap_reference/GBmap_metadata.csv.gz         (~5MB)
#
# 2. Human glioblastoma Visium spatial data (10x Genomics)
#       data/humanGlioblastoma/
#           Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5
#           spatial/  (extracted from spatial.tar.gz)
#               tissue_hires_image.png
#               tissue_lowres_image.png
#               tissue_positions_list.csv
#               scalefactors_json.json
#
# Total size: ~310MB
#
# After downloading, run:
#   python scripts/02_prepare_gbmap_reference.py
# to build the compressed h5ad from the raw GBmap files.
#
# Usage:
#   bash shell_scripts/download_gbm_data.sh
#   bash shell_scripts/download_gbm_data.sh --data-dir /custom/path/data
#   bash shell_scripts/download_gbm_data.sh --skip-reference   # Visium only
#   bash shell_scripts/download_gbm_data.sh --skip-visium      # reference only
# =============================================================================

set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/../data"
SKIP_REFERENCE=false
SKIP_VISIUM=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)       DATA_DIR="$2"; shift 2 ;;
        --skip-reference) SKIP_REFERENCE=true; shift ;;
        --skip-visium)    SKIP_VISIUM=true; shift ;;
        -h|--help)
            sed -n 's/^# \?//p' "$0" | head -40
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

DATA_DIR="$(realpath "${DATA_DIR}")"
REF_DIR="${DATA_DIR}/GBmap_reference"
GBM_DIR="${DATA_DIR}/humanGlioblastoma"

echo "============================================================"
echo " DOWNLOADING HUMAN GLIOBLASTOMA DATASETS"
echo "============================================================"
echo " Data directory : ${DATA_DIR}"
echo ""

# ── Part 1: GBmap scRNA-seq reference ─────────────────────────────────────────
if [[ "${SKIP_REFERENCE}" == "false" ]]; then
    echo "------------------------------------------------------------"
    echo " PART 1: GBmap reference (GSE211376)"
    echo " Source: NCBI GEO FTP"
    echo "------------------------------------------------------------"
    echo ""

    mkdir -p "${REF_DIR}"
    cd "${REF_DIR}"

    GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211376/suppl"

    COUNTS_URL="${GEO_BASE}/GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz"
    META_URL="${GEO_BASE}/GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz"

    # Counts matrix (~300MB)
    if [[ -f "GBmap_raw_counts.tsv.gz" && -s "GBmap_raw_counts.tsv.gz" ]]; then
        SIZE=$(du -h GBmap_raw_counts.tsv.gz | cut -f1)
        echo "  Skipping (exists): GBmap_raw_counts.tsv.gz  (${SIZE})"
    else
        echo "  Downloading: GBmap_raw_counts.tsv.gz  (~300MB, takes a few minutes)"
        wget --quiet --show-progress --continue \
            "${COUNTS_URL}" \
            -O "GBmap_raw_counts.tsv.gz"
        echo "  Done: $(du -h GBmap_raw_counts.tsv.gz | cut -f1)"
    fi

    # Metadata CSV (~5MB)
    if [[ -f "GBmap_metadata.csv.gz" && -s "GBmap_metadata.csv.gz" ]]; then
        SIZE=$(du -h GBmap_metadata.csv.gz | cut -f1)
        echo "  Skipping (exists): GBmap_metadata.csv.gz  (${SIZE})"
    else
        echo "  Downloading: GBmap_metadata.csv.gz  (~5MB)"
        wget --quiet --show-progress --continue \
            "${META_URL}" \
            -O "GBmap_metadata.csv.gz"
        echo "  Done: $(du -h GBmap_metadata.csv.gz | cut -f1)"
    fi

    # Quick integrity check — verify files are valid gzip
    echo ""
    echo "  Verifying file integrity..."
    for f in "GBmap_raw_counts.tsv.gz" "GBmap_metadata.csv.gz"; do
        if gzip -t "${f}" 2>/dev/null; then
            NLINES=$(zcat "${f}" | head -2 | wc -l)
            echo "    OK (valid gzip, ${NLINES} header lines readable): ${f}"
        else
            echo "    CORRUPTED: ${f} — re-run this script to re-download"
        fi
    done

    echo ""
    echo "  GBmap reference downloaded."
    echo "  Run: python scripts/02_prepare_gbmap_reference.py"
    echo "  to build GBmap_reference_compressed.h5ad"
    echo ""
else
    echo "  Skipping GBmap reference download (--skip-reference)"
    echo ""
fi

# ── Part 2: Human GBM Visium spatial data ────────────────────────────────────
if [[ "${SKIP_VISIUM}" == "false" ]]; then
    echo "------------------------------------------------------------"
    echo " PART 2: Human GBM Visium spatial data (10x Genomics)"
    echo " Source: 10x Genomics CDN"
    echo "------------------------------------------------------------"
    echo ""

    mkdir -p "${GBM_DIR}"
    cd "${GBM_DIR}"

    CDN_BASE="https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma"
    H5_FILE="Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
    TAR_FILE="Parent_Visium_Human_Glioblastoma_spatial.tar.gz"

    # Filtered count matrix (~18MB)
    if [[ -f "${H5_FILE}" && -s "${H5_FILE}" ]]; then
        SIZE=$(du -h "${H5_FILE}" | cut -f1)
        echo "  Skipping (exists): ${H5_FILE}  (${SIZE})"
    else
        echo "  Downloading: ${H5_FILE}  (~18MB)"
        wget --quiet --show-progress --continue \
            "${CDN_BASE}/${H5_FILE}" \
            -O "${H5_FILE}"
        echo "  Done: $(du -h "${H5_FILE}" | cut -f1)"
    fi

    # Spatial archive with H&E image and coordinates (~10MB)
    if [[ -f "${TAR_FILE}" && -s "${TAR_FILE}" ]]; then
        SIZE=$(du -h "${TAR_FILE}" | cut -f1)
        echo "  Skipping (exists): ${TAR_FILE}  (${SIZE})"
    else
        echo "  Downloading: ${TAR_FILE}  (~10MB)"
        wget --quiet --show-progress --continue \
            "${CDN_BASE}/${TAR_FILE}" \
            -O "${TAR_FILE}"
        echo "  Done: $(du -h "${TAR_FILE}" | cut -f1)"
    fi

    # Extract spatial archive
    if [[ -d "spatial" && -f "spatial/tissue_positions_list.csv" ]]; then
        echo "  Already extracted: spatial/"
    else
        echo "  Extracting: ${TAR_FILE} -> spatial/"
        tar -xzf "${TAR_FILE}"
        echo "  Extracted:"
        ls spatial/ | sed 's/^/    /'
    fi

    # Verify
    echo ""
    echo "  Verifying required files..."
    REQUIRED=(
        "${H5_FILE}"
        "spatial/tissue_hires_image.png"
        "spatial/tissue_lowres_image.png"
        "spatial/tissue_positions_list.csv"
        "spatial/scalefactors_json.json"
        "spatial/aligned_fiducials.jpg"
    )

    MISSING=0
    for f in "${REQUIRED[@]}"; do
        if [[ -f "$f" && -s "$f" ]]; then
            SIZE=$(du -h "$f" | cut -f1)
            printf "    OK      %-60s %s\n" "$f" "(${SIZE})"
        else
            printf "    MISSING %-60s\n" "$f"
            MISSING=$((MISSING + 1))
        fi
    done

    echo ""
    echo "  Total directory size: $(du -sh "${GBM_DIR}" | cut -f1)"
    echo ""

    if [[ $MISSING -gt 0 ]]; then
        echo "  WARNING: ${MISSING} required file(s) missing."
        exit 1
    fi

    echo "  Human GBM Visium data ready."
    echo ""
else
    echo "  Skipping GBM Visium download (--skip-visium)"
    echo ""
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo "============================================================"
echo " Download complete."
echo ""
echo " Next steps:"
echo "   1. python scripts/02_prepare_gbmap_reference.py"
echo "      (builds GBmap_reference_compressed.h5ad)"
echo ""
echo "   2. python scripts/05_gbm_spatial_mapping.py \\"
echo "          --n-iter-reg 100 --n-iter-spat 500 --run-name test"
echo "      (quick test run)"
echo ""
echo "   3. python scripts/06_gbm_parameter_sweep.py"
echo "      (full 5-configuration sweep)"
echo "============================================================"