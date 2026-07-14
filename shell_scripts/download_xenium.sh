#!/usr/bin/env bash
# =============================================================================
# download_xenium.sh
# =============================================================================
# Downloads the 10x Xenium mouse brain coronal subset dataset and
# saves it to data/xeniumMouseBrain/.
#
# Source:
#   10x Genomics datasets portal
#   Dataset: Xenium V1 FF Mouse Brain Coronal Subset (CTX + HP)
#   Version: 1.0.2
#   URL: https://www.10xgenomics.com/datasets/xenium-ff-mouse-brain-coronal-subset-ctx-hp-1-0-2
#
# Dataset specifications:
#   Cells:               36,602
#   Genes (panel):       248
#   Mean transcripts:    246.7 per cell
#   Region:              Cerebral cortex (CTX) + hippocampus (HP)
#   Resolution:          ~10 µm per location (near single-cell)
#
# Files downloaded (~zip archive, extracts to ~2GB):
#   cell_feature_matrix.h5        — cell x gene count matrix (USED)
#   cells.csv.gz                  — spatial coordinates, cell metrics (USED)
#   morphology.ome.tif            — fluorescence morphology image (USED)
#   morphology_focus.ome.tif      — focused morphology image
#   morphology_mip.ome.tif        — max intensity projection
#   gene_panel.json               — 275 targeted genes (USED)
#   experiment.xenium             — experiment metadata
#   metrics_summary.csv           — summary QC metrics
#   analysis_summary.html         — QC report
#   cell_boundaries.csv.gz        — cell polygon boundaries
#   cell_boundaries.parquet       — cell polygon boundaries (parquet)
#   nucleus_boundaries.csv.gz     — nucleus polygon boundaries
#   nucleus_boundaries.parquet    — nucleus boundaries (parquet)
#   transcripts.csv.gz            — individual transcript locations
#   transcripts.parquet           — individual transcripts (parquet)
#   analysis/                     — precomputed clustering
#
# Usage:
#   bash shell_scripts/download_xenium.sh
#   bash shell_scripts/download_xenium.sh --data-dir /custom/path/data
# =============================================================================

set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/../data"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir) DATA_DIR="$2"; shift 2 ;;
        -h|--help)
            sed -n 's/^# \?//p' "$0" | head -40
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

DATA_DIR="$(realpath "${DATA_DIR}")"
OUT_DIR="${DATA_DIR}/xeniumMouseBrain"

ZIP_URL="https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip"
ZIP_FILE="Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip"

echo "============================================================"
echo " DOWNLOADING XENIUM MOUSE BRAIN DATA"
echo "============================================================"
echo " Output directory : ${OUT_DIR}"
echo " Dataset          : Xenium V1 FF Mouse Brain Coronal Subset (CTX+HP)"
echo " Version          : 1.0.2"
echo " Cells            : 36,602"
echo " Genes            : 248"
echo ""

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# ── Download zip archive ──────────────────────────────────────────────────────
echo "[1/3] Downloading zip archive..."
echo ""

if [[ -f "${ZIP_FILE}" && -s "${ZIP_FILE}" ]]; then
    SIZE=$(du -h "${ZIP_FILE}" | cut -f1)
    echo "  Skipping (exists): ${ZIP_FILE}  (${SIZE})"
else
    echo "  URL: ${ZIP_URL}"
    echo "  Downloading... (this may take several minutes)"
    wget --quiet --show-progress --continue \
        "${ZIP_URL}" \
        -O "${ZIP_FILE}"
    echo "  Done: $(du -h "${ZIP_FILE}" | cut -f1)"
fi

# ── Extract zip archive ───────────────────────────────────────────────────────
echo ""
echo "[2/3] Extracting archive..."
echo ""

if [[ -f "cell_feature_matrix.h5" && -f "cells.csv.gz" ]]; then
    echo "  Already extracted (cell_feature_matrix.h5 and cells.csv.gz found)"
else
    echo "  Extracting: ${ZIP_FILE}"
    unzip -q "${ZIP_FILE}" -d "${OUT_DIR}" 2>/dev/null || \
    unzip "${ZIP_FILE}" -d "${OUT_DIR}"

    # Handle case where zip extracts into a subdirectory
    # (some versions of the dataset extract into an 'outs/' folder)
    if [[ -d "outs" && -f "outs/cell_feature_matrix.h5" ]]; then
        echo "  Files extracted into outs/ subfolder — moving to current directory..."
        mv outs/* . 2>/dev/null || true
        rmdir outs 2>/dev/null || true
    fi

    echo "  Extraction complete."
fi

echo ""
echo "  Directory contents:"
ls -lh "${OUT_DIR}" | grep -v "^total" | sed 's/^/  /'

# ── Verify required files ─────────────────────────────────────────────────────
echo ""
echo "[3/3] Verifying required pipeline files..."

REQUIRED=(
    "cell_feature_matrix.h5"
    "cells.csv.gz"
    "gene_panel.json"
    "morphology.ome.tif"
)

MISSING=0
for f in "${REQUIRED[@]}"; do
    if [[ -f "$f" && -s "$f" ]]; then
        SIZE=$(du -h "$f" | cut -f1)
        printf "  OK      %-45s %s\n" "$f" "(${SIZE})"
    else
        printf "  MISSING %-45s\n" "$f"
        MISSING=$((MISSING + 1))
    fi
done

# Check gene panel count
if [[ -f "gene_panel.json" ]]; then
    NGENES=$(python3 -c "
import json
with open('gene_panel.json') as f:
    panel = json.load(f)
print(len(panel['payload']['targets']))
" 2>/dev/null || echo "unknown")
    echo ""
    echo "  Genes in panel: ${NGENES}  (expected 275)"
fi

echo ""
echo "Total directory size: $(du -sh "${OUT_DIR}" | cut -f1)"

echo ""
echo "============================================================"
if [[ $MISSING -eq 0 ]]; then
    echo " Download and extraction complete."
    echo ""
    echo " Next step:"
    echo "   python scripts/07_xenium_spatial_mapping.py \\"
    echo "       --n-iter-reg 100 --n-iter-spat 500 --n-comb 20 \\"
    echo "       --run-name run1"
    echo ""
    echo " Note: The Xenium reference atlas (~700MB) is downloaded"
    echo " automatically by script 07 if not already present."
else
    echo " WARNING: ${MISSING} required file(s) missing."
    echo " The download may have been interrupted. Re-run this script."
    exit 1
fi
echo "============================================================"