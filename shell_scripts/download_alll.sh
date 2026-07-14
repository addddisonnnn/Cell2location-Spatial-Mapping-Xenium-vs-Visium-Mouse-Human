#!/usr/bin/env bash
# =============================================================================
# download_all.sh
# =============================================================================
# Master download script. Runs all individual download scripts in the correct
# order to prepare all datasets needed for the full pipeline.
#
# Downloads:
#   1. Mouse brain snRNA-seq (E-MTAB-11115)    ~2.8GB
#   2. Mouse brain Visium (E-MTAB-11114)       ~5.5GB
#   3. Human GBM Visium + GBmap reference      ~0.31GB
#   4. Xenium mouse brain coronal subset       ~2GB
#
# Total: ~10.6GB
#
# Usage:
#   bash shell_scripts/download_all.sh
#   bash shell_scripts/download_all.sh --data-dir /custom/path/data
#   bash shell_scripts/download_all.sh --skip-xenium     # skip Xenium download
#   bash shell_scripts/download_all.sh --only-reproduction  # mouse data only
# =============================================================================

set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/../data"
SKIP_XENIUM=false
ONLY_REPRODUCTION=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)         DATA_DIR="$2"; shift 2 ;;
        --skip-xenium)      SKIP_XENIUM=true; shift ;;
        --only-reproduction) ONLY_REPRODUCTION=true; shift ;;
        -h|--help)
            sed -n 's/^# \?//p' "$0" | head -35
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

DATA_DIR="$(realpath "${DATA_DIR}")"

echo "============================================================"
echo " CELL2LOCATION — MASTER DOWNLOAD SCRIPT"
echo "============================================================"
echo " Data directory: ${DATA_DIR}"
echo ""
echo " Datasets to download:"
echo "   [1] Mouse brain snRNA-seq (E-MTAB-11115)    ~2.8GB"
echo "   [2] Mouse brain Visium (E-MTAB-11114)       ~5.5GB"
if [[ "${ONLY_REPRODUCTION}" == "false" ]]; then
    echo "   [3] Human GBM Visium + GBmap reference      ~0.31GB"
    if [[ "${SKIP_XENIUM}" == "false" ]]; then
        echo "   [4] Xenium mouse brain                       ~2GB"
        echo ""
        echo " Total estimated: ~10.6GB"
    else
        echo ""
        echo " Total estimated: ~8.6GB"
    fi
else
    echo ""
    echo " Total estimated: ~8.3GB"
fi
echo ""
echo "============================================================"
echo ""

read -r -p "Continue? [y/N] " CONFIRM
if [[ ! "${CONFIRM}" =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 0
fi

# ── Track timing ──────────────────────────────────────────────────────────────
START_TIME=$(date +%s)

run_script() {
    local script="$1"
    shift
    echo ""
    echo "============================================================"
    echo " Running: ${script}"
    echo "============================================================"
    bash "${SCRIPT_DIR}/${script}" --data-dir "${DATA_DIR}" "$@"
    local EXIT_CODE=$?
    if [[ $EXIT_CODE -ne 0 ]]; then
        echo "ERROR: ${script} failed with exit code ${EXIT_CODE}"
        exit $EXIT_CODE
    fi
}

# ── Run downloads ─────────────────────────────────────────────────────────────
run_script "download_mouse_scrna.sh"
run_script "download_mouse_visium.sh"

if [[ "${ONLY_REPRODUCTION}" == "false" ]]; then
    run_script "download_gbm_data.sh"

    if [[ "${SKIP_XENIUM}" == "false" ]]; then
        run_script "download_xenium.sh"
    else
        echo ""
        echo "Skipping Xenium download (--skip-xenium)"
    fi
fi

# ── Final summary ─────────────────────────────────────────────────────────────
END_TIME=$(date +%s)
ELAPSED=$(( (END_TIME - START_TIME) / 60 ))

echo ""
echo "============================================================"
echo " ALL DOWNLOADS COMPLETE"
echo "============================================================"
echo ""
echo " Time elapsed: ${ELAPSED} minutes"
echo ""
echo " Directory summary:"
du -sh "${DATA_DIR}"/*/  2>/dev/null | sort -h | sed 's/^/   /'
echo ""
echo " Total data directory size:"
du -sh "${DATA_DIR}" | sed 's/^/   /'
echo ""
echo "============================================================"
echo " NEXT STEPS"
echo "============================================================"
echo ""
echo " 1. Build the GBmap reference h5ad:"
echo "      python scripts/02_prepare_gbmap_reference.py"
echo ""
echo " 2. Start the Singularity container:"
echo "      singularity exec --nv data/cell2location.sif bash"
echo "      jupyter lab --ip=0.0.0.0 --port=8888 --no-browser"
echo ""
echo " 3. Run the reproduction pipeline:"
echo "      python scripts/03_reproduction_scrna_umap.py"
echo "      python scripts/04_reproduction_visium_mapping.py --n-iter-spatial 30000"
echo ""
echo " 4. Run the GBM extension:"
echo "      python scripts/06_gbm_parameter_sweep.py"
echo ""
echo " 5. Run the Xenium extension:"
echo "      python scripts/07_xenium_spatial_mapping.py --sweep"
echo "============================================================"
