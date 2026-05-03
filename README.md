# Cell2Location: Spatial Mapping of Cell Types in Mouse Brain, Human Glioblastoma, and Xenium

Reproduction and extension of Kleshchevnikov et al. (2022), *Nature Biotechnology* —
"Cell2location maps fine-grained cell types in spatial transcriptomics."

This repository reproduces Figures 2b–2e from the original paper on mouse brain Visium data,
then extends the method to two new contexts: human glioblastoma Visium data and 10x Xenium
data at near single-cell resolution.

---

## Table of Contents

- [Background](#background)
- [Project Summary](#project-summary)
- [Complete Repository Structure](#complete-repository-structure)
- [Prerequisites](#prerequisites)
- [Environment Setup](#environment-setup)
- [Downloading the Datasets](#downloading-the-datasets)
  - [Quick Start — Download Everything](#quick-start--download-everything)
  - [Shell Scripts Reference](#shell-scripts-reference)
  - [Manual Download Instructions](#manual-download-instructions)
- [Pipeline Overview](#pipeline-overview)
- [Scripts Reference](#scripts-reference)
- [Running the Full Pipeline](#running-the-full-pipeline)
- [Parameter Reference](#parameter-reference)
- [Expected Outputs](#expected-outputs)
- [Known Issues and Limitations](#known-issues-and-limitations)
- [References](#references)

---

## Background

Spatial transcriptomics technologies like 10x Genomics Visium capture gene expression
from intact tissue sections. Each measurement spot contains multiple cells mixed together,
making it impossible to directly identify which cell types are present. Cell2location
solves this by using annotated single-cell RNA-seq data as a reference and computationally
deconvolving the spatial signal into per-spot cell type abundances.

The model is Bayesian and hierarchical — it borrows statistical strength across spatially
neighboring locations, which is its key advantage over simpler methods for detecting rare
and transcriptionally similar cell types.

---

## Project Summary

| Component | Description |
|---|---|
| **Reproduction** | Figures 2b–2e from the original paper using mouse brain Visium + snRNA-seq |
| **Extension 1** | Human glioblastoma Visium with GBmap reference (14 cell types, 5-config parameter sweep) |
| **Extension 2** | Xenium mouse brain at near single-cell resolution (5 parameter configurations) |
| **Tool version** | cell2location v0.3 — PyMC3/Theano backend (newer v0.1.3 exceeded storage limits) |
| **Hardware** | NVIDIA V100 GPU required — Theano incompatible with A100/L40s |
| **Environment** | Singularity container, desktop interactive session |

---

## Complete Repository Structure

```
Cell2Location-Project-Research-Reproduction-Extension/
│
├── shell_scripts/                        ← Dataset download scripts (run first)
│   ├── download_all.sh                   ← Master script — runs all downloads
│   ├── download_mouse_scrna.sh           ← Mouse brain snRNA-seq (E-MTAB-11115)
│   ├── download_mouse_visium.sh          ← Mouse brain Visium (E-MTAB-11114)
│   ├── download_gbm_data.sh             ← GBmap reference + GBM Visium
│   └── download_xenium.sh               ← Xenium mouse brain coronal subset
│
├── scripts/                             ← Analysis pipeline (Python)
│   ├── 00_gpu_setup.py                  ← GPU config — import first in all scripts
│   ├── 01_download_gbmap_reference.py   ← Download GBmap from NCBI GEO (Python alt)
│   ├── 02_prepare_gbmap_reference.py    ← Build GBmap h5ad from raw files
│   ├── 03_reproduction_scrna_umap.py    ← Load mouse snRNA-seq, produce Figure 2b
│   ├── 04_reproduction_visium_mapping.py ← Figures 2c–2e (regression + spatial model)
│   ├── 05_gbm_spatial_mapping.py        ← GBM extension, single configuration
│   ├── 06_gbm_parameter_sweep.py        ← GBM extension, all 5 configurations
│   ├── 07_xenium_spatial_mapping.py     ← Xenium extension, 1 or 5 configurations
│   └── utils.py                         ← Shared helpers used across all scripts
│
├── data/                                ← Raw input data (not tracked by git)
│   ├── cell2location.sif                ← Singularity container (~4.8GB)
│   │
│   ├── mousescRNAseq/                   ← E-MTAB-11115 (~2.8GB)
│   │   ├── 5705STDY8058280_filtered_feature_bc_matrix.h5   (27MB)
│   │   ├── 5705STDY8058281_filtered_feature_bc_matrix.h5   (28MB)
│   │   ├── 5705STDY8058282_filtered_feature_bc_matrix.h5   (22MB)
│   │   ├── 5705STDY8058283_filtered_feature_bc_matrix.h5   (22MB)
│   │   ├── 5705STDY8058284_filtered_feature_bc_matrix.h5   (15MB)
│   │   ├── 5705STDY8058285_filtered_feature_bc_matrix.h5   (34MB)
│   │   ├── cell_annotation.csv                             (2.8MB) ← key file
│   │   ├── E-MTAB-11115.idf.txt
│   │   └── E-MTAB-11115.sdrf.txt
│   │
│   ├── mouseVisium/                     ← E-MTAB-11114 (~5.5GB)
│   │   ├── ST8059048_filtered_feature_bc_matrix.h5  (18MB) ← used for mapping
│   │   ├── ST8059048_spatial.tar.gz                 (9.9MB)
│   │   ├── ST8059048_spatial/                       ← extracted from tar
│   │   │   └── spatial/                             ← INNER folder (important!)
│   │   │       ├── tissue_hires_image.png
│   │   │       ├── tissue_lowres_image.png
│   │   │       ├── tissue_positions_list.csv
│   │   │       ├── scalefactors_json.json
│   │   │       ├── aligned_fiducials.jpg
│   │   │       └── detected_tissue_image.jpg
│   │   ├── ST8059049_filtered_feature_bc_matrix.h5
│   │   ├── ST8059049_spatial/spatial/
│   │   ├── ST8059050_filtered_feature_bc_matrix.h5
│   │   ├── ST8059050_spatial/spatial/
│   │   ├── ST8059051_filtered_feature_bc_matrix.h5
│   │   ├── ST8059051_spatial/spatial/
│   │   ├── ST8059052_filtered_feature_bc_matrix.h5
│   │   ├── ST8059052_spatial/spatial/
│   │   ├── E-MTAB-11114.idf.txt
│   │   └── E-MTAB-11114.sdrf.txt
│   │
│   ├── GBmap_reference/                 ← GSE211376 (~305MB raw + h5ad)
│   │   ├── GBmap_raw_counts.tsv.gz      ← downloaded (~300MB, 27k genes x 39k cells)
│   │   ├── GBmap_metadata.csv.gz        ← downloaded (~5MB, cell type labels)
│   │   ├── GBmap_reference_compressed.h5ad  ← built by script 02 (fast reload)
│   │   └── reference_cell_type_distribution.png
│   │
│   ├── humanGlioblastoma/               ← 10x Genomics GBM Visium (~28MB)
│   │   ├── Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5 (18MB)
│   │   ├── Parent_Visium_Human_Glioblastoma_spatial.tar.gz (10MB)
│   │   └── spatial/                     ← extracted from tar
│   │       ├── tissue_hires_image.png
│   │       ├── tissue_lowres_image.png
│   │       ├── tissue_positions_list.csv
│   │       ├── scalefactors_json.json
│   │       ├── aligned_fiducials.jpg
│   │       └── detected_tissue_image.jpg
│   │
│   └── xeniumMouseBrain/                ← 10x Xenium coronal subset (~2GB)
│       ├── cell_feature_matrix.h5       ← 36,602 cells x 248 genes (USED)
│       ├── cells.csv.gz                 ← x/y centroids + cell metrics (USED)
│       ├── morphology.ome.tif           ← fluorescence tissue image (USED)
│       ├── morphology_focus.ome.tif
│       ├── morphology_mip.ome.tif
│       ├── gene_panel.json              ← 275 targeted genes (USED)
│       ├── experiment.xenium
│       ├── metrics_summary.csv
│       ├── analysis_summary.html
│       ├── cell_boundaries.csv.gz
│       ├── cell_boundaries.parquet
│       ├── nucleus_boundaries.csv.gz
│       ├── nucleus_boundaries.parquet
│       ├── transcripts.csv.gz
│       ├── transcripts.parquet
│       └── analysis/
│
├── results/                             ← All outputs (not tracked by git)
│   ├── reproduction/
│   │   ├── adata_snrna_processed.h5ad   ← processed snRNA-seq (script 03 → 04)
│   │   ├── figure_2b_umap.png           ← Figure 2b reproduction
│   │   ├── inf_aver.csv                 ← reference signature matrix
│   │   ├── adata_vis_mapped.h5ad        ← Visium with cell type abundances
│   │   ├── figure_2c_HE.png             ← Figure 2c reproduction
│   │   ├── figure_2d_major.png          ← Figure 2d reproduction
│   │   ├── figure_2e_sparse.png         ← Figure 2e reproduction
│   │   ├── regression_output/
│   │   └── spatial_output/
│   │
│   ├── human_gbm/
│   │   ├── single_run/                  ← script 05 output (one config)
│   │   │   ├── inf_aver.csv
│   │   │   ├── adata_gbm_mapped.h5ad
│   │   │   ├── cell_abundance_summary.csv
│   │   │   ├── cell_type_composition.png
│   │   │   ├── dominant_cell_type_map.png
│   │   │   ├── six_key_cell_types.png
│   │   │   └── all_cell_types_grid.png
│   │   └── parameter_runs/              ← script 06 output (5 configs)
│   │       ├── run1_test/
│   │       │   ├── inf_aver.csv
│   │       │   ├── adata_gbm_mapped.h5ad
│   │       │   ├── cell_abundance_summary.csv
│   │       │   └── plots/
│   │       │       ├── composition.png
│   │       │       ├── dominant_cell_type.png
│   │       │       ├── six_key_cell_types.png
│   │       │       └── all_cell_types_grid.png
│   │       ├── run2_more_samples/
│   │       ├── run3_higher_lr/
│   │       ├── run4_intermediate/
│   │       ├── run5_final/
│   │       └── comparison_across_runs.png
│   │
│   └── xenium/
│       ├── all_cells_20200625.h5ad      ← auto-downloaded reference (~700MB)
│       ├── snRNA_annotation_astro_subtypes_refined59_20200823.csv
│       ├── run1/
│       │   ├── inf_aver.csv
│       │   ├── adata_xenium_mapped.h5ad
│       │   ├── cell_type_abundances_8ct.png
│       │   └── all_cell_types_59ct.png
│       ├── run2/
│       ├── run3/
│       ├── run4/
│       └── run5/
│
├── project_environment.yml              ← Conda env record (not used at runtime)
└── README.md
```

---

## Prerequisites

- HPC cluster with Singularity installed
- NVIDIA **V100 (Volta) GPU** — Theano backend in cell2location v0.3 is incompatible
  with A100, L40s, and other modern architectures
- Desktop interactive session capability (not just login node access)
- ~25GB free storage total

---

## Environment Setup

### Pull the Singularity Container

Run once. Takes several minutes and produces a ~4.8GB `.sif` file.

```bash
cd /path/to/your/project/data/
singularity pull cell2location.sif docker://quay.io/vitkl/cell2location:latest
```

### Request a GPU Desktop Session

1. Open the SCC OnDemand web interface
2. Request a **Desktop** interactive app
3. Request form settings:
   - **GPU:** Volta/V100 specifically — newer architectures will error
   - **RAM:** 16GB minimum, 32GB recommended for GBmap loading
   - **Time:** Reproduction ~2h | GBM single run ~30min–17h | Xenium Run5 ~3h

> **Order matters:** GPU must be allocated before launching Singularity.
> You cannot add a GPU to an already-running container.

### Launch Container and Jupyter

```bash
# In terminal inside the desktop session:

# 1. Verify GPU is allocated
nvidia-smi   # must show a V100

# 2. Start container with GPU passthrough
singularity exec --nv /path/to/data/cell2location.sif bash

# 3. From inside the container, start Jupyter
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser
```

Copy the URL with token from terminal output and open it in the desktop browser.

### GPU Verification — MANDATORY First Step

```python
# This must be the ABSOLUTE FIRST cell in any notebook or first lines of any script
import os
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32"

import theano
print(theano.config.device)   # must print: cuda
```

> If you import `cell2location` or `pymc3` before setting `THEANO_FLAGS`, Theano
> initializes to CPU and the flag has no effect. You must restart the kernel.

---

## Downloading the Datasets

### Quick Start — Download Everything

```bash
# Make scripts executable
chmod +x shell_scripts/*.sh

# Download all four datasets (~10.6GB total)
bash shell_scripts/download_all.sh

# Or with a custom data directory:
bash shell_scripts/download_all.sh --data-dir /projectnb/yourproject/data

# Mouse brain data only (reproduction, no extensions):
bash shell_scripts/download_all.sh --only-reproduction

# Everything except Xenium:
bash shell_scripts/download_all.sh --skip-xenium
```

The master script runs all individual downloads in order, verifies file integrity,
asks for confirmation before starting, and prints a summary with next steps.

---

### Shell Scripts Reference

#### download_all.sh

Master download script. Orchestrates all other download scripts.

```bash
bash shell_scripts/download_all.sh [OPTIONS]

Options:
  --data-dir PATH         Override default data directory (default: ./data)
  --skip-xenium           Skip Xenium dataset download
  --only-reproduction     Download only mouse brain data (scripts 1 + 2)
  -h, --help              Show help

Estimated total: ~10.6GB | Time: 20–60 min depending on connection
```

#### download_mouse_scrna.sh

Downloads the mouse brain snRNA-seq reference (E-MTAB-11115) from EBI BioStudies.

```bash
bash shell_scripts/download_mouse_scrna.sh [OPTIONS]

Downloads to: data/mousescRNAseq/
Size:         ~2.8GB
Source:       ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/115/E-MTAB-11115/Files/

Key outputs:
  5705STDY8058280_filtered_feature_bc_matrix.h5  (x6 samples)
  cell_annotation.csv  (40,531 cell IDs with annotation_1 labels)

Options:
  --data-dir PATH         Override default data directory
```

**cell_annotation.csv format:**
```
Cell ID,sample,annotation_1,annotation_1_print
5705STDY8058283_AAACCCAAGCCTATTG-1,5705STDY8058283,Ext_L23,22_Ext_L23
```
- `Cell ID` format: `{sampleID}_{barcode}-1`
- `annotation_1` is the cell type label passed to cell2location

#### download_mouse_visium.sh

Downloads the mouse brain Visium data (E-MTAB-11114) from EBI BioStudies and
extracts all five spatial archives.

```bash
bash shell_scripts/download_mouse_visium.sh [OPTIONS]

Downloads to: data/mouseVisium/
Size:         ~5.5GB (skips .cloupe files, ~750MB each)
Source:       ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/114/E-MTAB-11114/Files/

Key outputs:
  ST8059048_filtered_feature_bc_matrix.h5
  ST8059048_spatial/spatial/   ← extracted from tar.gz, has INNER nesting

Options:
  --data-dir PATH         Override default data directory
  --section ST8059048     Download and extract a single section only
```

> **Spatial path nesting:** After extraction, the path is
> `ST8059048_spatial/spatial/` not `ST8059048_spatial/`. The script
> handles this automatically. Script 04 points to the inner `spatial/` directory.

#### download_gbm_data.sh

Downloads both the GBmap scRNA-seq reference (NCBI GEO GSE211376) and the
human glioblastoma Visium spatial data (10x Genomics).

```bash
bash shell_scripts/download_gbm_data.sh [OPTIONS]

Downloads to:
  data/GBmap_reference/        ~305MB (counts + metadata)
  data/humanGlioblastoma/      ~28MB  (count matrix + spatial archive)

Key outputs:
  GBmap_raw_counts.tsv.gz      (27,102 genes x 39,355 cells)
  GBmap_metadata.csv.gz        (cell type labels in predicted.high_hierarchy)
  Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5
  spatial/                     (extracted, flat — no inner nesting unlike mouse)

Options:
  --data-dir PATH         Override default data directory
  --skip-reference        Download Visium only, skip GBmap
  --skip-visium           Download GBmap only, skip Visium
```

After downloading, run:
```bash
python scripts/02_prepare_gbmap_reference.py
```
to build `GBmap_reference_compressed.h5ad` for fast reloading.

**GBmap cell type distribution (14 types):**

| Cell Type | Count | % |
|---|---|---|
| Oligodendrocyte | 12,716 | 32.3% |
| OPC-like | 5,201 | 13.2% |
| TAM-MG (microglia-derived macrophages) | 4,011 | 10.2% |
| TAM-BDM (bone marrow-derived macrophages) | 3,215 | 8.2% |
| MES-like (mesenchymal-like malignant) | 2,790 | 7.1% |
| AC-like (astrocyte-like malignant) | 2,752 | 7.0% |
| NPC-like (neural progenitor-like malignant) | 2,687 | 6.8% |
| OPC-like malignant, Pericyte, Endothelial, CD4/CD8, DC, Astrocyte, Neuron | ~5,983 | ~15.2% |

#### download_xenium.sh

Downloads the 10x Xenium mouse brain coronal subset and extracts the zip archive.

```bash
bash shell_scripts/download_xenium.sh [OPTIONS]

Downloads to: data/xeniumMouseBrain/
Size:         ~2GB (zip archive + extracted)
Source:       10x Genomics CDN (public, no account required)

Key outputs:
  cell_feature_matrix.h5    (36,602 cells x 248 genes)
  cells.csv.gz              (x_centroid, y_centroid per cell)
  gene_panel.json           (275 targeted genes)
  morphology.ome.tif        (fluorescence tissue image)

Options:
  --data-dir PATH         Override default data directory

Dataset specs:
  Cells:          36,602
  Genes:          248 (after QC) / 275 (panel)
  Mean transcripts: 246.7 per cell
  Shared with reference: 247 genes
  Resolution:     ~10 µm (near single-cell)
```

> **Note:** The Xenium reference atlas (~700MB) used by script 07 is downloaded
> automatically by the Python script, not by this shell script.

---

### Manual Download Instructions

If the shell scripts fail (network issues, different cluster), download manually:

**Mouse brain snRNA-seq:**
```bash
mkdir -p data/mousescRNAseq && cd data/mousescRNAseq
wget -r --no-parent --no-host-directories --cut-dirs=7 \
  ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/115/E-MTAB-11115/Files/
```

**Mouse brain Visium + extract:**
```bash
mkdir -p data/mouseVisium && cd data/mouseVisium
wget -r --no-parent --no-host-directories --cut-dirs=7 \
  ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/114/E-MTAB-11114/Files/

# Extract all spatial archives
for f in *.tar.gz; do
    sample="${f%%_spatial.tar.gz}"
    mkdir -p "${sample}_spatial"
    tar -xzf "$f" -C "${sample}_spatial"
done
```

**GBmap reference:**
```bash
mkdir -p data/GBmap_reference && cd data/GBmap_reference
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211376/suppl/GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz" -O GBmap_raw_counts.tsv.gz
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211376/suppl/GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz" -O GBmap_metadata.csv.gz
```

**GBM Visium:**
```bash
mkdir -p data/humanGlioblastoma && cd data/humanGlioblastoma
wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_spatial.tar.gz
tar -xzf Parent_Visium_Human_Glioblastoma_spatial.tar.gz
```

**Xenium:**
```bash
mkdir -p data/xeniumMouseBrain && cd data/xeniumMouseBrain
wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
unzip Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
```

---

## Pipeline Overview

```
PART 1: REPRODUCTION (mouse brain)
───────────────────────────────────────────────────────────────────────
download_mouse_scrna.sh ──► data/mousescRNAseq/
download_mouse_visium.sh ──► data/mouseVisium/ (+ extraction)
         │
         ▼
03_reproduction_scrna_umap.py
  Loads 6 h5 files, merges annotations, BBKNN UMAP
  ──► adata_snrna_processed.h5ad
  ──► figure_2b_umap.png
         │
         ▼
04_reproduction_visium_mapping.py
  Regression model (cell type signatures) +
  Spatial model (30k iterations, V100 GPU)
  ──► figure_2c_HE.png
  ──► figure_2d_major.png
  ──► figure_2e_sparse.png
  ──► adata_vis_mapped.h5ad

PART 2: GBM EXTENSION
───────────────────────────────────────────────────────────────────────
download_gbm_data.sh ──► data/GBmap_reference/ + data/humanGlioblastoma/
         │
         ▼
02_prepare_gbmap_reference.py
  Loads raw TSV + metadata, QC filtering, saves h5ad
  ──► GBmap_reference_compressed.h5ad
         │
         ├──► 05_gbm_spatial_mapping.py  (single config)
         │         ──► results/human_gbm/single_run/
         │
         └──► 06_gbm_parameter_sweep.py  (all 5 configs)
                   ──► results/human_gbm/parameter_runs/run1-5/
                   ──► comparison_across_runs.png

PART 3: XENIUM EXTENSION
───────────────────────────────────────────────────────────────────────
download_xenium.sh ──► data/xeniumMouseBrain/
         │
         ▼
07_xenium_spatial_mapping.py
  Auto-downloads reference atlas (~700MB)
  Gene intersection (247 shared), QC filtering
  Regression + spatial model per configuration
  ──► results/xenium/run1-5/
      cell_type_abundances_8ct.png
      all_cell_types_59ct.png
```

---

## Scripts Reference

### 00_gpu_setup.py

Sets Theano GPU flag. Must execute before any theano/pymc3/cell2location import.

```bash
python scripts/00_gpu_setup.py   # verify GPU detection
```

### 01_download_gbmap_reference.py

Python alternative to `download_gbm_data.sh --skip-visium` for downloading
the GBmap reference. Includes MD5 verification.

```bash
python scripts/01_download_gbmap_reference.py [--ref-dir PATH] [--force]
```

### 02_prepare_gbmap_reference.py

Loads raw GBmap files, creates AnnData, applies QC, saves compressed h5ad.

```bash
python scripts/02_prepare_gbmap_reference.py
python scripts/02_prepare_gbmap_reference.py --min-genes 200 --min-cells 3 --force
```

**Runtime:** ~5 minutes. **Output:** `GBmap_reference_compressed.h5ad`

### 03_reproduction_scrna_umap.py

Loads all 6 mouse snRNA-seq samples, merges annotations, runs BBKNN + UMAP,
produces Figure 2b.

```bash
python scripts/03_reproduction_scrna_umap.py
python scripts/03_reproduction_scrna_umap.py --n-hvg 3000 --n-pcs 50
```

**Runtime:** ~20–40 minutes. **Output:** `adata_snrna_processed.h5ad`, `figure_2b_umap.png`

**Important:** Raw counts are preserved in `layers["counts"]` before normalization.
Script 04 uses these raw counts as input to cell2location.

### 04_reproduction_visium_mapping.py

Trains regression model → loads Visium → runs spatial mapping → Figures 2c–2e.

Must run inside Singularity container with `--nv` GPU flag.

```bash
# Quick test (3,000 iterations):
python scripts/04_reproduction_visium_mapping.py --n-iter-spatial 3000

# Full run matching paper (30,000 iterations, ~25 min on V100):
python scripts/04_reproduction_visium_mapping.py --n-iter-spatial 30000

# All options:
python scripts/04_reproduction_visium_mapping.py \
    --n-iter-reg 100 \
    --n-iter-spatial 30000 \
    --cells-per-spot 8 \
    --section ST8059048
```

**Key paper parameters:**

| Parameter | Value |
|---|---|
| `cells_per_spot` | 8 |
| `alpha_mean` (detection_alpha) | 200 |
| spatial iterations | 30,000 |
| shared genes | 10,085 |
| learning rate | 0.001 |

**Known v0.3 API bug:** `covariate_col_names=[]` causes an error.
Fix (already in script): add `adata.obs["dummy_covar"] = "1"` and pass
`covariate_col_names=["dummy_covar"]`.

### 05_gbm_spatial_mapping.py

Single-configuration GBM extension. Use for testing before committing to the full sweep.

```bash
# Test run (~25 min):
python scripts/05_gbm_spatial_mapping.py \
    --n-iter-reg 100 --n-iter-spat 500 --n-comb 20 --run-name test

# Final quality (~17 hr on V100):
python scripts/05_gbm_spatial_mapping.py \
    --n-iter-reg 2000 --n-iter-spat 10000 --n-comb 30 \
    --n-samples 500 --cells-per-spot 15 --run-name final
```

`--cells-per-spot 15` is higher than reproduction (8) because glioblastoma
tissue has higher cellularity.

### 06_gbm_parameter_sweep.py

Runs all 5 GBM configurations sequentially. Loads reference and Visium once,
reuses across runs.

```bash
# Full sweep (asks for confirmation, shows time estimate):
python scripts/06_gbm_parameter_sweep.py

# Specific runs only:
python scripts/06_gbm_parameter_sweep.py --runs run1_test run5_final

# Preview without running:
python scripts/06_gbm_parameter_sweep.py --dry-run
```

**Configurations:**

| Run | Reg | Spat | n_comb | LR | Samples | Est. time |
|---|---|---|---|---|---|---|
| run1_test | 100 | 500 | 20 | 0.005 | 200 | ~25 min |
| run2_more_samples | 100 | 500 | 20 | 0.005 | 500 | ~30 min |
| run3_higher_lr | 100 | 500 | 20 | 0.025 | 200 | ~25 min |
| run4_intermediate | 500 | 5,000 | 30 | 0.005 | 500 | ~3.5 hr |
| run5_final | 2,000 | 10,000 | 30 | 0.005 | 500 | ~17 hr |

### 07_xenium_spatial_mapping.py

Xenium extension. Auto-downloads mouse brain reference atlas. Runs 1 or all 5 configs.

```bash
# Single quick test (~20 min on V100):
python scripts/07_xenium_spatial_mapping.py \
    --n-iter-reg 100 --n-iter-spat 500 --n-comb 20 --run-name run1

# Final quality run (~3 hr on V100):
python scripts/07_xenium_spatial_mapping.py \
    --n-iter-reg 4000 --n-iter-spat 20000 --n-comb 50 --run-name run5

# Full 5-configuration sweep:
python scripts/07_xenium_spatial_mapping.py --sweep

# Custom train/test split (Run 4):
python scripts/07_xenium_spatial_mapping.py \
    --train-split 0.7 --n-samples 400 --batch-size 400 --run-name run4
```

**Configurations:**

| Run | Reg | Spat | n_comb | LR | Split | Samples |
|---|---|---|---|---|---|---|
| run1 | 100 | 500 | 20 | 0.005 | 90/10 | 200 |
| run2 | 500 | 2,500 | 20 | 0.005 | 90/10 | 200 |
| run3 | 100 | 500 | 20 | 0.025 | 90/10 | 200 |
| run4 | 100 | 500 | 20 | 0.005 | 70/30 | 400 |
| run5 | 4,000 | 20,000 | 50 | 0.005 | 90/10 | 200 |

### utils.py

Shared helpers used across all pipeline scripts. Import as needed:

```python
from scripts.utils import (
    check_files,                 # verify required files exist before running
    ensure_dir,                  # create directories safely
    extract_signatures,          # pull cell type matrix from trained regression model
    extract_signatures_xenium,   # Xenium-specific variant
    add_spot_factors,            # attach spot_factors to AnnData.obs
    plot_spatial_panel,          # multi-panel spatial abundance grid (Visium or Xenium)
    plot_dominant_cell_type,     # dominant cell type per spot map
    Timer,                       # wall-clock timing
    print_section,               # formatted section headers
    print_summary,               # formatted key-value summary
)
```

---

## Running the Full Pipeline

### Reproduction Only

```bash
# 1. Download data
bash shell_scripts/download_mouse_scrna.sh
bash shell_scripts/download_mouse_visium.sh

# 2. Pull Singularity container (once)
singularity pull data/cell2location.sif docker://quay.io/vitkl/cell2location:latest

# 3. Start desktop session with V100 GPU, then:
singularity exec --nv data/cell2location.sif bash
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser

# 4. In Jupyter (or as scripts with GPU):
python scripts/03_reproduction_scrna_umap.py
python scripts/04_reproduction_visium_mapping.py --n-iter-spatial 30000
```

### GBM Extension

```bash
# 1. Download data
bash shell_scripts/download_gbm_data.sh

# 2. Prepare reference (no GPU needed — run outside container)
python scripts/02_prepare_gbmap_reference.py

# 3. Inside container (GPU required):
# Quick test:
python scripts/05_gbm_spatial_mapping.py \
    --n-iter-reg 100 --n-iter-spat 500 --run-name test

# Full sweep:
python scripts/06_gbm_parameter_sweep.py
```

### Xenium Extension

```bash
# 1. Download data
bash shell_scripts/download_xenium.sh

# 2. Inside container (GPU required):
# Quick test:
python scripts/07_xenium_spatial_mapping.py \
    --n-iter-reg 100 --n-iter-spat 500 --n-comb 20 --run-name run1

# Full sweep:
python scripts/07_xenium_spatial_mapping.py --sweep
```

---

## Parameter Reference

### Core Cell2location Model Parameters

| Parameter | Mouse Brain | GBM | Xenium |
|---|---|---|---|
| `cells_per_spot` (N̂) | 8 | 15 | 1–2 |
| `factors_per_spot` | 7 | 10 | 7 |
| `combs_per_spot` | 2.5 | 3 | 2.5 |
| `alpha_mean` | 200 | 200 | 200 |
| Spatial iterations (full) | 30,000 | 10,000 | 20,000 |
| Learning rate | 0.001 | 0.005 | 0.005 |
| Shared genes | 10,085 | 25,731 | 247 |
| Cell types | 59 | 14 | 59 |

`cells_per_spot` is the most important parameter to adjust per tissue:
- Derive from H&E histology by manually counting nuclei in 10–20 spots
- Tumor tissue (GBM) has higher cellularity → higher prior
- Xenium is near single-cell resolution → prior of ~1

### Gene Filtering Parameters

Three parameters in the regression model control which genes are retained:

| Parameter | Description | Effect of Increasing |
|---|---|---|
| `cell_count_cutoff` | Min cells a gene must appear in | Stricter — removes sparse genes |
| `cell_percentage_cutoff2` | Min fraction of cells expressing gene | Stricter — removes lowly expressed |
| `nonz_mean_cutoff` | Min mean expression among expressing cells | Removes weakly expressed |

In cell2location v0.3 these are applied via the `>5% expression` filter in script 03.

---

## Expected Outputs

### Reproduction

| Figure | File | Description |
|---|---|---|
| 2b | `figure_2b_umap.png` | UMAP of 59 cell subtypes, colored by subtype and broad class |
| 2c | `figure_2c_HE.png` | H&E histology of Visium section ST8059048 |
| 2d | `figure_2d_major.png` | Spatial abundance of 6 major regional cell types |
| 2e | `figure_2e_sparse.png` | Spatial abundance of 3 sparse inhibitory neuron subtypes |

The Allen Mouse Brain Atlas anatomical region overlay from original Figure 2c is
omitted — it requires tissue registration outside the scope of cell2location.

Training duration effect: at 3,000 iterations Excitatory L2/3 shows diffuse
cortex-wide signal. At 30,000 iterations it sharpens to the cortical ribbon only,
demonstrating posterior shrinkage with training.

### GBM Extension

| File | Description |
|---|---|
| `dominant_cell_type_map.png` | Highest-abundance cell type per spot |
| `six_key_cell_types.png` | MES-like, Astrocyte, CD4/CD8, OPC, AC-like, Pericyte |
| `all_cell_types_grid.png` | All 14 cell types in a grid |
| `comparison_across_runs.png` | Cell type composition across all 5 configurations |

MES-like cells were consistently identified as the dominant malignant state across
all five configurations. Test run (100/500 iter) matched top-3 rankings of the
final run (2,000/10,000 iter), indicating broad spatial patterns are robust.

### Xenium Extension

| File | Description |
|---|---|
| `cell_type_abundances_8ct.png` | 8 selected cell types with viridis colormap scatter |
| `all_cell_types_59ct.png` | All 59 cell types with magma colormap scatter |

Run5 (4,000/20,000 iterations, n_comb=50) produced the smoothest gradients.
Run4 (70/30 split) produced rectangular artifact patches across all 59 types —
a data handling artifact, not biology.

---

## Known Issues and Limitations

**GPU architecture:** Theano (cell2location v0.3) is incompatible with A100, L40s,
and other modern GPUs. Must use V100 specifically.

**GPU flag order:** `THEANO_FLAGS` must be set before any theano/pymc3/cell2location
import. Forgetting requires a kernel restart.

**`covariate_col_names` bug:** Passing `[]` or `[None]` to `run_regression` raises
an error in v0.3. Fixed in all scripts with a dummy constant covariate.

**Combined cell type figures:** The original paper overlays multiple cell types in
one spatial image. Cell2location v0.3 lacks this functionality and the code was not
found in the authors' GitHub. All figures here show each cell type separately.

**Storage:** Cell2location v0.1.3 dependencies exceed 10GB, requiring the v0.3
container. Full pipeline storage (container + data + results) is ~25GB.

**Session persistence:** Desktop session timeouts lose all in-memory variables.
Checkpoints saved: `adata_snrna_processed.h5ad` (script 03), `GBmap_reference_compressed.h5ad`
(script 02), `adata_vis_mapped.h5ad` (script 04), `adata_gbm_mapped.h5ad` (scripts 05/06),
`adata_xenium_mapped.h5ad` (script 07).

**Xenium Run4 artifacts:** The 70/30 train/test split configuration produced
rectangular artifact patches across all 59 cell types. Run5 is the recommended
high-quality configuration.

---

## References

Kleshchevnikov, V., Shmatko, A., Dann, E. et al. Cell2location maps fine-grained
cell types in spatial transcriptomics. *Nat Biotechnol* 40, 661–671 (2022).
https://doi.org/10.1038/s41587-021-01139-4

Ruiz-Moreno, C., et al. Harmonized single-cell atlas of glioblastoma.
GEO: GSE211376 (2022).

Neftel, C., Laffy, J., Filbin, M.G. et al. An Integrative Model of Cellular
States, Plasticity, and Genetics for Glioblastoma. *Cell* 178(4):835-849.e21 (2019).

10x Genomics. Parent Visium Human Glioblastoma dataset (2020).
https://www.10xgenomics.com/datasets

10x Genomics. Xenium V1 FF Mouse Brain Coronal Subset (CTX + HP) (2023).
https://www.10xgenomics.com/datasets/xenium-ff-mouse-brain-coronal-subset-ctx-hp-1-0-2
