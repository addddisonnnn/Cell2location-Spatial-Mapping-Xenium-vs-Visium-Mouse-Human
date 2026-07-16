"""
02_prepare_gbmap_reference.py
------------------------------
Loads the raw GBmap files downloaded by 01_download_gbmap_reference.py,
creates an AnnData object, applies QC filtering, and saves a compressed
.h5ad file ready for use in the human GBM cell2location pipeline.

Inputs (from REF_DIR):
    GBmap_raw_counts.tsv.gz   — raw count matrix
    GBmap_metadata.csv.gz     — cell metadata

Outputs (saved to REF_DIR):
    GBmap_reference_compressed.h5ad  — processed AnnData (cells x genes)
                                       with cell_type column from
                                       predicted.high_hierarchy

Run once before 05_gbm_regression_model.py or 06_gbm_spatial_mapping.py.

Usage:
    python scripts/02_prepare_gbmap_reference.py
    python scripts/02_prepare_gbmap_reference.py --ref-dir /path/to/ref
"""

import os
import time
import argparse
import warnings

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Prepare GBmap reference h5ad")
parser.add_argument(
    "--ref-dir",
    default=os.path.join(os.path.dirname(__file__), "..", "data", "GBmap_reference"),
    help="Directory containing GBmap raw files (default: data/GBmap_reference/)"
)
parser.add_argument(
    "--min-genes", type=int, default=200,
    help="Min genes per cell (default: 200)"
)
parser.add_argument(
    "--min-cells", type=int, default=3,
    help="Min cells per gene (default: 3)"
)
parser.add_argument(
    "--force",
    action="store_true",
    help="Re-build even if compressed h5ad already exists"
)
args = parser.parse_args()

REF_DIR = os.path.abspath(args.ref_dir)
COUNTS_FILE   = os.path.join(REF_DIR, "GBmap_raw_counts.tsv.gz")
METADATA_FILE = os.path.join(REF_DIR, "GBmap_metadata.csv.gz")
OUTPUT_H5AD   = os.path.join(REF_DIR, "GBmap_reference_compressed.h5ad")

# ── Main ──────────────────────────────────────────────────────────────────────
print("=" * 60)
print("PREPARING GBMAP REFERENCE")
print("=" * 60)

if os.path.exists(OUTPUT_H5AD) and not args.force:
    print(f"Already exists: {OUTPUT_H5AD}")
    print("Loading to verify... ", end="")
    adata = sc.read(OUTPUT_H5AD)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes, {adata.obs['cell_type'].nunique()} cell types")
    print("Use --force to rebuild.")
    exit(0)

# ── Load metadata ─────────────────────────────────────────────────────────────
print(f"\nLoading metadata from {METADATA_FILE}...")
metadata = pd.read_csv(METADATA_FILE, compression="gzip")
print(f"  Shape: {metadata.shape}")
print(f"  Columns: {metadata.columns.tolist()}")

print("\nCell type distribution (predicted.high_hierarchy):")
ct_counts = metadata["predicted.high_hierarchy"].value_counts()
for ct, count in ct_counts.items():
    print(f"  {ct:30s}: {count:6,} ({count/len(metadata)*100:.1f}%)")
print(f"  Total: {len(metadata):,} cells, {metadata['patient'].nunique()} patients")

# ── Load count matrix ─────────────────────────────────────────────────────────
print(f"\nLoading count matrix from {COUNTS_FILE}...")
print("  (This takes 2-5 minutes — matrix is ~27k genes x ~39k cells)")
t0 = time.time()
counts = pd.read_csv(COUNTS_FILE, compression="gzip", sep="\t", index_col=0)
print(f"  Shape: {counts.shape[0]:,} genes x {counts.shape[1]:,} cells")
print(f"  Load time: {time.time()-t0:.1f}s")

# ── Build AnnData (cells x genes) ─────────────────────────────────────────────
print("\nBuilding AnnData object...")
adata = ad.AnnData(
    X=counts.T.values.astype(np.int32),  # transpose: cells as rows
    obs=metadata,
    var=pd.DataFrame(index=counts.index),
)
adata.obs["cell_type"] = adata.obs["predicted.high_hierarchy"]
print(f"  Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

# Free the counts dataframe
del counts

# ── QC filtering ──────────────────────────────────────────────────────────────
print(f"\nQC filtering (min_genes={args.min_genes}, min_cells={args.min_cells})...")
print(f"  Before: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
sc.pp.filter_cells(adata, min_genes=args.min_genes)
print(f"  After cell filter: {adata.n_obs:,} cells")
sc.pp.filter_genes(adata, min_cells=args.min_cells)
print(f"  After gene filter: {adata.n_vars:,} genes")
print(f"  Final shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

# ── Cell type composition plot ────────────────────────────────────────────────
print("\nSaving cell type composition plot...")
ct_counts_filtered = adata.obs["cell_type"].value_counts()
fig, ax = plt.subplots(figsize=(14, 6))
colors = plt.cm.Set3(range(len(ct_counts_filtered)))
ax.bar(range(len(ct_counts_filtered)), ct_counts_filtered.values, color=colors)
ax.set_xticks(range(len(ct_counts_filtered)))
ax.set_xticklabels(ct_counts_filtered.index, rotation=45, ha="right", fontsize=10)
ax.set_ylabel("Number of cells")
ax.set_title("GBmap Reference: Cell Type Distribution (after QC)")
plt.tight_layout()
plot_path = os.path.join(REF_DIR, "reference_cell_type_distribution.png")
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {plot_path}")

# ── Save ──────────────────────────────────────────────────────────────────────
print(f"\nSaving compressed h5ad to {OUTPUT_H5AD}...")
adata.write(OUTPUT_H5AD)
size_mb = os.path.getsize(OUTPUT_H5AD) / 1024 / 1024
print(f"  Done: {size_mb:.1f} MB")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  Cells:       {adata.n_obs:,}")
print(f"  Genes:       {adata.n_vars:,}")
print(f"  Cell types:  {adata.obs['cell_type'].nunique()}")
print(f"  Patients:    {adata.obs['patient'].nunique()}")
print(f"  Output:      {OUTPUT_H5AD}")
print("\nReady for 05_gbm_regression_model.py")
print("=" * 60)
