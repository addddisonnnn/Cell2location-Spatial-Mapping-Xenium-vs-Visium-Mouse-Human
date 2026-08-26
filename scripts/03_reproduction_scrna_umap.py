"""
03_reproduction_scrna_umap.py
------------------------------
Loads the mouse brain snRNA-seq reference (E-MTAB-11115), merges cell type
annotations, preprocesses the data, and generates Figure 2b from the original
Cell2location paper: a UMAP of all 59 annotated cell subtypes colored by
fine-grained subtype and broad class.

Inputs:
    data/mousescRNAseq/
        5705STDY8058280_filtered_feature_bc_matrix.h5  (x6 samples)
        cell_annotation.csv

Outputs:
    results/reproduction/
        adata_snrna_processed.h5ad   — processed AnnData with UMAP coords
                                       and annotation columns; used by
                                       04_reproduction_regression_model.py
        figure_2b_umap.png           — reproduction of paper Figure 2b

Usage:
    python scripts/03_reproduction_scrna_umap.py
    python scripts/03_reproduction_scrna_umap.py --data-dir /path/to/data
"""

import os
import sys
import argparse
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import bbknn

warnings.filterwarnings("ignore")

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Load mouse snRNA-seq and produce Figure 2b UMAP")
parser.add_argument(
    "--data-dir",
    default=os.path.join(os.path.dirname(__file__), "..", "data"),
    help="Base data directory (default: data/)"
)
parser.add_argument(
    "--n-hvg", type=int, default=3000,
    help="Number of highly variable genes for UMAP (default: 3000)"
)
parser.add_argument(
    "--n-pcs", type=int, default=50,
    help="Number of PCA components (default: 50)"
)
args = parser.parse_args()

DATA_DIR    = os.path.abspath(args.data_dir)
SCRNA_DIR   = os.path.join(DATA_DIR, "mousescRNAseq")
RESULTS_DIR = os.path.join(DATA_DIR, "..", "results", "reproduction")
os.makedirs(RESULTS_DIR, exist_ok=True)

OUTPUT_H5AD = os.path.join(RESULTS_DIR, "adata_snrna_processed.h5ad")
OUTPUT_FIG  = os.path.join(RESULTS_DIR, "figure_2b_umap.png")

SCRNA_SAMPLES = [
    "5705STDY8058280",
    "5705STDY8058281",
    "5705STDY8058282",
    "5705STDY8058283",
    "5705STDY8058284",
    "5705STDY8058285",
]

# ── Load all 6 snRNA-seq samples ──────────────────────────────────────────────
print("=" * 60)
print("LOADING snRNA-seq DATA (E-MTAB-11115)")
print("=" * 60)

adatas = []
for sample in SCRNA_SAMPLES:
    path = os.path.join(SCRNA_DIR, f"{sample}_filtered_feature_bc_matrix.h5")
    if not os.path.exists(path):
        print(f"ERROR: Missing file: {path}")
        sys.exit(1)
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    adata.obs["sample"] = sample
    # Format cell IDs to match annotation file: {sampleID}_{barcode}
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
    adatas.append(adata)
    print(f"  Loaded {sample}: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

adata_snrna = sc.concat(adatas)
print(f"\nCombined: {adata_snrna.n_obs:,} cells x {adata_snrna.n_vars:,} genes")

# ── Load and merge annotations ────────────────────────────────────────────────
print("\nLoading cell annotations...")
annot_path = os.path.join(SCRNA_DIR, "cell_annotation.csv")
cell_annot  = pd.read_csv(annot_path, index_col="Cell ID")
print(f"  Annotation file: {cell_annot.shape[0]:,} cells, columns: {cell_annot.columns.tolist()}")

# Preserve raw counts before any normalization — required for cell2location
adata_snrna.layers["counts"] = adata_snrna.X.copy()

# Merge annotations (cells without annotation are dropped)
overlap = adata_snrna.obs_names.intersection(cell_annot.index)
adata_snrna.obs["annotation_1"]       = cell_annot.loc[overlap, "annotation_1"]
adata_snrna.obs["annotation_1_print"] = cell_annot.loc[overlap, "annotation_1_print"]
adata_snrna = adata_snrna[adata_snrna.obs["annotation_1"].notna()].copy()

print(f"\nAfter annotation merge: {adata_snrna.n_obs:,} cells  (expected ~40,531)")
print(f"Cell types (annotation_1): {adata_snrna.obs['annotation_1'].nunique()} unique")

# ── Broad class mapping ───────────────────────────────────────────────────────
unique_cts = adata_snrna.obs["annotation_1"].unique()
broad_map = {}
for ct in unique_cts:
    if "Astro" in ct:
        broad_map[ct] = "Astrocytes"
    elif "Ext" in ct:
        broad_map[ct] = "Excitatory"
    elif "Inh" in ct:
        broad_map[ct] = "Inhibitory"
    elif "Oligo" in ct or "OPC" in ct:
        broad_map[ct] = "Oligo-OPC"
    elif "Micro" in ct:
        broad_map[ct] = "Micro"
    elif "Endo" in ct:
        broad_map[ct] = "Endo"
    else:
        broad_map[ct] = "Other"

adata_snrna.obs["broad_class"] = adata_snrna.obs["annotation_1"].map(broad_map).fillna("Other")
print(f"Broad classes: {adata_snrna.obs['broad_class'].value_counts().to_dict()}")

# ── Preprocessing for UMAP ────────────────────────────────────────────────────
# Note: this normalizes adata_snrna.X for UMAP only.
# Raw counts are preserved in layers["counts"] for cell2location.
print(f"\nPreprocessing for UMAP ({args.n_hvg} HVGs, {args.n_pcs} PCs, BBKNN batch correction)...")
sc.pp.normalize_total(adata_snrna, target_sum=1e4)
sc.pp.log1p(adata_snrna)
sc.pp.highly_variable_genes(adata_snrna, n_top_genes=args.n_hvg, batch_key="sample")
sc.pp.scale(adata_snrna, max_value=10)
sc.tl.pca(adata_snrna, n_comps=args.n_pcs)
bbknn.bbknn(adata_snrna, batch_key="sample", neighbors_within_batch=3)
sc.tl.umap(adata_snrna)
print("  UMAP done.")

# ── Figure 2b ─────────────────────────────────────────────────────────────────
print(f"\nGenerating Figure 2b...")
fig, axes = plt.subplots(1, 2, figsize=(22, 8))

sc.pl.umap(
    adata_snrna,
    color="annotation_1",
    legend_loc="on data",
    legend_fontsize=5,
    title="Cell subtypes (59)",
    ax=axes[0],
    show=False,
)
sc.pl.umap(
    adata_snrna,
    color="broad_class",
    legend_loc="right margin",
    title="Broad class",
    ax=axes[1],
    show=False,
)
plt.suptitle("Figure 2b Reproduction — Mouse Brain snRNA-seq UMAP", fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {OUTPUT_FIG}")

# ── Save processed AnnData ────────────────────────────────────────────────────
print(f"\nSaving processed AnnData...")
# Restore raw counts to .X so cell2location receives correct input
adata_snrna.X = adata_snrna.layers["counts"].copy()
adata_snrna.write_h5ad(OUTPUT_H5AD)
size_mb = os.path.getsize(OUTPUT_H5AD) / 1024 / 1024
print(f"  Saved: {OUTPUT_H5AD} ({size_mb:.1f} MB)")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  Cells:       {adata_snrna.n_obs:,}")
print(f"  Genes:       {adata_snrna.n_vars:,}")
print(f"  Cell types:  {adata_snrna.obs['annotation_1'].nunique()}")
print(f"  UMAP:        {OUTPUT_FIG}")
print(f"  AnnData:     {OUTPUT_H5AD}")
print("\nRun 04_reproduction_visium_mapping.py next.")
print("=" * 60)
