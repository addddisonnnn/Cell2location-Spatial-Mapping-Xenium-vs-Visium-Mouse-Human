"""
05_gbm_spatial_mapping.py
--------------------------
Loads the GBmap reference (built by 02_prepare_gbmap_reference.py) and the
human glioblastoma Visium data, trains the regression model to learn cell type
signatures, runs spatial deconvolution, and saves results.

This script runs a single parameter configuration. To run the full parameter
sweep across five configurations, run 06_gbm_parameter_sweep.py instead.

IMPORTANT: Must be run inside the Singularity container with GPU (V100).

Inputs:
    data/GBmap_reference/GBmap_reference_compressed.h5ad  (from script 02)
    data/humanGlioblastoma/
        Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5
        spatial/

Outputs:
    results/human_gbm/
        inf_aver.csv                    — reference signature matrix
        adata_gbm_mapped.h5ad           — Visium AnnData with abundances
        cell_abundance_summary.csv      — per-cell-type summary statistics
        spatial_maps_gbm_states.png     — spatial maps of malignant + TME types
        dominant_cell_type_map.png      — dominant cell type per spot
        cell_type_composition.png       — bar chart of total abundances

Usage:
    # Quick test (fast, low quality):
    python scripts/05_gbm_spatial_mapping.py --n-iter-reg 100 --n-iter-spat 500 --n-comb 20

    # Final run (paper quality):
    python scripts/05_gbm_spatial_mapping.py --n-iter-reg 2000 --n-iter-spat 10000 --n-comb 30
"""

# GPU must be set before any theano import
import os
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32"

import sys
import argparse
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import torch
from sklearn.cluster import KMeans

import cell2location
from cell2location.models import RegressionGeneBackgroundCoverageTorch, LocationModelLinearDependentW

warnings.filterwarnings("ignore")

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Run cell2location on human GBM Visium")
parser.add_argument("--data-dir", default=os.path.join(os.path.dirname(__file__), "..", "data"))
parser.add_argument("--n-iter-reg",  type=int, default=100,   help="Regression iterations (default 100)")
parser.add_argument("--n-iter-spat", type=int, default=500,   help="Spatial iterations (default 500)")
parser.add_argument("--n-comb",      type=int, default=20,    help="n_comb parameter (default 20)")
parser.add_argument("--lr",          type=float, default=0.005, help="Learning rate (default 0.005)")
parser.add_argument("--n-samples",   type=int, default=200,   help="Posterior samples (default 200)")
parser.add_argument("--cells-per-spot", type=int, default=15, help="Cell number prior (default 15)")
parser.add_argument("--run-name",    default="single_run",    help="Name for output subdirectory")
args = parser.parse_args()

DATA_DIR    = os.path.abspath(args.data_dir)
REF_DIR     = os.path.join(DATA_DIR, "GBmap_reference")
GBM_DIR     = os.path.join(DATA_DIR, "humanGlioblastoma")
RESULTS_DIR = os.path.join(DATA_DIR, "..", "results", "human_gbm", args.run_name)
os.makedirs(RESULTS_DIR, exist_ok=True)

use_cuda = torch.cuda.is_available()
print(f"GPU available: {use_cuda}")
if use_cuda:
    print(f"  GPU: {torch.cuda.get_device_name(0)}")

# ── Load reference ────────────────────────────────────────────────────────────
print("=" * 60)
print("LOADING GBMAP REFERENCE")
print("=" * 60)

compressed_ref = os.path.join(REF_DIR, "GBmap_reference_compressed.h5ad")
if not os.path.exists(compressed_ref):
    print(f"ERROR: {compressed_ref} not found.")
    print("Run 02_prepare_gbmap_reference.py first.")
    sys.exit(1)

adata_ref = sc.read(compressed_ref)
print(f"  Loaded: {adata_ref.n_obs:,} cells x {adata_ref.n_vars:,} genes")
print(f"  Cell types: {adata_ref.obs['cell_type'].nunique()}")

# ── Load Visium spatial data ──────────────────────────────────────────────────
print("\nLoading GBM Visium spatial data...")
adata_spatial = sc.read_visium(
    GBM_DIR,
    count_file="Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5",
    load_images=True,
)
adata_spatial.var_names_make_unique()
adata_ref.var_names_make_unique()
print(f"  Loaded: {adata_spatial.n_obs:,} spots x {adata_spatial.n_vars:,} genes")

# ── Gene alignment ────────────────────────────────────────────────────────────
common_genes = np.intersect1d(adata_spatial.var_names, adata_ref.var_names)
print(f"\nShared genes: {len(common_genes):,} ({len(common_genes)/adata_spatial.n_vars*100:.1f}% of spatial)")
adata_ref     = adata_ref[:, common_genes].copy()
adata_spatial = adata_spatial[:, common_genes].copy()
adata_spatial.layers["counts"] = adata_spatial.X.copy()

# ── Regression model ──────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"REGRESSION MODEL (n_iter={args.n_iter_reg}, lr={args.lr})")
print("=" * 60)

adata_ref.obs["dummy_covar"] = "1"
cell2covar = pd.DataFrame({
    "patient":      adata_ref.obs["patient"].values,
    "cell_type":    adata_ref.obs["cell_type"].values,
    "dummy_covar":  adata_ref.obs["dummy_covar"].values,
}, index=adata_ref.obs_names)

X_ref = adata_ref.X.toarray() if hasattr(adata_ref.X, "toarray") else np.array(adata_ref.X)
X_ref = X_ref.astype(np.float32)

mod_ref = RegressionGeneBackgroundCoverageTorch(
    sample_id="patient",
    cell2covar=cell2covar,
    X_data=X_ref,
    data_type="float32",
    n_iter=args.n_iter_reg,
    learning_rate=args.lr,
    total_grad_norm_constraint=200,
    verbose=True,
    use_cuda=use_cuda,
    minibatch_size=2500,
    use_average_as_initial_value=True,
)
mod_ref.fit_advi_iterative(n=1, n_type="restart", n_iter=args.n_iter_reg, train_proportion=0.9)
mod_ref.sample_posterior(node="all", n_samples=500, save_samples=False, mean_field_slot="init_1")

inf_aver_raw    = mod_ref.samples["post_sample_means"]["gene_factors"]
fact_names      = mod_ref.fact_names
ct_indices      = [i for i, f in enumerate(fact_names) if f.startswith("cell_type_")]
ct_names        = [f.replace("cell_type_", "") for f in fact_names if f.startswith("cell_type_")]
inf_aver        = pd.DataFrame(inf_aver_raw[ct_indices, :].T, index=adata_ref.var_names, columns=ct_names)

inf_aver.to_csv(os.path.join(RESULTS_DIR, "inf_aver.csv"))
print(f"  Signatures: {inf_aver.shape}  (saved to {RESULTS_DIR}/inf_aver.csv)")

# ── Spatial model ─────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"SPATIAL MODEL (n_iter={args.n_iter_spat}, n_comb={args.n_comb})")
print("=" * 60)

common_sp = np.intersect1d(adata_spatial.var_names, inf_aver.index)
inf_aver_al = inf_aver.loc[common_sp, :].copy()
X_sp = adata_spatial[:, common_sp].X.toarray() if hasattr(adata_spatial.X, "toarray") else np.array(adata_spatial[:, common_sp].X)
X_sp = X_sp.astype(np.float32)

mod_spatial = LocationModelLinearDependentW(
    cell_state_mat=inf_aver_al.values.astype(np.float32),
    X_data=X_sp,
    n_comb=args.n_comb,
    data_type="float32",
    n_iter=args.n_iter_spat,
    learning_rate=args.lr,
    total_grad_norm_constraint=200,
    verbose=True,
    cell_number_prior={"cells_per_spot": args.cells_per_spot, "factors_per_spot": 10, "combs_per_spot": 3},
)
mod_spatial.fit_advi_iterative()
mod_spatial.sample_posterior(
    node="all", n_samples=args.n_samples, save_samples=False,
    mean_field_slot="init_1", batch_size=args.n_samples,
)

spot_factors = mod_spatial.samples["post_sample_means"]["spot_factors"]
cell_types   = inf_aver_al.columns.tolist()
abund_cols   = []
for i, ct in enumerate(cell_types):
    col = f"{ct}_mean"
    adata_spatial.obs[col] = spot_factors[:, i]
    abund_cols.append(col)

dom_idx = np.argmax(adata_spatial.obs[abund_cols].values, axis=1)
adata_spatial.obs["dominant_cell_type"] = [cell_types[i] for i in dom_idx]

# ── Save ──────────────────────────────────────────────────────────────────────
print("\nSaving results...")
adata_spatial.write_h5ad(os.path.join(RESULTS_DIR, "adata_gbm_mapped.h5ad"))

total_abund = adata_spatial.obs[abund_cols].sum().sort_values(ascending=False)
summary = pd.DataFrame({
    "Cell Type":        cell_types,
    "Total Abundance":  [adata_spatial.obs[f"{ct}_mean"].sum() for ct in cell_types],
    "Mean per Spot":    [adata_spatial.obs[f"{ct}_mean"].mean() for ct in cell_types],
    "Pct Spots Pos":    [(adata_spatial.obs[f"{ct}_mean"] > 0).mean() * 100 for ct in cell_types],
}).sort_values("Total Abundance", ascending=False)
summary.to_csv(os.path.join(RESULTS_DIR, "cell_abundance_summary.csv"), index=False)

# ── Figures ───────────────────────────────────────────────────────────────────
print("Generating figures...")

# Composition bar chart
fig, ax = plt.subplots(figsize=(14, 8))
ax.bar(range(len(total_abund)), total_abund.values, color=plt.cm.viridis(np.linspace(0, 1, len(total_abund))))
ax.set_xticks(range(len(total_abund)))
ax.set_xticklabels(total_abund.index, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Total Estimated Abundance")
ax.set_title("GBM Cell Type Composition")
ax.set_yscale("log")
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "cell_type_composition.png"), dpi=150, bbox_inches="tight")
plt.close()

# Dominant cell type map
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.spatial(adata_spatial, color="dominant_cell_type", ax=ax, size=1.2,
              title="Dominant Cell Type per Spot", legend_loc="right margin", show=False)
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "dominant_cell_type_map.png"), dpi=150, bbox_inches="tight")
plt.close()

# Six key cell types (2x3 grid)
highlight = ["MES-like", "Astrocyte", "CD4/CD8", "OPC", "AC-like", "Pericyte"]
available = [ct for ct in highlight if f"{ct}_mean" in adata_spatial.obs.columns]
if available:
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    for i, ct in enumerate(available):
        sc.pl.spatial(adata_spatial, color=f"{ct}_mean", cmap="winter", ax=axes[i],
                      title=ct, size=0.8, show=False)
        pct = total_abund.get(f"{ct}_mean", 0) / total_abund.sum() * 100
        axes[i].text(0.02, 0.98, f"{pct:.1f}%", transform=axes[i].transAxes,
                     fontsize=8, va="top", bbox=dict(boxstyle="round", fc="white", alpha=0.7))
    for j in range(len(available), len(axes)):
        axes[j].axis("off")
    plt.suptitle("Key GBM Cell Types", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "six_key_cell_types.png"), dpi=150, bbox_inches="tight")
    plt.close()

# All 14 cell types
n_cols = 4
n_rows = int(np.ceil(len(cell_types) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
axes = axes.flatten()
for i, ct in enumerate(cell_types):
    sc.pl.spatial(adata_spatial, color=f"{ct}_mean", cmap="winter", ax=axes[i],
                  title=ct, size=0.8, show=False)
for j in range(len(cell_types), len(axes)):
    axes[j].axis("off")
plt.suptitle("All GBM Cell Types — Spatial Distribution", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "all_cell_types_grid.png"), dpi=150, bbox_inches="tight")
plt.close()

print("\n" + "=" * 60)
print("GBM MAPPING COMPLETE")
print("=" * 60)
print(f"  Output directory: {RESULTS_DIR}")
for fn in sorted(os.listdir(RESULTS_DIR)):
    sz = os.path.getsize(os.path.join(RESULTS_DIR, fn)) / 1024
    print(f"    {fn}  ({sz:.0f} KB)")
print("=" * 60)
