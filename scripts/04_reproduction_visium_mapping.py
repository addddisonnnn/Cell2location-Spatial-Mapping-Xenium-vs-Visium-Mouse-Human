"""
04_reproduction_visium_mapping.py
----------------------------------
Loads the mouse brain Visium data (E-MTAB-11114), trains the cell2location
regression model on the snRNA-seq reference, runs spatial mapping on section
ST8059048, and generates Figures 2c, 2d, and 2e from the original paper.

IMPORTANT: Run 03_reproduction_scrna_umap.py first to generate the
           processed snRNA-seq AnnData used as input here.

IMPORTANT: Set THEANO_FLAGS before running (handled by 00_gpu_setup.py import
           at the top of this file). Must be run inside the Singularity
           container with --nv GPU access.

Inputs:
    results/reproduction/adata_snrna_processed.h5ad  (from script 03)
    data/mouseVisium/
        ST8059048_filtered_feature_bc_matrix.h5
        ST8059048_spatial/spatial/

Outputs:
    results/reproduction/
        inf_aver.csv             — reference signature matrix (genes x cell types)
        adata_vis_mapped.h5ad    — Visium AnnData with cell type abundances added
        figure_2c_HE.png         — H&E tissue image (Figure 2c)
        figure_2d_major.png      — major regional cell types (Figure 2d)
        figure_2e_sparse.png     — sparse inhibitory neurons (Figure 2e)

Key parameters (match paper Table 1):
    cells_per_spot  = 8
    alpha_mean      = 200   (detection_alpha)
    n_epochs        = 30000 (set N_ITER below; use 3000 for a quick test)
    learning_rate   = 0.001

Usage:
    python scripts/04_reproduction_visium_mapping.py
    python scripts/04_reproduction_visium_mapping.py --n-iter 3000  # quick test
"""

# GPU MUST be configured before any theano import
import os
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32"

import sys
import argparse
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import cell2location

warnings.filterwarnings("ignore")

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Run cell2location on mouse brain Visium")
parser.add_argument(
    "--data-dir",
    default=os.path.join(os.path.dirname(__file__), "..", "data"),
    help="Base data directory (default: data/)"
)
parser.add_argument(
    "--n-iter-reg", type=int, default=100,
    help="Regression model training epochs (default: 100; paper used 250)"
)
parser.add_argument(
    "--n-iter-spatial", type=int, default=30000,
    help="Spatial model training epochs (default: 30000; use 3000 for quick test)"
)
parser.add_argument(
    "--cells-per-spot", type=int, default=8,
    help="Expected cells per Visium spot prior (default: 8, matches paper)"
)
parser.add_argument(
    "--section", default="ST8059048",
    help="Which Visium section to map (default: ST8059048)"
)
args = parser.parse_args()

DATA_DIR    = os.path.abspath(args.data_dir)
VISIUM_DIR  = os.path.join(DATA_DIR, "mouseVisium")
RESULTS_DIR = os.path.join(DATA_DIR, "..", "results", "reproduction")
os.makedirs(RESULTS_DIR, exist_ok=True)

SCRNA_H5AD     = os.path.join(RESULTS_DIR, "adata_snrna_processed.h5ad")
INF_AVER_CSV   = os.path.join(RESULTS_DIR, "inf_aver.csv")
VIS_H5AD       = os.path.join(RESULTS_DIR, "adata_vis_mapped.h5ad")
FIG_2C         = os.path.join(RESULTS_DIR, "figure_2c_HE.png")
FIG_2D         = os.path.join(RESULTS_DIR, "figure_2d_major.png")
FIG_2E         = os.path.join(RESULTS_DIR, "figure_2e_sparse.png")

# ── Load snRNA-seq reference ──────────────────────────────────────────────────
print("=" * 60)
print("LOADING snRNA-seq REFERENCE")
print("=" * 60)

if not os.path.exists(SCRNA_H5AD):
    print(f"ERROR: {SCRNA_H5AD} not found.")
    print("Run 03_reproduction_scrna_umap.py first.")
    sys.exit(1)

adata_ref = sc.read_h5ad(SCRNA_H5AD)
print(f"  Loaded: {adata_ref.n_obs:,} cells x {adata_ref.n_vars:,} genes")
print(f"  Cell types: {adata_ref.obs['annotation_1'].nunique()}")

# Gene filtering — keep genes expressed in >5% of cells
n_cells = adata_ref.n_obs
gene_counts = (adata_ref.X > 0).sum(axis=0)
if hasattr(gene_counts, "A1"):
    gene_counts = gene_counts.A1
pct_expressed = gene_counts / n_cells
adata_ref = adata_ref[:, pct_expressed > 0.05].copy()
print(f"  After >5% expression filter: {adata_ref.n_vars:,} genes")

# ── Regression model — reference signatures ───────────────────────────────────
print(f"\n{'='*60}")
print(f"REGRESSION MODEL (n_iter={args.n_iter_reg})")
print("=" * 60)

# Required by cell2location v0.3 API: covariate_col_names must be a real column
adata_ref.obs["dummy_covar"] = "1"

results_regression = cell2location.run_regression(
    sc_data=adata_ref,
    train_args={
        "covariate_col_names": ["dummy_covar"],
        "sample_name_col": "sample",
        "tech_name_col": None,
        "n_epochs": args.n_iter_reg,
        "minibatch_size": 2500,
        "learning_rate": 0.01,
        "use_average_as_initial_value": True,
        "use_cuda": True,
        "train_proportion": 0.9,
        "l2_weight": True,
        "use_raw": False,
    },
    export_args={"path": os.path.join(RESULTS_DIR, "regression_output/"), "save_model": False},
)

# Extract reference signature matrix
inf_aver_raw = results_regression["mod"].samples["post_sample_means"]["per_cluster_mu_fg"]
inf_aver = pd.DataFrame(
    inf_aver_raw,
    index=adata_ref.var_names,
    columns=adata_ref.uns["mod_cell_types"],
)
inf_aver.to_csv(INF_AVER_CSV)
print(f"  Reference signatures shape: {inf_aver.shape}  (genes x cell types)")
print(f"  Saved: {INF_AVER_CSV}")

# ── Load Visium spatial data ──────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"LOADING VISIUM SECTION {args.section}")
print("=" * 60)

h5_path      = os.path.join(VISIUM_DIR, f"{args.section}_filtered_feature_bc_matrix.h5")
spatial_path = os.path.join(VISIUM_DIR, f"{args.section}_spatial")

if not os.path.exists(h5_path):
    print(f"ERROR: {h5_path} not found.")
    sys.exit(1)

adata_vis = sc.read_visium(
    path=spatial_path,
    count_file=h5_path,
    load_images=True,
)
adata_vis.var_names_make_unique()
adata_vis.obs["sample"] = args.section
print(f"  Loaded: {adata_vis.n_obs:,} spots x {adata_vis.n_vars:,} genes")

# Figure 2c — H&E image
print(f"\nGenerating Figure 2c (H&E)...")
fig, ax = plt.subplots(figsize=(10, 10))
sc.pl.spatial(adata_vis, img_key="hires", color=None, title=f"H&E — {args.section}", ax=ax, show=False)
plt.tight_layout()
plt.savefig(FIG_2C, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIG_2C}")

# ── Gene intersection ─────────────────────────────────────────────────────────
print("\nIntersecting genes between Visium and reference...")
# Preserve raw counts then align
adata_vis.layers["counts"] = adata_vis.X.copy()
adata_vis.raw = adata_vis

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis  = adata_vis[:, intersect].copy()
inf_aver   = inf_aver.loc[intersect, :].copy()
print(f"  Shared genes: {len(intersect):,}")

# ── Spatial model ─────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"SPATIAL MAPPING (n_iter={args.n_iter_spatial})")
print("=" * 60)

results_c2l = cell2location.run_cell2location(
    sc_ref=inf_aver,
    sp_data=adata_vis,
    train_args={
        "n_epochs": args.n_iter_spatial,
        "minibatch_size": None,
        "learning_rate": 0.001,
        "use_cuda": True,
    },
    model_kwargs={
        "cell_number_prior": {
            "cells_per_spot": args.cells_per_spot,
            "factors_per_spot": 7,
            "combs_per_spot": 2.5,
        },
        "cell_number_var": {
            "cells_per_spot_dtype": "float32",
            "alpha_polynomial": 1,
            "alpha_mean": 200,
        },
    },
    export_args={"path": os.path.join(RESULTS_DIR, "spatial_output/"), "save_model": False},
)

adata_vis = results_c2l["sp_data"]
print(f"  Mapping complete. Obs columns: {len(adata_vis.obs.columns)}")

# ── Save mapped AnnData ───────────────────────────────────────────────────────
adata_vis.write_h5ad(VIS_H5AD)
print(f"  Saved: {VIS_H5AD}")

# ── Figure 2d — major regional subtypes ──────────────────────────────────────
print("\nGenerating Figure 2d (major regional subtypes)...")
fig2d_types = [ct for ct in ["Oligo_2", "Inh_Meis2_3", "Inh_4", "Ext_Thal_1", "Ext_L23", "Ext_L56"]
               if ct in adata_vis.obs.columns]

n_cols = 3
n_rows = int(np.ceil(len(fig2d_types) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
axes = axes.flatten()
for i, ct in enumerate(fig2d_types):
    sc.pl.spatial(adata_vis, color=ct, img_key="hires", cmap="magma",
                  size=1.3, vmin=0, vmax="p99", ax=axes[i], title=ct, show=False)
for j in range(len(fig2d_types), len(axes)):
    axes[j].axis("off")
plt.suptitle(f"Figure 2d — Major Regional Cell Types ({args.n_iter_spatial} iter)", fontsize=13)
plt.tight_layout()
plt.savefig(FIG_2D, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIG_2D}")

# ── Figure 2e — sparse inhibitory neurons ────────────────────────────────────
print("Generating Figure 2e (sparse inhibitory neurons)...")
fig2e_types = [ct for ct in ["Inh_Sst", "Inh_Lamp5", "Inh_Vip"]
               if ct in adata_vis.obs.columns]

fig, axes = plt.subplots(1, len(fig2e_types), figsize=(6 * len(fig2e_types), 6))
if len(fig2e_types) == 1:
    axes = [axes]
for i, ct in enumerate(fig2e_types):
    sc.pl.spatial(adata_vis, color=ct, img_key="hires", cmap="magma",
                  size=1.3, vmin=0, vmax="p99", ax=axes[i], title=ct, show=False)
plt.suptitle(f"Figure 2e — Sparse Inhibitory Neurons ({args.n_iter_spatial} iter)", fontsize=13)
plt.tight_layout()
plt.savefig(FIG_2E, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIG_2E}")

print("\n" + "=" * 60)
print("REPRODUCTION COMPLETE")
print("=" * 60)
print(f"  Figure 2b:  results/reproduction/figure_2b_umap.png  (from script 03)")
print(f"  Figure 2c:  {FIG_2C}")
print(f"  Figure 2d:  {FIG_2D}")
print(f"  Figure 2e:  {FIG_2E}")
print(f"  Mapped h5ad:{VIS_H5AD}")
print("=" * 60)
