"""
07_xenium_spatial_mapping.py
-----------------------------
Applies cell2location to 10x Xenium mouse brain data at near single-cell
resolution (~10 µm per location). Uses the same mouse brain snRNA-seq reference
atlas from the original paper's tutorial.

Key difference from Visium: each Xenium location contains ~1 cell (not 5-15),
so N_cells_per_location and other priors must be adjusted accordingly. The gene
panel is also targeted (248 genes) rather than whole-transcriptome, which reduces
the number of shared genes between reference and spatial data.

This script runs one configuration. Adjust the SPEED SETTINGS block for the
desired quality/runtime tradeoff.

IMPORTANT: Must be run inside the Singularity container with GPU (V100).

Inputs (auto-downloaded if not present):
    results/xenium/all_cells_20200625.h5ad             — reference atlas (~700MB)
    results/xenium/snRNA_annotation_astro_subtypes_refined59_20200823.csv

Inputs (must be downloaded manually):
    data/xeniumMouseBrain/cell_feature_matrix.h5
    data/xeniumMouseBrain/cells.csv.gz

Outputs:
    results/xenium/
        adata_xenium_mapped.h5ad     — Xenium AnnData with cell type abundances
        inf_aver.csv                 — reference signature matrix (247 x 59)
        cell_type_abundances_8ct.png — spatial plots of 8 selected cell types
        all_cell_types_59ct.png      — spatial plots of all 59 cell types

Usage:
    # Quick test (~20 min on V100):
    python scripts/07_xenium_spatial_mapping.py --n-iter-reg 100 --n-iter-spat 500 --n-comb 20

    # Final run (~3 hr on V100):
    python scripts/07_xenium_spatial_mapping.py --n-iter-reg 4000 --n-iter-spat 20000 --n-comb 50

    # Full parameter sweep (all 5 runs):
    python scripts/07_xenium_spatial_mapping.py --sweep
"""

# GPU must be set before any theano import
import os
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32"

import sys
import argparse
import warnings
import urllib.request

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import torch

from cell2location.models import (
    RegressionGeneBackgroundCoverageTorch,
    LocationModelLinearDependentWMultiExperiment,
)

warnings.filterwarnings("ignore")

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Run cell2location on Xenium mouse brain")
parser.add_argument("--data-dir", default=os.path.join(os.path.dirname(__file__), "..", "data"))
parser.add_argument("--n-iter-reg",  type=int, default=100,   help="Regression iterations (default 100)")
parser.add_argument("--n-iter-spat", type=int, default=500,   help="Spatial iterations (default 500)")
parser.add_argument("--n-comb",      type=int, default=20,    help="n_comb parameter (default 20)")
parser.add_argument("--lr",          type=float, default=0.005)
parser.add_argument("--n-samples",   type=int, default=200)
parser.add_argument("--batch-size",  type=int, default=200)
parser.add_argument("--train-split", type=float, default=0.9, help="Train/test split (default 0.9)")
parser.add_argument("--run-name",    default=None,            help="Output subdirectory name")
parser.add_argument(
    "--sweep", action="store_true",
    help="Run all 5 parameter configurations sequentially"
)
args = parser.parse_args()

DATA_DIR    = os.path.abspath(args.data_dir)
XENIUM_DIR  = os.path.join(DATA_DIR, "xeniumMouseBrain")
RESULTS_DIR = os.path.join(DATA_DIR, "..", "results", "xenium")
os.makedirs(RESULTS_DIR, exist_ok=True)

use_cuda = torch.cuda.is_available()
print(f"GPU available: {use_cuda}")
if use_cuda:
    print(f"  GPU: {torch.cuda.get_device_name(0)}")

# ── Reference URLs ────────────────────────────────────────────────────────────
REF_H5AD_URL = "https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_snrna/all_cells_20200625.h5ad"
REF_CSV_URL  = (
    "https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_snrna/"
    "snRNA_annotation_astro_subtypes_refined59_20200823.csv"
)
REF_H5AD_PATH = os.path.join(RESULTS_DIR, "all_cells_20200625.h5ad")
REF_CSV_PATH  = os.path.join(RESULTS_DIR, "snRNA_annotation_astro_subtypes_refined59_20200823.csv")

# ── Step 1: Download reference if needed ──────────────────────────────────────
def download_if_missing(url, path, label):
    if os.path.exists(path):
        print(f"  {label}: already present ({os.path.getsize(path)/1024/1024:.1f} MB)")
        return
    print(f"  Downloading {label} from {url}...")
    urllib.request.urlretrieve(url, path)
    print(f"  Done: {os.path.getsize(path)/1024/1024:.1f} MB")

print("=" * 60)
print("DOWNLOADING REFERENCE (if not present)")
print("=" * 60)
download_if_missing(REF_H5AD_URL, REF_H5AD_PATH, "reference h5ad (~700MB)")
download_if_missing(REF_CSV_URL,  REF_CSV_PATH,  "annotation CSV")

# ── Step 2: Load and prepare reference ────────────────────────────────────────
print("\n" + "=" * 60)
print("LOADING AND PREPARING REFERENCE")
print("=" * 60)

adata_ref = sc.read_h5ad(REF_H5AD_PATH)
ann = pd.read_csv(REF_CSV_PATH, index_col=0)
adata_ref.obs = adata_ref.obs.join(ann, how="left")

# Swap Ensembl IDs to gene symbols for matching with Xenium gene names
adata_ref.var["SYMBOL"] = adata_ref.var["SYMBOL"].astype(str)
adata_ref.var_names = adata_ref.var["SYMBOL"]
adata_ref.var_names_make_unique()
print(f"  Reference: {adata_ref.n_obs:,} cells x {adata_ref.n_vars:,} genes")
print(f"  Cell types: {adata_ref.obs['annotation_1'].nunique()}")
print(f"  Reference gene examples: {adata_ref.var_names[:5].tolist()}")

# ── Step 3: Load Xenium spatial data ──────────────────────────────────────────
print("\n" + "=" * 60)
print("LOADING XENIUM SPATIAL DATA")
print("=" * 60)

h5_path  = os.path.join(XENIUM_DIR, "cell_feature_matrix.h5")
csv_path = os.path.join(XENIUM_DIR, "cells.csv.gz")
for p in [h5_path, csv_path]:
    if not os.path.exists(p):
        print(f"ERROR: Missing: {p}")
        print("Download the Xenium dataset first. See README.md Part 3.")
        sys.exit(1)

adata_xenium = sc.read_10x_h5(h5_path)
adata_xenium.var_names_make_unique()

cells = pd.read_csv(csv_path).set_index("cell_id")
cells.index = cells.index.astype(str)
adata_xenium.obs = adata_xenium.obs.join(cells[["x_centroid", "y_centroid"]], how="left")
adata_xenium.obsm["spatial"] = adata_xenium.obs[["x_centroid", "y_centroid"]].values
print(f"  Xenium: {adata_xenium.n_obs:,} cells x {adata_xenium.n_vars:,} genes")
print(f"  Gene examples: {adata_xenium.var_names[:5].tolist()}")

# ── Step 4: Gene intersection and QC ─────────────────────────────────────────
print("\n" + "=" * 60)
print("GENE INTERSECTION AND QC")
print("=" * 60)

shared_genes = adata_xenium.var_names.intersection(adata_ref.var_names)
print(f"  Shared genes: {len(shared_genes):,}")
print(f"  Xenium genes NOT in reference: {set(adata_xenium.var_names) - set(adata_ref.var_names)}")

adata_xenium = adata_xenium[:, shared_genes].copy()
adata_ref    = adata_ref[:, shared_genes].copy()

n_before = adata_xenium.n_obs
sc.pp.filter_cells(adata_xenium, min_counts=10)
sc.pp.filter_genes(adata_xenium, min_cells=5)
print(f"  Cells removed by QC: {n_before - adata_xenium.n_obs:,}")
print(f"  After QC: {adata_xenium.n_obs:,} cells x {adata_xenium.n_vars:,} genes")


# ── Core pipeline function ────────────────────────────────────────────────────
def run_xenium_pipeline(run_name, cfg, adata_ref, adata_xenium):
    run_dir = os.path.join(RESULTS_DIR, run_name)
    os.makedirs(run_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"RUN: {run_name}")
    print(f"  n_iter_reg={cfg['n_iter_reg']}, n_iter_spat={cfg['n_iter_spat']}, "
          f"n_comb={cfg['n_comb']}, lr={cfg['lr']}, train_split={cfg['train_split']}")
    print("=" * 60)

    # Regression model
    cell2covar = pd.DataFrame({
        "sample":       adata_ref.obs["sample"].values,
        "annotation_1": adata_ref.obs["annotation_1"].values,
    }, index=adata_ref.obs_names)

    X_ref = adata_ref.X.toarray() if hasattr(adata_ref.X, "toarray") else np.array(adata_ref.X)

    mod = RegressionGeneBackgroundCoverageTorch(
        sample_id="sample",
        cell2covar=cell2covar,
        X_data=X_ref,
        n_iter=cfg["n_iter_reg"],
        learning_rate=cfg["lr"],
        use_cuda=use_cuda,
        var_names=adata_ref.var_names.tolist(),
        obs_names=adata_ref.obs_names.tolist(),
    )
    mod.fit_advi_iterative(n=1, n_iter=cfg["n_iter_reg"], train_proportion=cfg["train_split"])
    mod.sample_posterior(node="all", n_samples=1000, save_samples=False, mean_field_slot="init_1")

    inf_aver = mod.samples["post_sample_means"]["gene_factors"]
    inf_aver = pd.DataFrame(inf_aver, index=mod.fact_names, columns=mod.var_names).T
    inf_aver = inf_aver.loc[:, ~mod.which_sample.values]
    inf_aver.columns = [c.replace("annotation_1_", "") for c in inf_aver.columns]
    inf_aver.to_csv(os.path.join(run_dir, "inf_aver.csv"))
    print(f"  Signatures: {inf_aver.shape}")

    # Spatial model
    X_xenium = adata_xenium.X.toarray() if hasattr(adata_xenium.X, "toarray") else np.array(adata_xenium.X)

    mod_spatial = LocationModelLinearDependentWMultiExperiment(
        cell_state_mat=inf_aver.values,
        X_data=X_xenium,
        n_comb=cfg["n_comb"],
        n_iter=cfg["n_iter_spat"],
        learning_rate=cfg["lr"],
        sample_id=np.array(["xenium_sample"] * adata_xenium.n_obs),
        var_names=adata_xenium.var_names.tolist(),
        obs_names=adata_xenium.obs_names.tolist(),
        fact_names=inf_aver.columns.tolist(),
    )
    mod_spatial.fit_advi(n=1, n_type="restart")
    mod_spatial.sample_posterior(
        node="all", n_samples=cfg["n_samples"], save_samples=False,
        mean_field_slot="init_1", batch_size=cfg["batch_size"],
    )

    spot_factors = mod_spatial.samples["post_sample_means"]["spot_factors"]
    adata_out = adata_xenium.copy()
    adata_out.obsm["cell_abundance"] = spot_factors
    for i, ct in enumerate(inf_aver.columns):
        adata_out.obs[ct] = spot_factors[:, i]

    adata_out.write_h5ad(os.path.join(run_dir, "adata_xenium_mapped.h5ad"))
    print(f"  Saved: {os.path.join(run_dir, 'adata_xenium_mapped.h5ad')}")

    # Figures
    x = adata_out.obs["x_centroid"]
    y = adata_out.obs["y_centroid"]

    # 8 selected cell types
    selected_8 = ["Astro_CTX", "Oligo_2", "Ext_L23", "Ext_Hpc_CA1", "Micro", "Inh_Sst", "OPC_1", "Ext_L56"]
    available_8 = [ct for ct in selected_8 if ct in adata_out.obs.columns]
    n_cols, n_rows = 4, int(np.ceil(len(available_8) / 4))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
    axes = axes.flatten()
    for i, ct in enumerate(available_8):
        ab = adata_out.obs[ct]
        sc_plot = axes[i].scatter(x, y, c=ab, cmap="viridis", s=0.5, vmin=0, vmax=np.percentile(ab, 99))
        axes[i].set_title(ct, fontsize=11)
        axes[i].axis("off")
        plt.colorbar(sc_plot, ax=axes[i], shrink=0.7)
    for j in range(len(available_8), len(axes)):
        axes[j].set_visible(False)
    plt.suptitle(f"Xenium Mouse Brain — 8 Cell Types ({run_name})", fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(run_dir, "cell_type_abundances_8ct.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # All 59 cell types
    all_cts = inf_aver.columns.tolist()
    n_cols_all = 6
    n_rows_all = int(np.ceil(len(all_cts) / n_cols_all))
    fig, axes = plt.subplots(n_rows_all, n_cols_all, figsize=(4 * n_cols_all, 4 * n_rows_all))
    axes = axes.flatten()
    for i, ct in enumerate(all_cts):
        if ct not in adata_out.obs.columns:
            axes[i].set_visible(False)
            continue
        ab = adata_out.obs[ct]
        axes[i].scatter(x, y, c=ab, cmap="magma", s=0.3, vmin=0, vmax=np.percentile(ab, 99))
        axes[i].set_title(ct, fontsize=7)
        axes[i].axis("off")
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle(f"Xenium Mouse Brain — All 59 Cell Types ({run_name})", fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(run_dir, "all_cell_types_59ct.png"), dpi=120, bbox_inches="tight")
    plt.close()

    print(f"  Figures saved to {run_dir}/")
    return adata_out


# ── Configuration definitions ─────────────────────────────────────────────────
SWEEP_CONFIGS = {
    "run1": {"n_iter_reg": 100,  "n_iter_spat": 500,   "n_comb": 20, "lr": 0.005, "train_split": 0.9, "n_samples": 200, "batch_size": 200},
    "run2": {"n_iter_reg": 500,  "n_iter_spat": 2500,  "n_comb": 20, "lr": 0.005, "train_split": 0.9, "n_samples": 200, "batch_size": 200},
    "run3": {"n_iter_reg": 100,  "n_iter_spat": 500,   "n_comb": 20, "lr": 0.025, "train_split": 0.9, "n_samples": 200, "batch_size": 200},
    "run4": {"n_iter_reg": 100,  "n_iter_spat": 500,   "n_comb": 20, "lr": 0.005, "train_split": 0.7, "n_samples": 400, "batch_size": 400},
    "run5": {"n_iter_reg": 4000, "n_iter_spat": 20000, "n_comb": 50, "lr": 0.005, "train_split": 0.9, "n_samples": 200, "batch_size": 200},
}

# ── Execute ───────────────────────────────────────────────────────────────────
if args.sweep:
    for run_key, cfg in SWEEP_CONFIGS.items():
        try:
            run_xenium_pipeline(run_key, cfg, adata_ref.copy(), adata_xenium.copy())
        except Exception as e:
            print(f"  ERROR in {run_key}: {e}")
else:
    run_name = args.run_name or f"reg{args.n_iter_reg}_spat{args.n_iter_spat}_comb{args.n_comb}"
    cfg = {
        "n_iter_reg": args.n_iter_reg, "n_iter_spat": args.n_iter_spat,
        "n_comb": args.n_comb, "lr": args.lr,
        "train_split": args.train_split, "n_samples": args.n_samples,
        "batch_size": args.batch_size,
    }
    run_xenium_pipeline(run_name, cfg, adata_ref, adata_xenium)

print("\n" + "=" * 60)
print("XENIUM ANALYSIS COMPLETE")
print(f"  Results in: {RESULTS_DIR}")
print("=" * 60)
