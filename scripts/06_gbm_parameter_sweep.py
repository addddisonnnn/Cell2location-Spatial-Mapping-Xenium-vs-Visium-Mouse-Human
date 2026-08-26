"""
06_gbm_parameter_sweep.py
--------------------------
Runs the human GBM cell2location pipeline across five parameter configurations
to assess sensitivity of spatial mapping results to training duration,
learning rate, and posterior sampling.

Each configuration runs the full pipeline (regression + spatial model) and
saves results, figures, and a summary CSV to a separate subdirectory under
results/human_gbm/parameter_runs/.

This script wraps 05_gbm_spatial_mapping.py logic into a loop.
All five configurations share the same reference and Visium data loaded once
at the start to avoid repeated disk I/O.

Configurations (ordered fastest to slowest):
    run1_test          — 100 reg / 500 spat / n_comb=20 / lr=0.005  (~25 min)
    run2_more_samples  — same + 500 posterior samples                (~30 min)
    run3_higher_lr     — 100 reg / 500 spat / n_comb=20 / lr=0.025  (~25 min)
    run4_intermediate  — 500 reg / 5000 spat / n_comb=30             (~3.5 hr)
    run5_final         — 2000 reg / 10000 spat / n_comb=30           (~17 hr)

IMPORTANT: Must be run inside the Singularity container with GPU (V100).

Inputs:
    data/GBmap_reference/GBmap_reference_compressed.h5ad
    data/humanGlioblastoma/Parent_Visium_Human_Glioblastoma_*.h5

Outputs (one subdirectory per run):
    results/human_gbm/parameter_runs/{run_name}/
        inf_aver.csv
        adata_gbm_mapped.h5ad
        cell_abundance_summary.csv
        plots/
            composition.png
            dominant_cell_type.png
            six_key_cell_types.png
            all_cell_types_grid.png
            summary_table.png
    results/human_gbm/parameter_runs/
        comparison_across_runs.png

Usage:
    # Run all five configurations:
    python scripts/06_gbm_parameter_sweep.py

    # Run only specific configurations:
    python scripts/06_gbm_parameter_sweep.py --runs run1_test run5_final

    # Preview without running:
    python scripts/06_gbm_parameter_sweep.py --dry-run
"""

# GPU must be set before any theano import
import os
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32"

import sys
import argparse
import time
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import torch

from cell2location.models import RegressionGeneBackgroundCoverageTorch, LocationModelLinearDependentW

warnings.filterwarnings("ignore")

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="GBM parameter sweep across 5 configurations")
parser.add_argument("--data-dir", default=os.path.join(os.path.dirname(__file__), "..", "data"))
parser.add_argument(
    "--runs", nargs="+",
    default=["run1_test", "run2_more_samples", "run3_higher_lr", "run4_intermediate", "run5_final"],
    help="Which runs to execute"
)
parser.add_argument("--dry-run", action="store_true", help="Print plan without running")
parser.add_argument("--colormap", default="winter", help="Colormap for spatial plots (default: winter)")
args = parser.parse_args()

DATA_DIR     = os.path.abspath(args.data_dir)
REF_DIR      = os.path.join(DATA_DIR, "GBmap_reference")
GBM_DIR      = os.path.join(DATA_DIR, "humanGlioblastoma")
SWEEP_DIR    = os.path.join(DATA_DIR, "..", "results", "human_gbm", "parameter_runs")
os.makedirs(SWEEP_DIR, exist_ok=True)

use_cuda = torch.cuda.is_available()

# ── Configuration definitions ─────────────────────────────────────────────────
CONFIGS = {
    "run1_test": {
        "display": "Test Run (Baseline)",
        "n_iter_reg": 100,   "n_iter_spat": 500,   "n_comb": 20,
        "lr": 0.005, "n_samples": 200, "batch_size": 200,
        "est_min": 25,
    },
    "run2_more_samples": {
        "display": "More Posterior Samples",
        "n_iter_reg": 100,   "n_iter_spat": 500,   "n_comb": 20,
        "lr": 0.005, "n_samples": 500, "batch_size": 500,
        "est_min": 30,
    },
    "run3_higher_lr": {
        "display": "Higher Learning Rate",
        "n_iter_reg": 100,   "n_iter_spat": 500,   "n_comb": 20,
        "lr": 0.025, "n_samples": 200, "batch_size": 200,
        "est_min": 25,
    },
    "run4_intermediate": {
        "display": "Intermediate Run",
        "n_iter_reg": 500,   "n_iter_spat": 5000,  "n_comb": 30,
        "lr": 0.005, "n_samples": 500, "batch_size": 500,
        "est_min": 210,
    },
    "run5_final": {
        "display": "Final Run (Paper Quality)",
        "n_iter_reg": 2000,  "n_iter_spat": 10000, "n_comb": 30,
        "lr": 0.005, "n_samples": 500, "batch_size": 500,
        "est_min": 1020,
    },
}

# ── Dry run ───────────────────────────────────────────────────────────────────
selected = [r for r in args.runs if r in CONFIGS]
total_min = sum(CONFIGS[r]["est_min"] for r in selected)

print("=" * 70)
print("GBM PARAMETER SWEEP PLAN")
print("=" * 70)
print(f"{'Run':<22} {'Reg':>6} {'Spat':>7} {'n_comb':>7} {'LR':>7} {'Samples':>8} {'Est.':>8}")
print("-" * 70)
for r in selected:
    c = CONFIGS[r]
    print(f"{r:<22} {c['n_iter_reg']:>6} {c['n_iter_spat']:>7} {c['n_comb']:>7} "
          f"{c['lr']:>7.3f} {c['n_samples']:>8} {c['est_min']:>7}m")
print("-" * 70)
print(f"{'Total estimated time':>57}: {total_min}m ({total_min/60:.1f}h)")
print(f"GPU available: {use_cuda}")
print("=" * 70)

if args.dry_run:
    print("\nDry run — exiting without training.")
    sys.exit(0)

confirm = input("\nType 'yes' to start the sweep: ").strip().lower()
if confirm != "yes":
    print("Aborted.")
    sys.exit(0)

# ── Load data once ────────────────────────────────────────────────────────────
print("\nLoading shared data (done once for all runs)...")

compressed_ref = os.path.join(REF_DIR, "GBmap_reference_compressed.h5ad")
if not os.path.exists(compressed_ref):
    print(f"ERROR: {compressed_ref} not found. Run 02_prepare_gbmap_reference.py first.")
    sys.exit(1)

adata_ref_base = sc.read(compressed_ref)
print(f"  Reference: {adata_ref_base.n_obs:,} cells x {adata_ref_base.n_vars:,} genes")

adata_spatial_base = sc.read_visium(
    GBM_DIR,
    count_file="Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5",
    load_images=True,
)
adata_spatial_base.var_names_make_unique()
adata_ref_base.var_names_make_unique()
print(f"  Spatial:   {adata_spatial_base.n_obs:,} spots x {adata_spatial_base.n_vars:,} genes")

common_genes = np.intersect1d(adata_spatial_base.var_names, adata_ref_base.var_names)
print(f"  Shared genes: {len(common_genes):,}")

# ── Helper: run one configuration ─────────────────────────────────────────────
def run_one(config_key, config, adata_ref_base, adata_spatial_base, common_genes):
    run_dir = os.path.join(SWEEP_DIR, config_key)
    plots_dir = os.path.join(run_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print(f"\n{'='*70}")
    print(f"RUNNING: {config['display']}")
    print(f"{'='*70}")
    t0 = time.time()

    adata_ref     = adata_ref_base[:, common_genes].copy()
    adata_spatial = adata_spatial_base[:, common_genes].copy()
    adata_spatial.layers["counts"] = adata_spatial.X.copy()

    # Regression
    adata_ref.obs["dummy_covar"] = "1"
    cell2covar = pd.DataFrame({
        "patient":     adata_ref.obs["patient"].values,
        "cell_type":   adata_ref.obs["cell_type"].values,
        "dummy_covar": adata_ref.obs["dummy_covar"].values,
    }, index=adata_ref.obs_names)

    X_ref = adata_ref.X.toarray() if hasattr(adata_ref.X, "toarray") else np.array(adata_ref.X)
    X_ref = X_ref.astype(np.float32)

    mod_ref = RegressionGeneBackgroundCoverageTorch(
        sample_id="patient", cell2covar=cell2covar, X_data=X_ref,
        data_type="float32", n_iter=config["n_iter_reg"], learning_rate=config["lr"],
        total_grad_norm_constraint=200, verbose=False, use_cuda=use_cuda,
        minibatch_size=2500, use_average_as_initial_value=True,
    )
    mod_ref.fit_advi_iterative(n=1, n_type="restart", n_iter=config["n_iter_reg"], train_proportion=0.9)
    mod_ref.sample_posterior(node="all", n_samples=500, save_samples=False, mean_field_slot="init_1")

    inf_aver_raw = mod_ref.samples["post_sample_means"]["gene_factors"]
    fact_names   = mod_ref.fact_names
    ct_idx  = [i for i, f in enumerate(fact_names) if f.startswith("cell_type_")]
    ct_names = [f.replace("cell_type_", "") for f in fact_names if f.startswith("cell_type_")]
    inf_aver = pd.DataFrame(inf_aver_raw[ct_idx, :].T, index=adata_ref.var_names, columns=ct_names)
    inf_aver.to_csv(os.path.join(run_dir, "inf_aver.csv"))

    # Spatial model
    common_sp = np.intersect1d(adata_spatial.var_names, inf_aver.index)
    inf_aver_al = inf_aver.loc[common_sp, :].copy()
    X_sp = adata_spatial[:, common_sp].X.toarray() if hasattr(adata_spatial.X, "toarray") else np.array(adata_spatial[:, common_sp].X)
    X_sp = X_sp.astype(np.float32)

    mod_spatial = LocationModelLinearDependentW(
        cell_state_mat=inf_aver_al.values.astype(np.float32), X_data=X_sp,
        n_comb=config["n_comb"], data_type="float32", n_iter=config["n_iter_spat"],
        learning_rate=config["lr"], total_grad_norm_constraint=200, verbose=False,
        cell_number_prior={"cells_per_spot": 15, "factors_per_spot": 10, "combs_per_spot": 3},
    )
    mod_spatial.fit_advi_iterative()
    mod_spatial.sample_posterior(
        node="all", n_samples=config["n_samples"], save_samples=False,
        mean_field_slot="init_1", batch_size=config["batch_size"],
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

    # Save h5ad and summary
    adata_spatial.write_h5ad(os.path.join(run_dir, "adata_gbm_mapped.h5ad"))
    total_abund = adata_spatial.obs[abund_cols].sum().sort_values(ascending=False)
    summary = pd.DataFrame({
        "Cell Type": ct_names,
        "Total Abundance": [adata_spatial.obs[f"{ct}_mean"].sum() for ct in ct_names],
        "Pct of Total":    [adata_spatial.obs[f"{ct}_mean"].sum() / adata_spatial.obs[abund_cols].sum().sum() * 100 for ct in ct_names],
    }).sort_values("Total Abundance", ascending=False)
    summary.to_csv(os.path.join(run_dir, "cell_abundance_summary.csv"), index=False)

    # --- Figures ---
    CMAP = args.colormap

    # Composition
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.bar(range(len(total_abund)), total_abund.values)
    ax.set_xticks(range(len(total_abund)))
    ax.set_xticklabels(total_abund.index, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Total Abundance")
    ax.set_title(f"{config['display']}: Cell Type Composition")
    ax.set_yscale("log")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "composition.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # Dominant cell type map
    fig, ax = plt.subplots(figsize=(11, 9))
    sc.pl.spatial(adata_spatial, color="dominant_cell_type", ax=ax, size=1.0,
                  title=f"{config['display']}: Dominant Cell Type", legend_loc="right margin", show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "dominant_cell_type.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # Six key cell types
    highlight = ["MES-like", "Astrocyte", "CD4/CD8", "OPC", "AC-like", "Pericyte"]
    avail = [ct for ct in highlight if f"{ct}_mean" in adata_spatial.obs.columns]
    if avail:
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        for i, ct in enumerate(avail):
            sc.pl.spatial(adata_spatial, color=f"{ct}_mean", cmap=CMAP, ax=axes[i],
                          title=ct, size=0.8, show=False)
            pct = total_abund.get(f"{ct}_mean", 0) / total_abund.sum() * 100
            axes[i].text(0.02, 0.98, f"{pct:.1f}%", transform=axes[i].transAxes,
                         fontsize=8, va="top", bbox=dict(boxstyle="round", fc="white", alpha=0.7))
        for j in range(len(avail), len(axes)):
            axes[j].axis("off")
        plt.suptitle(f"{config['display']}: Key GBM Cell Types", fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "six_key_cell_types.png"), dpi=150, bbox_inches="tight")
        plt.close()

    # All cell types
    n_cols = 4
    n_rows = int(np.ceil(len(cell_types) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
    axes = axes.flatten()
    for i, ct in enumerate(cell_types):
        sc.pl.spatial(adata_spatial, color=f"{ct}_mean", cmap=CMAP, ax=axes[i],
                      title=ct, size=0.8, show=False)
    for j in range(len(cell_types), len(axes)):
        axes[j].axis("off")
    plt.suptitle(f"{config['display']}: All Cell Types", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "all_cell_types_grid.png"), dpi=150, bbox_inches="tight")
    plt.close()

    elapsed = (time.time() - t0) / 60
    print(f"  Completed in {elapsed:.1f} min | Top 3: {summary.head(3)['Cell Type'].tolist()}")
    return summary


# ── Execute sweep ─────────────────────────────────────────────────────────────
all_summaries = {}
for run_key in selected:
    try:
        summary = run_one(CONFIGS[run_key], CONFIGS[run_key], adata_ref_base, adata_spatial_base, common_genes)
        all_summaries[run_key] = summary
    except Exception as e:
        print(f"  ERROR in {run_key}: {e}")
        continue

# ── Cross-run comparison plot ─────────────────────────────────────────────────
if len(all_summaries) > 1:
    print("\nGenerating cross-run comparison plot...")
    fig, axes = plt.subplots(len(all_summaries), 1, figsize=(14, 5 * len(all_summaries)))
    if len(all_summaries) == 1:
        axes = [axes]
    for ax, (run_key, summary) in zip(axes, all_summaries.items()):
        top10 = summary.head(10)
        ax.bar(range(len(top10)), top10["Pct of Total"].values)
        ax.set_xticks(range(len(top10)))
        ax.set_xticklabels(top10["Cell Type"], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("% of Total Abundance")
        ax.set_title(f"{CONFIGS[run_key]['display']}")
    plt.suptitle("Parameter Sweep: Cell Type Composition Across Runs", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(SWEEP_DIR, "comparison_across_runs.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.join(SWEEP_DIR, 'comparison_across_runs.png')}")

print("\n" + "=" * 70)
print("PARAMETER SWEEP COMPLETE")
print("=" * 70)
print(f"  Results in: {SWEEP_DIR}")
for run_key in selected:
    status = "done" if run_key in all_summaries else "FAILED"
    print(f"  {run_key}: {status}")
print("=" * 70)
