"""
utils.py
--------
Shared utility functions used across the cell2location pipeline scripts.

Includes helpers for:
    - Verifying file existence before pipeline steps
    - Extracting cell type signatures from trained regression models
    - Adding spot factor abundances to AnnData objects
    - Generating standard spatial plot panels
    - Printing formatted progress/summary messages

Import this in any script:
    from scripts.utils import (
        check_files, extract_signatures, add_spot_factors,
        plot_spatial_panel, print_summary
    )
"""

import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


# ── File utilities ────────────────────────────────────────────────────────────

def check_files(file_dict, context=""):
    """
    Verify that all required files exist before starting a pipeline step.

    Parameters
    ----------
    file_dict : dict
        {label: filepath} pairs to check
    context : str
        Message to display in error output (e.g. 'Run script 02 first.')

    Returns
    -------
    bool
        True if all files exist, exits with error if any are missing.
    """
    missing = []
    for label, path in file_dict.items():
        if not os.path.exists(path):
            missing.append(f"  MISSING [{label}]: {path}")
    if missing:
        print(f"\nERROR: Required files not found:")
        for m in missing:
            print(m)
        if context:
            print(f"\n{context}")
        sys.exit(1)
    return True


def ensure_dir(*paths):
    """Create directories if they do not exist."""
    for p in paths:
        os.makedirs(p, exist_ok=True)


# ── Cell type signature extraction ────────────────────────────────────────────

def extract_signatures(mod, var_names, prefix="cell_type_"):
    """
    Extract cell type expression signatures from a trained
    RegressionGeneBackgroundCoverageTorch model.

    The regression model's gene_factors matrix contains columns for cell types,
    patient/sample effects, and any other covariates. This function filters to
    only the cell type columns (those starting with `prefix`).

    Parameters
    ----------
    mod : RegressionGeneBackgroundCoverageTorch
        Trained regression model with sampled posterior.
    var_names : list or Index
        Gene names (used as row index for the output DataFrame).
    prefix : str
        Column prefix that identifies cell type factors (default: 'cell_type_').

    Returns
    -------
    inf_aver : pd.DataFrame
        Shape (n_genes, n_cell_types). Values are posterior mean expression
        of each gene in each cell type.
    cell_type_names : list
        Cell type names (prefix stripped).
    """
    inf_aver_raw = mod.samples["post_sample_means"]["gene_factors"]
    fact_names   = mod.fact_names

    ct_indices   = [i for i, f in enumerate(fact_names) if f.startswith(prefix)]
    ct_names     = [f.replace(prefix, "") for f in fact_names if f.startswith(prefix)]

    if not ct_indices:
        raise ValueError(
            f"No factors found with prefix '{prefix}'. "
            f"Available factors: {fact_names[:10]}"
        )

    inf_aver = pd.DataFrame(
        inf_aver_raw[ct_indices, :].T,
        index=var_names,
        columns=ct_names,
    )
    return inf_aver, ct_names


def extract_signatures_xenium(mod):
    """
    Extract signatures from a Xenium-style regression model where
    cell type columns are identified by mod.which_sample mask.

    Parameters
    ----------
    mod : RegressionGeneBackgroundCoverageTorch

    Returns
    -------
    inf_aver : pd.DataFrame
    """
    inf_aver = mod.samples["post_sample_means"]["gene_factors"]
    inf_aver = pd.DataFrame(
        inf_aver,
        index=mod.fact_names,
        columns=mod.var_names,
    ).T
    # Keep only cell type columns (not sample/batch columns)
    inf_aver = inf_aver.loc[:, ~mod.which_sample.values]
    inf_aver.columns = [c.replace("annotation_1_", "") for c in inf_aver.columns]
    return inf_aver


# ── Spot factor utilities ──────────────────────────────────────────────────────

def add_spot_factors(adata, spot_factors, cell_type_names, col_suffix="_mean"):
    """
    Add spatial deconvolution results to an AnnData object.

    Adds one column per cell type to adata.obs (named '{ct}{col_suffix}')
    and a 'dominant_cell_type' column identifying the highest-abundance
    cell type at each spot.

    Parameters
    ----------
    adata : AnnData
    spot_factors : np.ndarray, shape (n_spots, n_cell_types)
    cell_type_names : list
    col_suffix : str

    Returns
    -------
    adata : AnnData (modified in place)
    abund_cols : list of column names added
    """
    abund_cols = []
    for i, ct in enumerate(cell_type_names):
        col = f"{ct}{col_suffix}"
        adata.obs[col] = spot_factors[:, i]
        abund_cols.append(col)

    dom_idx = np.argmax(adata.obs[abund_cols].values, axis=1)
    adata.obs["dominant_cell_type"] = [cell_type_names[i] for i in dom_idx]
    return adata, abund_cols


# ── Visualization ─────────────────────────────────────────────────────────────

def plot_spatial_panel(
    adata,
    cell_types,
    abund_col_suffix="_mean",
    n_cols=4,
    cmap="winter",
    title="",
    size=0.8,
    vmax_percentile=99,
    save_path=None,
    dpi=150,
    use_scanpy=True,
):
    """
    Plot a multi-panel grid of spatial cell type abundance maps.

    Works for Visium (using sc.pl.spatial) or Xenium (using scatter).
    Automatically uses scatter for Xenium data (detected by presence of
    x_centroid column instead of spatial image).

    Parameters
    ----------
    adata : AnnData
    cell_types : list of str
        Cell type names (without suffix).
    abund_col_suffix : str
    n_cols : int
    cmap : str
    title : str
    size : float
        Point size.
    vmax_percentile : int
    save_path : str or None
    dpi : int
    use_scanpy : bool
        If True use sc.pl.spatial (Visium); if False use scatter (Xenium).
    """
    n_rows = int(np.ceil(len(cell_types) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
    axes = axes.flatten()

    for i, ct in enumerate(cell_types):
        col = f"{ct}{abund_col_suffix}"
        if col not in adata.obs.columns:
            axes[i].set_visible(False)
            continue

        abundance = adata.obs[col].values
        vmax = np.percentile(abundance, vmax_percentile)

        if use_scanpy:
            sc.pl.spatial(
                adata, color=col, cmap=cmap, ax=axes[i], title=ct,
                size=size, vmin=0, vmax=vmax, show=False,
            )
        else:
            x = adata.obs["x_centroid"]
            y = adata.obs["y_centroid"]
            sc_plot = axes[i].scatter(x, y, c=abundance, cmap=cmap, s=size, vmin=0, vmax=vmax)
            axes[i].axis("off")
            plt.colorbar(sc_plot, ax=axes[i], shrink=0.7, label="Abundance")

        axes[i].set_title(ct, fontsize=10, fontweight="bold")

    for j in range(len(cell_types), len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(title, fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {save_path}")
    else:
        plt.show()


def plot_dominant_cell_type(adata, title="", save_path=None, dpi=150, use_scanpy=True):
    """
    Plot the dominant cell type at each spatial location.

    Parameters
    ----------
    adata : AnnData
        Must have 'dominant_cell_type' in obs.
    title : str
    save_path : str or None
    dpi : int
    use_scanpy : bool
    """
    if "dominant_cell_type" not in adata.obs.columns:
        print("WARNING: 'dominant_cell_type' not found in adata.obs. Run add_spot_factors first.")
        return

    fig, ax = plt.subplots(figsize=(12, 10))
    if use_scanpy:
        sc.pl.spatial(
            adata, color="dominant_cell_type", ax=ax, title=title or "Dominant Cell Type",
            size=1.0, legend_loc="right margin", show=False,
        )
    else:
        cell_types = adata.obs["dominant_cell_type"].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))
        ct_to_color = dict(zip(cell_types, colors))
        x = adata.obs["x_centroid"]
        y = adata.obs["y_centroid"]
        for ct in cell_types:
            mask = adata.obs["dominant_cell_type"] == ct
            ax.scatter(x[mask], y[mask], c=[ct_to_color[ct]], s=0.5, label=ct)
        ax.legend(loc="right", markerscale=5, fontsize=7, bbox_to_anchor=(1.15, 0.5))
        ax.set_title(title or "Dominant Cell Type")
        ax.axis("off")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {save_path}")
    else:
        plt.show()


# ── Formatting / progress ──────────────────────────────────────────────────────

class Timer:
    """Simple wall-clock timer for logging step durations."""

    def __init__(self, label=""):
        self.label = label
        self.start = time.time()

    def elapsed(self):
        return time.time() - self.start

    def log(self, step=""):
        t = self.elapsed()
        tag = f"[{self.label}] " if self.label else ""
        step_tag = f"{step}: " if step else ""
        print(f"  {tag}{step_tag}{t:.1f}s elapsed")


def print_section(title, width=60):
    print("\n" + "=" * width)
    print(title)
    print("=" * width)


def print_summary(items, title="SUMMARY", width=60):
    print("\n" + "=" * width)
    print(title)
    print("=" * width)
    for key, val in items.items():
        print(f"  {key:<30} {val}")
    print("=" * width)
