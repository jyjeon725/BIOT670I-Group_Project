#!/usr/bin/env python3
"""
This script loads a preprocessed and clustered AnnData object (.h5ad),
optionally performs differential expression (DE), and generates a variety
of plots commonly used in single‑cell RNA‑seq analysis.

It is designed for PBMC 3k–style datasets but works for any clustered AnnData.
"""

import argparse
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Ensures a "figures" directory exists before saving plots
os.makedirs("figures", exist_ok=True)

# -------------------------------------------------------------------------
# Command‑line interface
# -------------------------------------------------------------------------
# Users can choose which analyses and plots to run.
parser = argparse.ArgumentParser(description="Run DE and selected plots on a preprocessed clustered AnnData.")

parser.add_argument("--input", "-i", required=True,
                    help="Path to preprocessed clustered AnnData (.h5ad).")

# Optional analysis steps
parser.add_argument("--run-de", action="store_true",
                    help="Run differential expression (rank_genes_groups).")

# Optional plotting steps
parser.add_argument("--umap", action="store_true", help="Generate UMAP colored by clusters and markers.")
parser.add_argument("--heatmap", action="store_true", help="Generate heatmap of marker genes across clusters.")
parser.add_argument("--matrixplot", action="store_true", help="Generate matrixplot of marker genes by cluster.")
parser.add_argument("--dotplot", action="store_true", help="Generate dotplot of marker genes by cluster.")
parser.add_argument("--tracksplot", action="store_true", help="Generate tracksplot of marker genes ordered by cluster.")
parser.add_argument("--stacked-violin", action="store_true", help="Generate stacked violin plot for marker genes.")
parser.add_argument("--violin", action="store_true", help="Generate multi-panel violin plots for marker genes.")
parser.add_argument("--volcano", action="store_true", help="Generate volcano plots for each cluster.")  # NEW

# Convenience flag to run everything
parser.add_argument("--all-plots", action="store_true", help="Enable all plotting options.")

# Whether to show plots interactively instead of saving them
parser.add_argument("--show", action="store_true", help="Show plots interactively instead of saving.")

# Number of DE genes to compute per cluster
parser.add_argument("--top-n-de", type=int, default=50,
                    help="Number of top DE genes per cluster to compute/save.")

args = parser.parse_args()

# If user requests all plots, enable all plot flags
if args.all_plots:
    args.umap = args.heatmap = args.matrixplot = args.dotplot = True
    args.tracksplot = args.stacked_violin = args.violin = True
    args.volcano = True


# -------------------------------------------------------------------------
# Interactive menu (only used if user provides no flags)
# -------------------------------------------------------------------------
def interactive_menu():
    """
    Provides a simple numbered menu so users can choose actions interactively.
    """
    print("\nNo action flags provided. Choose operations to run (enter numbers separated by commas):")
    menu = [
        "Run differential expression (rank_genes_groups)",
        "UMAP (clusters and markers)",
        "Heatmap (marker genes by cluster)",
        "Matrixplot (mean expression + fraction)",
        "Dotplot (compact overview)",
        "Tracksplot (per-cell ordered expression)",
        "Stacked violin (many genes compactly)",
        "Violin (multi-panel per-gene)",
        "Volcano plots (per-cluster DE visualization)",
        "Show plots interactively (toggle)"
    ]

    for i, item in enumerate(menu, 1):
        print(f"  {i}. {item}")

    choice = input("Selection (e.g., 1,2,5) or press Enter to run nothing: ").strip()
    if not choice:
        return

    sel = {int(x.strip()) for x in choice.split(",") if x.strip().isdigit()}

    args.run_de = 1 in sel
    args.umap = 2 in sel
    args.heatmap = 3 in sel
    args.matrixplot = 4 in sel
    args.dotplot = 5 in sel
    args.tracksplot = 6 in sel
    args.stacked_violin = 7 in sel
    args.violin = 8 in sel
    args.volcano = 9 in sel
    if 10 in sel:
        args.show = True


# If user provided no flags, open the menu
action_flags = any([
    args.run_de, args.umap, args.heatmap, args.matrixplot, args.dotplot,
    args.tracksplot, args.stacked_violin, args.violin, args.volcano
])
if not action_flags:
    interactive_menu()


# -------------------------------------------------------------------------
# Load the AnnData object
# -------------------------------------------------------------------------
print(f"Loading preprocessed AnnData from: {args.input}")
try:
    adata = sc.read(args.input)
except Exception as e:
    print(f"Error loading file: {e}", file=sys.stderr)
    sys.exit(1)

# Check that clustering labels exist
if "leiden" not in adata.obs and "louvain" not in adata.obs:
    print("Warning: no clustering column ('leiden' or 'louvain') found in adata.obs.", file=sys.stderr)

# Ensure raw counts exist for DE
if adata.raw is None:
    adata.raw = adata.copy()


# -------------------------------------------------------------------------
# Helper: run DE only when needed
# -------------------------------------------------------------------------
def compute_de_if_needed():
    """
    Runs differential expression if:
    - the user requested it, OR
    - volcano plots require it and DE has not been computed yet.
    """
    if "rank_genes_groups" in adata.uns:
        return  # Already computed

    print("Computing differential expression (required for volcano plots)...")
    groupby_key = "leiden" if "leiden" in adata.obs else "louvain"

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby_key,
        method="t-test",
        n_genes=args.top_n_de,
        use_raw=True,
    )

    print("DE computation complete.")


# -------------------------------------------------------------------------
# Run DE if requested
# -------------------------------------------------------------------------
if args.run_de:
    compute_de_if_needed()

    # Convert Scanpy DE results into a tidy DataFrame
    def rank_genes_to_df(adata):
        result = adata.uns["rank_genes_groups"]
        groups = result["names"].dtype.names
        dfs = []
        for g in groups:
            df = pd.DataFrame({
                "group": g,
                "gene": result["names"][g],
                "logfoldchange": result.get("logfoldchanges", {g: [None]})[g],
                "pvals_adj": result.get("pvals_adj", {g: [None]})[g],
                "scores": result.get("scores", {g: [None]})[g],
            })
            dfs.append(df)
        return pd.concat(dfs, ignore_index=True)

    de_df = rank_genes_to_df(adata)
    de_df.to_csv("DE_results_by_cluster.csv", index=False)
    print("Saved DE results to DE_results_by_cluster.csv")


# -------------------------------------------------------------------------
# Volcano plot helper function
# -------------------------------------------------------------------------
def plot_volcano(de_df, group, save_fig=True, show=False):
    """
    Creates a volcano plot for a single cluster.
    X-axis: log2 fold change
    Y-axis: -log10 adjusted p-value
    """
    df = de_df[de_df["group"] == group].copy()

    # Clean up infinite or missing values
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchange", "pvals_adj"])

    # Compute -log10(p)
    df["neglog10_p"] = -np.log10(df["pvals_adj"] + 1e-300)

    plt.figure(figsize=(7, 6))
    plt.scatter(
        df["logfoldchange"],
        df["neglog10_p"],
        c=(df["pvals_adj"] < 0.05),  # Color significant genes differently
        cmap="coolwarm",
        alpha=0.7,
        edgecolor="none"
    )

    plt.axvline(0, color="black", linestyle="--", linewidth=1)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 Adjusted p-value")
    plt.title(f"Volcano Plot – Cluster {group}")

    if save_fig:
        plt.savefig(f"figures/volcano_cluster_{group}.png", dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


# -------------------------------------------------------------------------
# Marker genes for PBMC visualization
# -------------------------------------------------------------------------
canonical_markers = [
    "IL7R", "CCR7", "LTB",
    "NKG7", "GNLY", "GZMB", "PRF1",
    "MS4A1", "CD79A", "CD79B",
    "LYZ", "S100A8", "S100A9", "LGALS3",
    "FCER1A", "CST3",
    "PPBP"
]

# Only keep markers present in the dataset
marker_genes = [g for g in canonical_markers if g in adata.var_names]
print(f"Using {len(marker_genes)} canonical marker genes for visualization.")

save_fig = not args.show
groupby_key = "leiden" if "leiden" in adata.obs else "louvain"


# -------------------------------------------------------------------------
# Volcano plots (optional)
# -------------------------------------------------------------------------
if args.volcano:
    # Ensure DE results exist
    compute_de_if_needed()

    # Build DE DataFrame from Scanpy results
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    de_df = pd.concat([
        pd.DataFrame({
            "group": g,
            "gene": result["names"][g],
            "logfoldchange": result.get("logfoldchanges", {g: [None]})[g],
            "pvals_adj": result.get("pvals_adj", {g: [None]})[g],
            "scores": result.get("scores", {g: [None]})[g],
        })
        for g in groups
    ], ignore_index=True)

    print("Generating volcano plots...")
    for g in sorted(de_df["group"].unique()):
        plot_volcano(de_df, group=g, save_fig=save_fig, show=args.show)
    print("Volcano plots saved.")


# -------------------------------------------------------------------------
# Other optional plots (UMAP, heatmap, dotplot, etc.)
# -------------------------------------------------------------------------
if args.umap:
    print("UMAP: plotting clusters and markers...")
    sc.pl.umap(adata, color=[groupby_key], save=("_clusters.png" if save_fig else None), show=args.show)
    if marker_genes:
        sc.pl.umap(adata, color=marker_genes[:6], ncols=3,
                   save=("_markers.png" if save_fig else None), show=args.show)

if args.heatmap and marker_genes:
    print("Heatmap: generating marker heatmap...")
    sc.tl.dendrogram(adata, groupby=groupby_key)
    sc.pl.heatmap(
        adata,
        var_names=marker_genes,
        groupby=groupby_key,
        use_raw=False,
        cmap="viridis",
        dendrogram=True,
        swap_axes=True,
        show=args.show,
        save=("_marker_heatmap.png" if save_fig else None),
    )

if args.matrixplot and marker_genes:
    print("Matrixplot: generating matrixplot...")
    sc.pl.matrixplot(
        adata,
        var_names=marker_genes,
        groupby=groupby_key,
        use_raw=False,
        cmap="viridis",
        standard_scale="var",
        show=args.show,
        save=("_marker_matrixplot.png" if save_fig else None),
    )

if args.dotplot and marker_genes:
    print("Dotplot: generating dotplot...")
    sc.pl.dotplot(
        adata,
        var_names=marker_genes,
        groupby=groupby_key,
        standard_scale="var",
        show=args.show,
        save=("_marker_dotplot.png" if save_fig else None),
    )

if args.tracksplot and marker_genes:
    print("Tracksplot: generating tracksplot...")
    sc.pl.tracksplot(
        adata,
        var_names=marker_genes,
        groupby=groupby_key,
        use_raw=False,
        show=args.show,
        save=("_marker_tracksplot.png" if save_fig else None),
    )

if args.stacked_violin and marker_genes:
    print("Stacked violin: generating stacked violin plot...")
    sc.pl.stacked_violin(
        adata,
        var_names=marker_genes,
        groupby=groupby_key,
        use_raw=False,
        swap_axes=False,
        show=args.show,
        save=("_marker_stacked_violin.png" if save_fig else None),
    )

if args.violin and marker_genes:
    print("Violin: generating multi-panel violin plots...")
    sc.pl.violin(
        adata,
        keys=marker_genes,
        groupby=groupby_key,
        rotation=90,
        multi_panel=True,
        show=args.show,
        save=("_marker_violin.png" if save_fig else None),
    )

print("Selected operations completed.")
