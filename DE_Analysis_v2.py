#!/usr/bin/env python

"""
Differential expression and visualization on 10x Genomics PBMC 3k

Requirements:
    pip install scanpy anndata matplotlib seaborn leidenalg
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# 1. Load data (built‑in PBMC 3k from Scanpy)
# ---------------------------------------------------------------------
print("Loading PBMC3k dataset...")
adata = sc.datasets.pbmc3k()  # 2700 PBMCs from 10x Genomics

# ---------------------------------------------------------------------
# 2. Basic preprocessing and clustering
# ---------------------------------------------------------------------
print("Preprocessing...")

# Filter low-quality genes/cells (light filtering for demo)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Mitochondrial content
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# (Optional) filter high‑mito or low‑count cells
adata = adata[adata.obs["pct_counts_mt"] < 15, :].copy()

# Normalize and log‑transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor="seurat")

# Scale and regress out unwanted sources of variation
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver="arpack")

# Neighborhood graph and clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden")

# ---------------------------------------------------------------------
# 3. Differential expression (cluster markers)
# ---------------------------------------------------------------------
print("Running differential expression (rank_genes_groups)...")

# Use Wilcoxon rank-sum test to find markers per cluster
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    n_genes=50,  # top 50 per cluster
    use_raw=False,
)

# Save DE results to a table (optional)
def rank_genes_to_df(adata, groupby="leiden"):
    """Convert rank_genes_groups result to a tidy DataFrame."""
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    dfs = []
    for g in groups:
        df = pd.DataFrame(
            {
                "group": g,
                "gene": result["names"][g],
                "logfoldchange": result["logfoldchanges"][g],
                "pvals_adj": result["pvals_adj"][g],
                "scores": result["scores"][g],
            }
        )
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


de_df = rank_genes_to_df(adata)
de_df.to_csv("pbmc3k_DE_results_by_cluster.csv", index=False)
print("Saved DE results to pbmc3k_DE_results_by_cluster.csv")

# ---------------------------------------------------------------------
# 3B. Cell type annotation (Luxurie Mills add)
# ---------------------------------------------------------------------

print("Annotating clusters.....")

celltype_annotations = {
    "0": "T cells",
    "1": "B cells",
    "2": "Cytotoxic T/NK cells",
    "3": "Dendritic cells",
    "4": "Monocytes",
    "5": "Platelets/ Megakaryocytes"
}

# Map Leiden clusters to cell type labels
adata.obs["cell_type"] = adata.obs["leiden"].map(celltype_annotations)

print("Cell type annotation added to dataset.")

# -------------------------------------------------------p
# CHECK OUTPUT (this is the part for step 2)
# -------------------------------------------------------

print("\nCell type annotated data review:")
print(adata.obs.head())

#Luxurie END
# ---------------------------------------------------------------------
# 4. Choose marker genes for visualization
#    You can either:
#      - use canonical PBMC markers, or
#      - take top DE genes from rank_genes_groups.
# ---------------------------------------------------------------------

# Canonical PBMC markers (example set)
canonical_markers = [
    # T cells
    "IL7R", "CCR7", "LTB",
    # Cytotoxic T / NK
    "NKG7", "GNLY", "GZMB", "PRF1",
    # B cells
    "MS4A1", "CD79A", "CD79B",
    # Monocytes
    "LYZ", "S100A8", "S100A9", "LGALS3",
    # Dendritic cells
    "FCER1A", "CST3",
    # Megakaryocytes / platelets
    "PPBP"
]

# Filter to genes that actually exist in the dataset
marker_genes = [g for g in canonical_markers if g in adata.var_names]
print(f"Using {len(marker_genes)} marker genes for visualization.")

# Alternatively, you could pick top N DE genes per cluster:
# top_de_genes = (
#     de_df.sort_values(["group", "scores"], ascending=[True, False])
#     .groupby("group")
#     .head(5)["gene"]
#     .unique()
#     .tolist()
# )

# ---------------------------------------------------------------------
# 5. UMAP colored by clusters and marker genes (quick overview)
# ---------------------------------------------------------------------
print("Plotting UMAPs...")
sc.pl.umap(adata, color=["leiden"], save="_clusters.png", show=False)
if marker_genes:
    sc.pl.umap(adata, color=marker_genes[:6], ncols=3, save="_markers.png", show=False)

# ---------------------------------------------------------------------
# 6. Heatmap of marker genes across clusters
# ---------------------------------------------------------------------
print("Generating heatmap...")

if marker_genes:
    sc.pl.heatmap(
        adata,
        var_names=marker_genes,
        groupby="leiden",
        use_raw=False,
        cmap="viridis",
        dendrogram=True,
        swap_axes=True,
        show=False,
        save="_marker_heatmap.png",
    )

# ---------------------------------------------------------------------
# 7. Violin plots of marker genes by cluster
# ---------------------------------------------------------------------
print("Generating violin plots...")

if marker_genes:
    sc.pl.violin(
        adata,
        keys=marker_genes,
        groupby="leiden",
        rotation=90,
        multi_panel=True,
        show=False,
        save="_marker_violin.png",
    )

# ---------------------------------------------------------------------
# 8. Dot plot of marker genes by cluster
# ---------------------------------------------------------------------
print("Generating dot plots...")

if marker_genes:
    sc.pl.dotplot(
        adata,
        var_names=marker_genes,
        groupby="leiden",
        standard_scale="var",  # scale per gene
        show=False,
        save="_marker_dotplot.png",
    )

print("All plots saved in the current directory (figures/ if Scanpy uses default). Done.")
