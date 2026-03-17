
# --------------------------------
# Cell Type Annotation
# --------------------------------

celltype_annotations = {
    "0": "CD4 T cells",
    "1": "NK cells",
    "2": "B cells",
    "3": "CD14+ Monocytes",
    "4": "FCGR3A+ Monocytes",
    "5": "Dendritic cells",
    "6": "Platelets",
    "7": "CD8 T cells"
}

# Map Leiden clusters to cell types
adata.obs["cell type"] = adata.obs["leiden"].astype(str).map(celltype_annotations)

# Fill clusters not annotated
adata.obs["cell type"] = adata.obs["cell type"].fillna("Unknown").astype("category")

# Plot annotated clusters
sc.pl.umap(adata, color=["leiden", "cell type"])

# --------------------------------
# Save annotated dataset
# --------------------------------

output_file = "write/pbmc3k_annotated.h5ad"
adata.write(output_file)

# Print completion message
print("Cell annotation complete.")
print("Annotated file saved as:", output_file)

