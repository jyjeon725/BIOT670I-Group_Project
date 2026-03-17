

# ---------------------------------------------------------------------
# 3B. Cell Type Annotation 
# ---------------------------------------------------------------------
print("Annotating clusters with biological cell types...")

celltype_annotations = {
    "0": "CD4 T cells",
    "1": "B cells",
    "2": "CD14+ Monocytes",
    "3": "NK cells",
    "4": "CD8 T cells",
    "5": "FCGR3A+ Monocytes",
    "6": "Dendritic cells",
    "7": "Megakaryocytes"
}

# Map cluster IDs → cell type names
adata.obs["cell_type"] = adata.obs["leiden"].map(celltype_annotations)

print("Cell type annotation complete.")

