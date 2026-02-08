# Group Project: Project 2 Group 2b: Single-Cell RNA Analysis
  - Team members: YoungJu Jeon (Preprocessing and Quality Control), Berke Sahbazoglu (Dimensionality Reduction and Clustering), Luxurie Mills (Cell Type Annotation), Christian Gifueroa-Perez (Differential Expression (DE) Analysis)


## Preprocessing and Quality Control
- **Enviroment Setup**

```bash
/usr/bin/python3 -m venv scanpy_env
source scanpy_env/bin/activate
pip install --upgrade pip setuptools wheel
pip install scanpy
pip install scikit-misc
```
- **Data Download**
```bash
mkdir -p data write
cd data
test -f pbmc3k_filtered_gene_bc_matrices.tar.gz || curl https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -o pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

- **QC command execution**
```bash
python3 scanpy_script_qc.py
```

- **Outputs**
1. Highest expressed genes across cells
- Highly expressed genes were inspected to ensure no single gene dominated total counts across cells, which could indicate technical artifacts or low-quality cells.
- This plot shows the fraction of total counts contributed by the most highly expressed genes across cells. It is used to detect potential technical artifacts or low-quality cells dominated by a small number of genes.
![QC highest_expr_genes](write/figures/highest_expr_genes.png)

2. Violin plot
- Violin plots show the distribution of the number of detected genes per cell (n_genes_by_counts), total UMI counts per cell (total_counts), and the percentage of mitochondrial gene expression (pct_counts_mt). These metrics were used to assess cell quality and define filtering thresholds.
- Violin plots allow visualization of the full distribution of QC metrics across cells, making it easier to identify outliers and define appropriate filtering thresholds.
![QC violin plot](write/figures/violin_plot.png)

3. Scatter plot
- The left panel shows the relationship between total UMI counts and the percentage of mitochondrial gene expression (pct_counts_mt), while the right panel shows the relationship between total UMI counts and the number of detected genes (n_genes_by_counts). These plots were used to identify low-quality cells, doublets, and to guide quality control filtering thresholds.
- These scatter plots allow us to visually assess relationships between QC metrics and to justify filtering thresholds for removing low-quality cells and potential doublets.
![QC scatter](write/figures/scatter.png)

4. HVG (Highly Variable Genes) Filter genes dispersion
- Mean gene expression is plotted against gene expression variance. Highly variable genes (black) are highlighted, as they capture the most informative biological variation and are used for downstream dimensionality reduction and clustering.
- Highly variable genes capture the most biologically informative variation while reducing noise from ubiquitously expressed or low-variance genes.
![QC filter_genes_dispersion](write/figures/filter_genes_dispersion.png)

