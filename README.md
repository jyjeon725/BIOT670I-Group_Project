# Group Project: Project 2 Group 2b: Single-Cell RNA Analysis

## Preparation
- **GitHub Setup**
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

- **Command execution**
```bash
python3 scanpy_script_qc.py

python3

python3

#PCA uses principal components to plot genes with most significant variation and usually will only show a couple dozen PCs. #
#Each dot is a single cell and physical proximity means the cells have similar gene expressions. #
#Clusters hint to different cell types or metabolic states and overlapping cells hint to similar gene expressions #

sc.pp.pca(adata)
sc.pl.pca(adata)


#If you are curious which genes are leading to the major principal components use the following to rank them#

sc.pl.pca_loadings(adata)


#Another view of the data involves seeing how much each PC contributes to the overall variance for the number of PCs we had previously preselected. #

sc.pl.pca_variance_ratio(adata, log=True)

#Build a neighborhood graph connecting each cell to its nearest neighbors based on PCA coordinates #
#Cells with similar gene expression are connected. Here each cell is connected to 10 most similar other cells and 40 principal components are used. These numbers can be tweaked if needed. #

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)


# Cells are clustered using the Leiden algorithm which uses the nearest neighbor graph from the previous step#
#Similar cells are grouped-clustered together. #
#Higher than 1 lead to numerous small clusters while closer to zero leads to few large clusters. #
#UMAP converts cell and gene data into x and y coordinates to make a scatter plot. Then each cluster gets a different color for easier visualization. #


sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden')

#. tSNE can be preferred for smaller datasets with few close clusters. UMAP can be used for larger datasets and is faster. #

sc.tl.leiden(adata, resolution=0.5)
sc.tl.tsne(adata)
sc.pl.umap(adata, color='leiden')


python3


```
