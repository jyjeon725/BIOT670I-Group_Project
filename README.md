# Group Project: Project 2 Group 2b: Single-Cell RNA Analysis

## Preparation
- **GitHub Setup**
  - Team members: YoungJu Jeon (Preprocessing and Quality Control), Berke Sahbazoglu (Dimensionality Reduction and Clustering), Luxurie Mills (Cell Type Annotation), Christian Gifueroa-Perez (Differential Expression (DE) Analysis)

- **Enviroment Setup Windows**
```windows
python3 -m venv scanpy_env
#Activate environment
source scanpy_env/bin/activate
#Install dependencies
pip install -r requirements.txt
```
- **Enviroment Setup**
```bash
/usr/bin/python3 -m venv scanpy_env
#Activate environment
source scanpy_env/bin/activate
#Install dependencies
pip install --upgrade pip setuptools wheel
pip install scanpy
pip install scikit-misc

#Install dependecies options 2:
pip install -r requirements.txt
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

python3 cell_annotation.py (LM)

python3

```
