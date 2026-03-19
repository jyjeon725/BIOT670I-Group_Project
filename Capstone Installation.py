Install Python 3 and pip if not already installed sudo apt update sudo apt install python3 python3-pip

# Install the required libraries
pip install scanpy pandas matplotlib

# Install the required libraries
pip install --user scikit-misc
pip install leidenalg

# Or with conda:
conda install -c conda-forge scanpy pandas matplotlib

# Create the required directories
mkdir -p /home/StudentFirst/git/BIOT670I-Group_Project/write/figures
mkdir -p /home/StudentFirst/git/BIOT670I-Group_Project/data/filtered_gene_bc_matrices/hg19/
mkdir -p /home/StudentFirst/git/BIOT670I-Group_Project/data/filtered_gene_bc_matrices/hg19/
#Put py script in /home/StudentFirst/git/BIOT670I-Group_Project
# Create the directory 
#mkdir -p data/filtered_gene_bc_matrices/hg19/
#Put matrix.mtx , barcodes.tsv , and features.tsv files in this folder


# Navigate to the script directory
cd /home/StudentFirst/git/BIOT670I-Group_Project/

# Run the script
# python3 your_script_name.py

