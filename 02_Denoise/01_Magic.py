import argparse
import scanpy as sc
from magic import MAGIC
import scipy.sparse as sp


# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run MAGIC on single-cell RNA-seq data.")
parser.add_argument("--input_path", type=str, required=True, help="Path to input .h5ad file.")
parser.add_argument("--output_path", type=str, required=True, help="Path to save MAGIC-processed .h5ad file.")
args = parser.parse_args()

# Load the dataset
print(f"Loading data from {args.input_path}...")
adata = sc.read_h5ad(args.input_path)
adata.X = adata.raw.X.copy()
adata.var.index.name = 'genes'
adata.raw.var.columns = ['genes']


# Apply MAGIC
print("Applying MAGIC...")
magic_op = MAGIC(solver='approximate')
adata.layers['magic_approximate'] = magic_op.fit_transform(adata.X)
# adata.layers['magic_approximate'] = adata.layers['magic_approximate'].tocsr()  # Convert to CSR format
# adata.layers['magic_approximate'] = scipy.sparse.csr_matrix(adata.layers['magic_approximate'])


# Check var index issue
print("Columns in adata.var before saving:", adata.var.columns)

# Drop `_index` if it exists
# if '_index' in adata.var.columns:
#     print("Dropping `_index` from adata.var")
#     adata.var = adata.var.drop(columns=['_index_'])

# Reset index safely
# adata.var.reset_index(inplace=True, drop=True)

# Save the MAGIC-processed dataset
print(f"Saving MAGIC-processed data to {args.output_path}...")
adata.write_h5ad(args.output_path)

print("MAGIC processing completed successfully.")