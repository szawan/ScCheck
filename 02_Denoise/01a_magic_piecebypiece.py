import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from magic import MAGIC

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run MAGIC on single-cell RNA-seq data.")
parser.add_argument("--input_path", type=str, required=True, help="Path to input .h5ad file.")
parser.add_argument("--output_path", type=str, required=True, help="Path to save MAGIC-processed .h5ad file.")
args = parser.parse_args()

# Load the dataset
print(f"Loading data from {args.input_path}...")
adata = sc.read_h5ad(args.input_path)
adata.var.index.name = 'genes'
adata.raw.var.columns = ['genes']
# Convert the count matrix to a Pandas DataFrame
print("Extracting raw counts...")

# # Convert to dense if sparse
if isinstance(adata.X, np.ndarray):
    df_counts = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
else:
    df_counts = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)

# df_counts = pd.DataFrame(adata.raw.X, index=adata.obs_names, columns=adata.var_names)

# Apply MAGIC
print("Applying MAGIC...")
magic_op = MAGIC(solver='approximate')
df_magic = magic_op.fit_transform(df_counts)

# df_magic.to_csv(args.output_path)
# Convert back to numpy for storage in AnnData
magic_matrix = df_magic.to_numpy()

# If the original matrix was sparse, convert back to sparse
if isinstance(adata.X, scipy.sparse.spmatrix):
    adata.X = scipy.sparse.csr_matrix(adata.X)
# Verify update
print("Updated adata.X shape:", adata.X.shape)

# Save the updated AnnData object
print(f"Saving updated AnnData to {args.output_path}...")
adata.write_h5ad(args.output_path)

print("MAGIC processing completed successfully.")