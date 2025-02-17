import scanpy as sc
import argparse
import os
import sys
import time
import tracemalloc  # For memory profiling
from scipy import sparse
import numpy as np


def process_anndata(input_path, output_path):
    """
    Process and save a raw AnnData object by removing unwanted columns and structures.
    Tracks execution time and memory usage.
    """
    try:
        # Start execution timer and memory tracking
        start_time = time.time()
        tracemalloc.start()

        # Check if input file exists
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Error: Input file not found at {input_path}")

        print(f"📂 Loading data from: {input_path}")
        adata = sc.read_h5ad(input_path)

        # Ensure raw data exists
        if adata.raw is None:
            print("⚠️ Warning: No .raw data found! Using .X instead.")
            adata_raw = sc.AnnData(adata.X.copy())  # Fallback to X if raw isn't available
        else:
            adata_raw = sc.AnnData(adata.raw.X.copy())  # Use raw counts

        # Ensure raw data is in sparse format
        # if not sparse.issparse(adata_raw.X):
        #     adata_raw.X = sparse.csr_matrix(adata_raw.X)

        # # Convert to float32 to save space
        # adata_raw.X = adata_raw.X.astype(np.float32)

        # Keep all obs columns except the specified ones
        columns_to_remove = ['nCount_RNA', 'nFeature_RNA', 'Percent.mt', 'Seurat_clusters']
        adata_raw.obs = adata.obs.drop(columns=columns_to_remove, errors='ignore')  # Remove specific columns

        # Keep the original gene metadata (var)
        adata_raw.var = adata.var.copy()  # Preserve all gene information

        # Remove processed data that isn't needed
        adata_raw.uns = {}  # Remove unstructured data
        adata_raw.obsm = {}  # Remove dimensionality reductions (PCA, UMAP)
        adata_raw.obsp = {}  # Remove pairwise distances
        adata_raw.varm = {}  # Remove additional variable metadata

        # Print summary to confirm
        print(adata_raw)

        # Ensure output directory exists before saving
        output_dir = os.path.dirname(output_path)
        os.makedirs(output_dir, exist_ok=True)

        # Save the modified AnnData object
        adata_raw.write_h5ad(output_path)

        print(f"✅ Processed data successfully saved to: {output_path}")

        # End timing and memory tracking
        end_time = time.time()
        current_mem, peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        # Print execution time and memory usage
        print("\n🕒 Execution Time: {:.2f} seconds".format(end_time - start_time))
        print("🧠 Peak Memory Usage: {:.2f} MB".format(peak_mem / (1024 * 1024)))

    except FileNotFoundError as e:
        print(f"\n🚨 File Error: {e}", file=sys.stderr)
        sys.exit(1)

    except OSError as e:
        print(f"\n🚨 OS Error: Could not create/save file. Check directory permissions: {e}", file=sys.stderr)
        sys.exit(1)

    except Exception as e:
        print(f"\n🚨 Unexpected Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process and save raw AnnData object with memory profiling.")
    parser.add_argument("--input_path", type=str, required=True, help="Path to the input .h5ad file")
    parser.add_argument("--output_path", type=str, required=True, help="Path to save the output .h5ad file")
    
    args = parser.parse_args()

    # Run the function
    process_anndata(args.input_path, args.output_path)