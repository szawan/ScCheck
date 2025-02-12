import os
import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse
import matplotlib.pyplot as plt

def quality_check(adata, method_name, dataset_name, output_dir="quality_checks"):
    """
    Perform quality checks on denoised data and save results.

    Parameters:
    - adata: AnnData object
    - method_name: Name of denoising method (e.g., MAGIC, DCA)
    - dataset_name: Name of dataset (e.g., SRP173393)
    - output_dir: Directory where results will be saved
    """

    # Create structured output folder
    save_path = os.path.join(output_dir, f"{dataset_name}_{method_name}")
    os.makedirs(save_path, exist_ok=True)

    # Extract Data
    raw_data = adata.raw.X if adata.raw is not None else None
    denoised_data = adata.X if "magic_approximate" in adata.layers else adata.X

# Convert sparse matrix to dense (if applicable)
    if scipy.sparse.issparse(denoised_data):
        denoised_data = denoised_data.toarray()
        


    # Compute Statistics
    stats = {
        "Mean": np.mean(denoised_data),
        "Standard Deviation": np.std(denoised_data),
        "Min": np.min(denoised_data),
        "Max": np.max(denoised_data),
        "Are All Values Integers": np.all(denoised_data == denoised_data.astype(int))
    }

    # Save statistics to file
    stats_df = pd.DataFrame(stats, index=[0])
    stats_df.to_csv(os.path.join(save_path, "denoising_stats.csv"), index=False)
    
    # Print summary
    print(f"Quality Check Results for {dataset_name} with {method_name}:")
    for k, v in stats.items():
        print(f"  {k}: {v}")

    # Create Histogram
    plt.figure(figsize=(8, 5))
    plt.hist(denoised_data.flatten(), bins=50, alpha=0.7, color="blue", label="Denoised")
    plt.title(f"{method_name} Denoised Data Distribution")
    plt.xlabel("Expression Values")
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(os.path.join(save_path, "histogram.png"))
    plt.close()

    # Violin Plot for top genes
    top_genes = np.argsort(np.var(denoised_data, axis=0))[-10:]  # Top 10 most variable genes
    plt.figure(figsize=(10, 6))
    sc.pl.violin(adata, keys=adata.var_names[top_genes], jitter=0.4, multi_panel=True, save=f"_violin_{method_name}.png")
    
    print(f"Saved results in: {save_path}")

if __name__ == "__main__":
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Run quality checks on a denoised dataset.")
    parser.add_argument("--input", type=str, required=True, help="Path to denoised .h5ad file")
    parser.add_argument("--method", type=str, required=True, help="Denoising method (e.g., MAGIC, DCA)")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name (e.g., SRP173393)")
    args = parser.parse_args()

    # Load dataset
    adata = sc.read_h5ad(args.input)
    
    # Run Quality Checks
    quality_check(adata, args.method, args.dataset)