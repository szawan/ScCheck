import scanpy as sc
import argparse
import os
import sys
import time
import tracemalloc  # For memory profiling
from scipy import sparse
import numpy as np
import doubletdetection
import warnings
import logging
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")
# Set up logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
log = logging.getLogger("qualityControl")
log.setLevel(logging.INFO)



def filter_anndata(context_path,input_file, output_file):
    """
    Filters raw counts and removes doublets, empty droplets, and low-quality cells.
    Tracks execution time and memory usage.
    """
    original_dir = os.getcwd()
    os.chdir(context_path)
    log.info(f"Current working directory: {os.getcwd()}")
    try:
        # Start execution timer and memory tracking
        start_time = time.time()
        tracemalloc.start()

        # Check if input file exists
        input_path = input_file
        if not os.path.exists(input_path):
            log.error(f"Error: Input file not found at {input_path}")   
            raise FileNotFoundError(f"Error: Input file not found at {input_path}")
        
        log.info(f"📂 Loading data from: {input_path}")
        adata_raw = sc.read_h5ad(input_path)
        
        # Identify mitochondrial and plastid genes
        mt_genes = [gene for gene in adata_raw.var_names if gene.startswith("ATMG")]
        pt_genes = [gene for gene in adata_raw.var_names if gene.startswith("ATCG")]
        # Print the number of detected genes
        log.info(f"Found {len(mt_genes)} mitochondrial genes and {len(pt_genes)} plastid genes.")

        # Compute the fraction of mitochondrial gene counts per cell
        adata_raw.obs["percent_mt"] = (
            np.sum(adata_raw[:, mt_genes].X, axis=1) / np.sum(adata_raw.X, axis=1)
        ) * 100

        # Compute the fraction of plastid gene counts per cell
        adata_raw.obs["percent_pt"] = (
            np.sum(adata_raw[:, pt_genes].X, axis=1) / np.sum(adata_raw.X, axis=1)
        ) * 100

        adata_raw.obs["n_genes_by_counts"] = (adata_raw.X > 0).sum(axis=1)
        # Compute total UMI counts per cell (nCount_RNA equivalent)
        adata_raw.obs["total_counts"] = adata_raw.X.sum(axis=1)


        # Visualize distributions before filtering
        sc.pl.violin(adata_raw, ["percent_mt", "percent_pt"], 
                        jitter=0.4, log=True, \
                        ylabel="Percentage of Counts", save="before_ptmt.png")
        sc.pl.violin(adata_raw,["n_genes_by_counts","total_counts"], jitter=0.4, 
                        log=True, ylabel="Counts",save="before_counts.png")

        # Step 6: Apply IQR-based filtering
        Q1 = np.percentile(adata_raw.obs["n_genes_by_counts"], 25)  # First quartile
        Q3 = np.percentile(adata_raw.obs["n_genes_by_counts"], 75)  # Third quartile
        IQR = Q3 - Q1

        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        adata_raw = adata_raw[
                    (adata_raw.obs['percent_mt'] <= 10) &
                    (adata_raw.obs['percent_pt'] <= 5) &
                    (adata_raw.obs['n_genes_by_counts'] > lower_bound) &
                    (adata_raw.obs['n_genes_by_counts'] < upper_bound) &
                    (adata_raw.obs['total_counts'] > 200) &
                    (adata_raw.obs['total_counts'] < 10000)]
        # Visualize distributions before filtering
        sc.pl.violin(adata_raw, ["percent_mt", "percent_pt"], 
                        jitter=0.4, log=True, \
                        ylabel="Percentage of Counts", save="after_ptmt.png")
        sc.pl.violin(adata_raw,["n_genes_by_counts","total_counts"], jitter=0.4, 
                        log=True, ylabel="Counts",save="after_counts.png")

        log.info(f"Remaining cells after filtering: {adata_raw.n_obs}")
        
        # Step 8: Remove empty droplets and doublets
        clf = doubletdetection.BoostClassifier(
        n_iters=10, 
        clustering_algorithm="louvain", 
        standard_scaling=True,
        pseudocount=0.1,
        n_jobs=-1)
        
        doublets = clf.fit(adata_raw.X).predict(p_thresh=1e-16, voter_thresh=0.5)
        doublet_score = clf.doublet_score()
        adata_raw.obs["doublet"] = doublets
        adata_raw.obs["doublet_score"] = doublet_score
        print(adata_raw.obs["doublet"].value_counts())
        
        sc.tl.pca(adata_raw)  # Ensure PCA is computed
        sc.pp.neighbors(adata_raw)  # Compute neighbors graph
        sc.tl.umap(adata_raw)
        sc.pl.umap(adata_raw, color=["doublet", "doublet_score"],save="doublet_umap.png")
        adata_raw = adata_raw[adata_raw.obs['doublet'] == 0, :]
        log.info(f"Remaining cells after doubletRemoval: {adata_raw.n_obs}")
        
        adata_raw.layers['counts'] = adata_raw.X.copy()
        # Normalizing to median total counts
        sc.pp.normalize_total(adata_raw)
        # Logarithmize the data:
        sc.pp.log1p(adata_raw)
        
        #Selecting HVG
        sc.pp.highly_variable_genes(adata_raw, n_top_genes=6000, batch_key="Orig.ident")
        sc.pl.highly_variable_genes(adata_raw, save = "HVG.png")
        # TODO(Sania): Add flag to keep all genes
        #Filter to keep only highly variable genes
        adata_raw = adata_raw[:, adata_raw.var["highly_variable"]]
        
        log.info("Scaling data to normalize variance...")
        sc.pp.scale(adata_raw, max_value=10)
        
        output_path = output_file
        log.info(f"📂 Saving filtered AnnData object to: {output_path}")
        adata_raw.write_h5ad(output_path)
        log.info("✅ Filtering and saving completed successfully.")
        
        # End timing and memory tracking
        end_time = time.time()
        current_mem, peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        # Print execution time and memory usage
        log.info("\n🕒 Execution Time: {:.2f} seconds".format(end_time - start_time))
        log.info("🧠 Peak Memory Usage: {:.2f} MB".format(peak_mem / (1024 * 1024)))
        

    except FileNotFoundError as e:
        log.error(f"🚨 FileNotFoundError: {e}")
        sys.exit(1)
    except Exception as e:
        log.error(f"🚨 An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process and save raw AnnData object with memory profiling.")
    parser.add_argument("--context_path", type=str, required=True, help="Path to all inputs and outputs")
    parser.add_argument("--input_file", type=str, required=True, help="input .h5ad file")
    parser.add_argument("--output_file", type=str, required=True, help="output .h5ad file")
    
    args = parser.parse_args()

    # Run the function
    filter_anndata(args.context_path, args.input_file, args.output_file)