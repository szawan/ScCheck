import scanpy as sc
import scipy.sparse as sp
import scanpy.external as sce
import numpy as np
import time
import tracemalloc 
# import scvi
from sklearn.metrics import silhouette_score
# from scibi.metrics import silhouette_batch, kBET, ilisi_graph, ari, nmi
import argparse
import os
import warnings
import logging
warnings.filterwarnings("ignore", category=UserWarning)
# Set up logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
log = logging.getLogger("BatchCorrection")
log.setLevel(logging.INFO)

batch_key = "Orig.ident"
# Function to calculate clustering and batch correction metrics
def compute_metrics(adata, batch_key="batch", cluster_key="Celltype"):
    metrics = {}

    # Silhouette Score for clustering
    metrics["Silhouette Score"] = silhouette_score(adata.obsm["X_pca"], adata.obs[cluster_key])

    # Batch mixing: kBET (0 = perfect mixing, 1 = poor mixing)
    metrics["kBET"] = kBET(adata, batch_key=batch_key, use_rep="X_pca")

    # Clustering preservation: ASW score
    metrics["ASW Score"] = silhouette_batch(adata, batch_key=batch_key, embed="X_pca")

    # Clustering Index: Adjusted Rand Index (ARI)
    metrics["ARI Score"] = ari(adata, cluster_key=cluster_key)

    # Normalized Mutual Information (NMI)
    metrics["NMI Score"] = nmi(adata, cluster_key=cluster_key)

    return metrics

# Function to save figures
def save_figure(adata, filename, color=[batch_key]):
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.pca(adata, n_comps=50)
    sc.tl.leiden(adata, resolution=0.5)
    # sc.tl.tsne(adata)
    # sc.tl.draw_graph(adata)
    sc.tl.louvain(adata, resolution=0.5)
    
    
    sc.pl.umap(adata, color=color, save=filename, show=False)
    sc.pl.pca(adata, color=color, save=filename, show=False)
    # sc.pl.draw_graph(adata, color=color, save=filename, show=False)
    sc.pl.umap(adata, color=['leiden'], save=filename+'leiden', show=False)
    sc.pl.umap(adata, color=['louvain'], save=filename+'louvain', show=False)
    


def correctBatch_Anndata(method, context_path, input_file, output_file):
    """ 
    Corrects batch effects in an AnnData object using Harmony, BKNN and scanroma."""
    original_dir = os.getcwd()
    os.chdir(context_path)
    log.info(f"Current working directory: {os.getcwd()}")

    try:
        # Start execution timer and memory tracking
        start_time = time.time()
        tracemalloc.start()
        # Check if input file exists
        if not os.path.exists(input_file):
            log.error(f"Error: Input file not found at {input_file}")   
            raise FileNotFoundError(f"Error: Input file not found at {input_file}")
        
        log.info(f"📂 Loading data from: {input_file}")
        adata_raw = sc.read_h5ad(input_file)
        
        #plotting before batch correction
        log.info("Plotting before batch correction...")
        save_figure(adata_raw, "before_batch_correction.png", color=[batch_key])
        
        if method == "bbknn":
            sce.pp.bbknn(adata_raw, batch_key=batch_key)
        elif method == "harmony":
            sce.pp.harmony_integrate(adata_raw, batch_key)
        elif method == "scanorama":
            sce.pp.scanorama_integrate(adata_raw, batch_key=batch_key, verbose = True)
        elif method == "combat":
            sc.pp.combat(adata_raw, key=batch_key)
        elif method == "mnn":
            sce.pp.mnn_correct(adata_raw, batch_key=batch_key)
        # elif method == "scvi":
        #     scvi.model.SCVI.setup_anndata(adata_raw, batch_key=batch_key)
        #     vae = scvi.model.SCVI(adata_raw)
        #     vae.train()
        #     adata_raw.obsm["X_scVI"] = vae.get_latents()
        else:
            raise ValueError(f"Batch correction method '{method}' is not supported.")
        
        #plotting after batch correction
        log.info("Plotting after batch correction...")
        save_figure(adata_raw, "after_batch_correction.png", color=[batch_key])
        
        # Compute metrics
        # log.info("Computing metrics...")
        # metrics = compute_metrics(adata_raw, batch_key=batch_key)
        # log.info("Metrics computed:")
        # for metric, value in metrics.items():
            # log.info(f"{metric}: {value}")
        # Print memory usage
        current, peak = tracemalloc.get_traced_memory()
        log.info(f"Current memory usage: {current / 10**6}MB; Peak: {peak / 10**6}MB")
        # Stop memory tracking
        tracemalloc.stop()
        
        # Save the corrected data
        adata_raw.write(output_file)
        log.info(f"Corrected data saved to: {output_file}")

     
    except Exception as e:
        log.error(f"An error occurred: {e}")
    finally:
        os.chdir(original_dir)
    
    
    return 1


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process and save raw AnnData object with memory profiling.")
    parser.add_argument("--context_path", type=str, required=True, help="Path to all inputs and outputs")
    parser.add_argument("--input_file", type=str, required=True, help="input .h5ad file")
    parser.add_argument("--output_file", type=str, required=True, help="output .h5ad file")
    
    args = parser.parse_args()
    method = "bbknn"  # Default method, can be changed to "harmony", "scanorama", "combat", "mnn", or "scvi"    
    # Run the function
    correctBatch_Anndata(method, args.context_path, args.input_file, args.output_file)