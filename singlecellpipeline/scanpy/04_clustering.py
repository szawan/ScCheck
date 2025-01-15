import scanpy as sc
import pandas as pd
import anndata as ad
from anndata import AnnData
import argparse


def reduce_dimensions(adata_obj:AnnData) -> AnnData:
    """ Calculates principle component analysis"""
    sc.tl.pca(adata_obj)
    sc.pl.pca_variance_ratio(adata_obj, n_pcs=50, log=True, save="04_clustering_PCA.png")

    return adata_obj

def visualize_data(adata_obj:AnnData) -> AnnData:
    """Visualizes the data using tsne and umap"""
    sc.pp.neighbors(adata_obj)
    sc.tl.umap(adata_obj)
    sc.tl.tsne(adata_obj)
    sc.pl.umap(adata_obj, color='sample', save="04_clustering_umap.pdf")
    sc.pl.tsne(adata_obj, color='sample',save="04_clustering_tsne.pdf")
    return adata_obj

def cluster_data(adata_obj:AnnData)-> AnnData:
    """This function uses louvain and leiden to generate clusters"""
    li_res = [0.01, 0.02, 0.05, 0.1, 0.5] 
    for res in li_res:
        sc.tl.leiden(adata_obj,  key_added='leiden_res'+str(res),resolution=res)
        sc.tl.louvain(adata_obj,  key_added='louvain_res'+str(res),resolution=res)

    sc.pl.umap(adata_obj, color=['louvain_res0.01','louvain_res0.02','louvain_res0.05','louvain_res0.1','louvain_res0.5'],save="04_clustering_louvain_umap.pdf")
    sc.pl.tsne(adata_obj, color=['louvain_res0.01','louvain_res0.02','louvain_res0.05','louvain_res0.1','louvain_res0.5'],save="04_clustering_louvain_tsne.pdf")
    
    sc.pl.umap(adata_obj, color=['leiden_res0.01','leiden_res0.02','leiden_res0.05','leiden_res0.1','leiden_res0.5'],save="04_clustering_leiden_umap.pdf")
    sc.pl.tsne(adata_obj, color=['leiden_res0.01','leiden_res0.02','leiden_res0.05','leiden_res0.1','leiden_res0.5'],save="04_clustering_leiden_tsne.pdf")

    return adata_obj

def reassess_QC(adata_obj:AnnData) -> AnnData:
    """Reassess the quanlity and filter doublets if any"""
    
    # TODO (Sania): make selected_res an argument in the script
    
    selected_res = "louvain_res0.1"
    adata_obj.obs["predicted_doublet"] = adata_obj.obs["predicted_doublet"].astype("category")
    sc.pl.umap(
    adata_obj,
    color=[selected_res, "predicted_doublet", "doublet_score"],
    wspace=0.1,
    save="04_clustering_predicted_doublets.pdf")
    
    # Subset the Anndata to exclude cells predicted as doublets
    adata_obj =adata_obj[~adata_obj.obs['predicted_doublet'].to_numpy()].copy()
    
    sc.pl.umap(
    adata_obj, 
    color=["louvain_res0.1", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"], 
    wspace=0.5, ncols=2,
    save = "04_clustering_final_clusters.pdf"
)
    
    return adata_obj
    
def main(input_path:str, output_path:str)-> None:
    
    adata_obj = ad.read_h5ad(input_path)
    adata_obj = reduce_dimensions(adata_obj)
    adata_obj = visualize_data(adata_obj)
    adata_obj = cluster_data(adata_obj)
    adata_obj = reassess_QC(adata_obj)
    adata_obj.write(output_path+"04_clusters.h5ad")
    print(f"Clustering complete! File written to {output_path}04_clusters.h5ad")

    pass
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clustering completed AnnData object written to output.")
    parser.add_argument('--input_path', type=str, required=True, help='Path to input AnnData file in h5ad format')
    parser.add_argument('--output_path', type=str, required=True, help='Output file path for combined AnnData')
    args = parser.parse_args()
    
    main(args.input_path, args.output_path)
   