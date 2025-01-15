import scanpy as sc
import pandas as pd
import anndata as ad
from anndata import AnnData
import argparse

def preprocess_data(adata_obj)->AnnData:
    """ Performs preprocessing on the data"""
    
    #doublet detection using scrublet
    sc.pp.scrublet(adata_obj, batch_key="sample")
    
    #Normaliza the data
    adata_obj.layers['counts'] = adata_obj.X.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(adata_obj)
    # Logarithmize the data:
    sc.pp.log1p(adata_obj)
    
    return adata_obj

def select_features(adata_obj, output_path) -> AnnData:
    """Selecting 2000 highly variable genes"""
    
    sc.pp.highly_variable_genes(adata_obj, n_top_genes=2000, batch_key="sample")
    sc.pl.highly_variable_genes(adata_obj, save = "HVG.png")
    adata_obj.write(output_path+"03_HVG_adata.h5ad")
    
    return adata_obj

def main(input_path: str, output_path: str) -> None:
    adata_obj = ad.read_h5ad(input_path)
    adata_obj = preprocess_data(adata_obj)
    adata_obj = select_features(adata_obj, output_path)
    print(f"Feature Engineering complete! File written to {output_path}03_HVG_adata.h5ad")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Feature Engineering completed AnnData object written to output.")
    parser.add_argument('--input_path', type=str, required=True, help='Path to input AnnData file in h5ad format')
    parser.add_argument('--output_path', type=str, required=True, help='Output file path for combined AnnData')
    args = parser.parse_args()
    
    main(args.input_path, args.output_path)
    