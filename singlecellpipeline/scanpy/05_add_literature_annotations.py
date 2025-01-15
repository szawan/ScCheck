"""This script refers to annotation which were provided by GSE131907"""
import scanpy as sc
import numpy as np
import anndata as ad
import pandas as pd
from anndata import AnnData
import argparse

def assign_literature_annotations(adata_obj:AnnData, df_annotations) -> AnnData:
    df_nlung = df_annotations
    combined_adata = adata_obj
    dict_lungn = df_nlung.set_index('Barcode')['Cell_type'].to_dict()
    dict_lungn_refined = df_nlung.set_index('Barcode')['Cell_type.refined'].to_dict()
    dict_lungn_subtype = df_nlung.set_index('Barcode')['Cell_subtype'].to_dict()
    combined_adata.obs['lcell_type'] = combined_adata.obs['barcode'].map(dict_lungn)
    combined_adata.obs['lcell_type.refined'] = combined_adata.obs['barcode'].map(dict_lungn_refined)
    combined_adata.obs['lcell_subtype'] = combined_adata.obs['barcode'].map(dict_lungn_subtype)
    # sc.pl.umap(combined_adata, color=["lcell_subtype"], ncols=1, save="05_annotations_subtypes.pdf")
    # sc.pl.umap(combined_adata, color=["lcell_type"], ncols=1, save="05_annotations_types.pdf")
    return combined_adata

def main(input_path:str, input_path_ann:str, output_path:str)-> None:
    
    adata_obj = ad.read_h5ad(input_path)
    df_ann = pd.read_csv(input_path_ann,sep='\t')
    adata_obj = assign_literature_annotations(adata_obj, df_ann)
    adata_obj.write(output_path+"05a_literature_annotations.h5ad")
    print(f"Annotation from literature complete! File written to {output_path}05a_literature_annotations.h5ad")

    pass
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="literature Annotation completed AnnData object written to output.")
    parser.add_argument('--input_path', type=str, required=True, help='Path to input AnnData file in h5ad format')
    parser.add_argument('--input_path_ann', type=str, required=True, help='Path to input annotation file in csv/txt format')
    parser.add_argument('--output_path', type=str, required=True, help='Output file path for combined AnnData')
    args = parser.parse_args()
    
    main(args.input_path, args.input_path_ann, args.output_path)
   