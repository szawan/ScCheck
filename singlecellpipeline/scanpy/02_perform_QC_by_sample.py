import scanpy as sc
import pandas as pd
import anndata as ad
from anndata import AnnData
import argparse

def split_by_sample(adata_obj: AnnData) -> list:
    """Split an AnnData object by sample name and tissue type."""
    # Extract sample, tissue, and barcode information from observation names
    adata_obj.obs['barcode'] = adata_obj.obs_names.str.split('_').str[0]
    adata_obj.obs['Tissue'] = adata_obj.obs_names.str.split('_').str[1]
    adata_obj.obs['sample'] = adata_obj.obs_names.str.split('_').str[2]
    adata_obj.obs['phenotype'] = adata_obj.obs['sample'].apply(lambda x: "normal" if "N" in x else "tumour")
    
    # Split the AnnData object by unique sample names
    return [adata_obj[adata_obj.obs['sample'] == sample].copy() for sample in adata_obj.obs['sample'].unique()]

def perform_QC(adata_by_sample: list) -> list:
    """Perform quality control on a list of AnnData objects by sample."""
    for i, adata_sample in enumerate(adata_by_sample):
        # Identify mitochondrial, ribosomal, and hemoglobin genes
        adata_sample.var['mt'] = adata_sample.var_names.str.startswith('MT-')
        adata_sample.var["ribo"] = adata_sample.var_names.str.startswith(("RPS", "RPL"))
        adata_sample.var["hb"] = adata_sample.var_names.str.contains("^HB[^(P)]")  # Adjust prefix if necessary
        
        # Calculate quality control metrics
        sc.pp.calculate_qc_metrics(adata_sample, qc_vars=['mt', "ribo", "hb"], inplace=True)
        
        # Filter cells based on QC metrics
        adata_sample = adata_sample[
            (adata_sample.obs['pct_counts_mt'] <= 20) &
            (adata_sample.obs['n_genes_by_counts'] > 100) &
            (adata_sample.obs['n_genes_by_counts'] < 150000) &
            (adata_sample.obs['total_counts'] > 200) &
            (adata_sample.obs['total_counts'] < 10000)
        ]
        
        # Remove unwanted QC columns
        remove_columns = ['total_counts_mt', 'log1p_total_counts_mt', 
                          'total_counts_ribo', 'log1p_total_counts_ribo', 
                          'total_counts_hb', 'log1p_total_counts_hb']
        adata_sample.obs = adata_sample.obs.drop(columns=remove_columns, errors='ignore')
        
        # Store the filtered AnnData object back
        adata_by_sample[i] = adata_sample
    
    return adata_by_sample

def combine_data(adata_by_sample: list, output_path: str) -> None:
    """Combine a list of AnnData objects into one and save to file."""
    combined_adata = adata_by_sample[0].concatenate(adata_by_sample[1:])
    combined_adata.write(f"{output_path}02_QCPerformed_adata.h5ad")

def main(input_path: str, output_path: str) -> None:
    adata_obj = ad.read_h5ad(input_path)
    adata_split = split_by_sample(adata_obj)
    adata_split = perform_QC(adata_split)
    combine_data(adata_split, output_path)
    print(f"QC complete! File written to {output_path}02_QCPerformed_adata.h5ad")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quality Check by Sample and write combined AnnData object.")
    parser.add_argument('--input_path', type=str, required=True, help='Path to input AnnData file in h5ad format')
    parser.add_argument('--output_path', type=str, required=True, help='Output file path for combined AnnData')
    args = parser.parse_args()
    
    main(args.input_path, args.output_path)