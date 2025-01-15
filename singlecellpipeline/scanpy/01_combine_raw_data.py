import scanpy as sc
import pandas as pd
import anndata as ad
import argparse

def process_data(path_1,path_2, output_path):
    
    # context_path = "/scratch/sah2p/projects/G2PDeep-sc-LLM/data/GSE131907/"
    # input_1_path = context_path+"GSE131907_raw_UMI_T_lung.txt"
    # input_2_path = context_path+"GSE131907_raw_UMI_N_lung.txt"

    output_file = "01_TN.h5ad"

    counts_T = pd.read_csv(path_1, sep="\t", index_col=0)
    counts_N = pd.read_csv(path_2, sep="\t", index_col=0)

    counts = pd.concat([counts_N,counts_T],axis=1)

    print(counts.columns)

    adata = ad.AnnData(X=counts.T.values)
    adata.var_names = counts.index
    adata.obs_names = counts.columns

    adata.write(output_path+output_file)



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Process and save AnnData from input files")
    parser.add_argument('--input_path_1',type=str, required=True, help='Path to the first file')
    parser.add_argument('--input_path_2',type=str, required=True, help='Path to the second file')
    parser.add_argument('--output_path',type=str, required=True, help='output file path')
    
    args = parser.parse_args()
    
    process_data(args.input_path_1, args.input_path_2, args.output_path)
