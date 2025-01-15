import scanpy as sc
import numpy as np
import anndata as ad
import argparse
from anndata import AnnData
from sklearn.model_selection import train_test_split

RANDOM_SEED = 43

def split_into_traintest(label:str, adata:AnnData):
    """Splits the data into train test splits"""
    test_size = 0.3
    # if duplication found
    # adata = adata[~adata.obs.index.duplicated(keep='first')]
    
    # Create train and test indices from the obs names
    train_indices, test_indices = train_test_split(adata.obs_names, 
                                                   test_size=0.3,
                                                   stratify=adata.obs[label],
                                                   random_state=RANDOM_SEED)

    # Split the AnnData object based on the indices
    adata_train = adata[train_indices].copy()
    adata_test = adata[test_indices].copy()
    
    print("Train set size:", adata_train.shape)
    print("Test set size:", adata_test.shape)
    return adata_train, adata_test
   
def main(input_path: str, output_path: str) -> None:
    
    adata_obj = ad.read_h5ad(input_path)
    adata_train, adata_test = split_into_traintest('phenotype',adata_obj)
    adata_train.write(output_path+"train.h5ad")
    adata_test.write(output_path+"test.h5ad")
    
    print(f"Train Test split complete! File written to {output_path}")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Train Test Split completed AnnData object written to output.")
    parser.add_argument('--input_path', type=str, required=True, help='Path to input AnnData file in h5ad format')
    parser.add_argument('--output_path', type=str, required=True, help='Output file path for combined AnnData')
    args = parser.parse_args()
    
    main(args.input_path, args.output_path)
    