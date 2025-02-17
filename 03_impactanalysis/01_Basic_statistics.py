import scanpy as sc
# from magic import MAGIC
import scipy.sparse as sp
import scanpy.external as sce
import numpy as np
import pandas as pd
import os

context_path = "/home/sah2p/ondemand/singlecell_data/"
organism = 'Maize/'
dataset = 'SRP335180'
denoised_adata = sc.read_h5ad(context_path+organism+dataset+"_magic_denoised.h5ad")
raw_adata = sc.read_h5ad(context_path+organism+"SRP335180.h5ad")

 # Extract Data
raw_data = raw_adata.raw.X
denoised_data = denoised_adata.layers['magic_approximate']


# Convert sparse matrix to dense (if applicable)
if sp.issparse(denoised_data):
    denoised_data = denoised_data.toarray()
        
if sp.issparse(raw_data):
    raw_data = raw_data.toarray()
        
# Compute Statistics
stats_d = {
    "Mean": np.mean(denoised_data),
    "Standard Deviation": np.std(denoised_data),
    "Min": np.min(denoised_data),
    "Max": np.max(denoised_data),
    "Are All Values Integers": np.all(denoised_data == denoised_data.astype(int))
}

stats_r = {"Mean": np.mean(raw_data),
    "Standard Deviation": np.std(raw_data),
    "Min": np.min(raw_data), 
    "Max": np.max(raw_data),
    "Are All Values Integers": np.all(raw_data == raw_data.astype(int))
}

print("Denoised Data Statistics:")
for key, value in stats_d.items():
    print(f"{key}: {value}")
print("\nRaw Data Statistics:")
for key, value in stats_r.items():
    print(f"{key}: {value}")

raw_stats_df = pd.DataFrame(stats_r, index=[0])
raw_stats_df.to_csv(os.path.join(context_path+organism, (dataset+"_raw_stats.csv")), index=False)

denoised_stats_df = pd.DataFrame(stats_d, index=[0]) 
denoised_stats_df.to_csv(os.path.join(context_path+organism, (dataset+"_denoised_stats.csv")), index=False)
