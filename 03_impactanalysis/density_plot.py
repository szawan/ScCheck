import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
# from magic import MAGIC
import scipy.sparse as sp
import scanpy.external as sce

context_path = "/home/sah2p/ondemand/singlecell_data/"
organism = 'Arabidopsis/'
adata = sc.read_h5ad(context_path+organism+"SRP166333_magic_denoised.h5ad")

# Assign MAGIC data as the new main matrix (if needed)
adata.X= adata.layers['magic_approximate']

# Set figure size
plt.figure(figsize=(10, 5))

# Plot KDE distribution for raw and denoised data
sns.kdeplot(np.ravel(adata.raw.X.A), label="Raw Counts", fill=True, alpha=0.5)
sns.kdeplot(np.ravel(adata.X.A), label="Denoised", fill=True, alpha=0.5)

# Labels and title
plt.xlabel("Gene Expression Values")
plt.ylabel("Density")
plt.legend()
plt.title("Gene Expression Distribution: Raw vs Denoised")

# Save the figure
output_path = context_path+organism+"gene_expression_distribution.png"  # Change to your desired save path
plt.savefig(output_path, dpi=300, bbox_inches="tight")

# Show the figure (optional)
plt.show()