#!/bin/bash
#SBATCH --partition=requeue
#SBATCH --time=1-21:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --job-name=2024_10_29_02_perform_QC_by_sample
#SBATCH --output=2024_10_29_02_perform_QC_by_sample_%A-%a.out

srun cd $(pwd)
source activate scanpy

python3 02_perform_QC_by_sample.py \
--input_path /cluster/pixstor/xudong-lab/sah2p/projects/sc_data/GSE131907/01_TN.h5ad \
--output_path /cluster/pixstor/xudong-lab/sah2p/projects/sc_data/GSE131907/