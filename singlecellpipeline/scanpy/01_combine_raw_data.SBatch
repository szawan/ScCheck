#!/bin/bash
#SBATCH --partition=requeue
#SBATCH --time=1-21:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --job-name=2024_10_129_01_combine_raw_data
#SBATCH --output=2024_10_129_01_combine_raw_data_%A-%a.out

srun cd $(pwd)

source activate scanpy

python3 01_combine_raw_data.py --input_path_1 /cluster/pixstor/xudong-lab/sah2p/projects/sc_data/GSE131907/GSE131907_T_lung.txt \
--input_path_2 /cluster/pixstor/xudong-lab/sah2p/projects/sc_data/GSE131907/GSE131907_N_lung.txt --output_path /cluster/pixstor/xudong-lab/sah2p/projects/sc_data/GSE131907/