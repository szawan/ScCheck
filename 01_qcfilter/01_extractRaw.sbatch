#!/bin/bash
#SBATCH --partition=general,requeue,logical_cpu,logical_cpu2,interactive
#SBATCH --time=1-21:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --job-name=SRP166333_01_extractRaw
#SBATCH --output=%x_%j.out  # job name and job ID in output

# Move to current working directory
cd "$(pwd)"

# Activate conda environment
source /home/sah2p/ondemand/tools/miniconda3/bin/activate scanpy

# Define paths
INPUT_PATH="/home/sah2p/ondemand/singlecell_data/Arabidopsis/SRP166333.h5ad"
OUTPUT_PATH="/home/sah2p/ondemand/singlecell_data/Arabidopsis/SRP166333_raw.h5ad"

# Run the Python script
python 01_extractRaw.py --input_path $INPUT_PATH --output_path $OUTPUT_PATH

# Print job completion message
echo "✅ Job Completed at $(date)"