#!/bin/bash
#SBATCH --partition=general,requeue,logical_cpu,logical_cpu2,interactive
#SBATCH --time=1-21:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --job-name=SRP166333_03_BatchCorrection
#SBATCH --output=%x_%j.out  # job name and job ID in output

# Move to current working directory
cd "$(pwd)"

# Activate conda environment
source /home/sah2p/ondemand/tools/miniconda3/bin/activate scanpy

# Define paths
CONTEXT_PATH="/home/sah2p/ondemand/singlecell_data/Arabidopsis/"
INPUT_FILE="SRP166333_hvg_filtered.h5ad"
OUTPUT_FILE="SRP166333_batchCorrected.h5ad"

# Run the Python script
srun python3 03_batchCorrection.py --context_path $CONTEXT_PATH --input_file $INPUT_FILE --output_file $OUTPUT_FILE

# Print job completion message
echo "✅ Job Completed at $(date)"