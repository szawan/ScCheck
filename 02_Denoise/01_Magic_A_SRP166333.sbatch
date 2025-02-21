#!/bin/bash
#SBATCH --partition=general,requeue
#SBATCH --time=1-21:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --job-name=01_magic_processing_SRP173393
#SBATCH --output=%x_%j.out  # job name and job ID in output

# Move to current working directory
cd "$(pwd)"

# Activate conda environment
source /home/sah2p/ondemand/tools/miniconda3/bin/activate scanpy

# Input and output paths
INPUT_FILE="/home/sah2p/ondemand/singlecell_data/Arabidopsis/SRP166333_batchCorrected.h5ad"
OUTPUT_FILE="/home/sah2p/ondemand/singlecell_data/Arabidopsis/SRP166333_magic_denoised.h5ad"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist."
    exit 1
fi

# Log environment information
echo "Job started on $(hostname) at $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Job Name: $SLURM_JOB_NAME"
echo "Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Run the Python script for MAGIC
srun python3 01_Magic.py \
    --input_path "$INPUT_FILE" \
    --output_path "$OUTPUT_FILE"

# Check if MAGIC output file was created
if [ -f "$OUTPUT_FILE" ]; then
    echo "MAGIC processing completed successfully. Output file: $OUTPUT_FILE"
else
    echo "Error: MAGIC output file not found."
    exit 1
fi

# Confirm completion
echo "Job completed at $(date)"