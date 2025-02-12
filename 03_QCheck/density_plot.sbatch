#!/bin/bash
#SBATCH --partition=general,requeue
#SBATCH --time=1-21:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --job-name=DenistyPlot_SRP166333
#SBATCH --output=%x_%j.out  # job name and job ID in output

# Move to current working directory
cd "$(pwd)"

# Activate conda environment
source /home/sah2p/ondemand/tools/miniconda3/bin/activate scanpy

# Log environment information
echo "Job started on $(hostname) at $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Job Name: $SLURM_JOB_NAME"
echo "Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Run Python Script
python density_plot.py
# Confirm completion
echo "Job completed at $(date)"