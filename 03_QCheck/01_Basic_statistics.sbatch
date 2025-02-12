#!/bin/bash
#SBATCH --partition=general,requeue,logical_cpu,logical_cpu2,interactive
#SBATCH --time=1-21:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --job-name=01_Basic_statistics_SRP166333
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
python 01_Basic_statistics.py
# Confirm completion
echo "Job completed at $(date)"