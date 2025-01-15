#!/bin/bash
#SBATCH --partition=requeue
#SBATCH --time=1-21:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --job-name=05_add_literature_annotations_2024_10_29
#SBATCH --output=%x_%j.out  # job name and job ID in output

# Move to current working directory
cd "$(pwd)"

# Activate conda environment
source /cluster/pixstor/xudong-lab/sah2p/tools/miniconda3/bin scanpy

# Check if input file exists
INPUT_FILE="/cluster/pixstor/xudong-lab/sah2p/projects/sc_data/GSE131907/03_HVG_adata.h5ad"
INPUT_ANN_FILE="/cluster/pixstor/xudong-lab/sah2p/projects/sc_data/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt"
OUTPUT_DIR="/cluster/pixstor/xudong-lab/sah2p/projects/sc_data/GSE131907/"
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist."
    exit 1
fi

# Log environment information
echo "Job started on $(hostname) at $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Job Name: $SLURM_JOB_NAME"
echo "Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Run the Python script
srun python3 05_add_literature_annotations.py \
    --input_path "$INPUT_FILE" \
    --input_path_ann "$INPUT_ANN_FILE" \
    --output_path "$OUTPUT_DIR"

# Confirm completion
echo "Job completed at $(date)"