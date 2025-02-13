#!/bin/bash
#SBATCH --partition=general     # use partition Gpu
#SBATCH --time=1-00:00:00             # Max runtime (1 day)
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=8            # Number of CPUs per task
#SBATCH --mem=48G                     # Memory (64 GB)
#SBATCH --ntasks-per-node 2
#SBATCH --job-name=dca_gpu_denoising  # Job name
#SBATCH --output=%x_%j.out            # Output file (%x=job-name, %j=job-id)

# Move to the working directory where the script is submitted
cd "$(pwd)"
module load miniconda3/4.10.3_gcc_9.5.0
module load cuda/11.8.0_gcc_9.5.0

# Activate conda environment
source /home/sah2p/ondemand/tools/miniconda3/bin/activate scanpy

# Input and output paths
INPUT_FILE="/home/sah2p/ondemand/singlecell_data/SRP267870.h5ad"
OUTPUT_FILE="/home/sah2p/ondemand/singlecell_data"

# Run the Python script
echo "Starting DCA denoising with GPU at $(date)"
srun python3 02_dca.py --input_file "$INPUT_FILE" --output_file "$OUTPUT_FILE" --epochs 300 --hidden_size 64

# Check if the output file was created
if [ -f "$OUTPUT_FILE" ]; then
    echo "DCA processing completed successfully. Output saved at $OUTPUT_FILE."
else
    echo "Error: DCA processing failed. Output file not found."
    exit 1
fi

echo "Job finished at $(date)"