import argparse
from dca.api import dca
import os

def main(input_file, output_file, epochs, hidden_size):
    # Ensure input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")
    
    print(f"Running DCA with GPU support on file: {input_file}")
    print(f"Output will be saved to: {output_file}")
    
    # Run DCA
    dca(
        input_file,
        output_file,
        epochs=epochs,  # Number of epochs for training
        hidden_size=hidden_size,  # Size of hidden layers
        gpu=True  # Enable GPU
    )
    print("DCA processing completed successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DCA on single-cell RNA-seq data using GPU.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input .h5ad file.")
    parser.add_argument("--output_file", type=str, required=True, help="Path to save the denoised .h5ad file.")
    parser.add_argument("--epochs", type=int, default=300, help="Number of training epochs (default: 300).")
    parser.add_argument("--hidden_size", type=int, default=64, help="Hidden layer size (default: 64).")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.epochs, args.hidden_size)