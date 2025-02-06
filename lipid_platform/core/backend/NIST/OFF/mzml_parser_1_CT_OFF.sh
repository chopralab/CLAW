#!/bin/bash

# SLURM configuration
#SBATCH --account=gchopra
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/NIST/OFF/mzml_parser_CT_%j_output.txt
#SBATCH --error=logs/NIST/OFF/mzml_parser_CT_%j_err.txt

# Configuration
INPUT_DIR="Projects/NIST/mzml/OFF"
OUTPUT_DIR="Projects/NIST/mzml_parsed/OFF"
OUTPUT_PREFIX="mzml_parser_1_OFF"
PYTHON_SCRIPT="core/python/NIST/OFF/mzml_parser_1_CT_OFF.py"

# Environment setup
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Create log directory
mkdir -p logs

# Log start time and working directory
echo "Starting job at $(date)"
echo "Working directory: $(pwd)"

# Run parser with arguments
SECONDS=0
python "$PYTHON_SCRIPT" \
    "$INPUT_DIR" \
    "$OUTPUT_DIR" \
    --prefix "$OUTPUT_PREFIX"

# Log completion time
echo "Job completed in $SECONDS seconds"
