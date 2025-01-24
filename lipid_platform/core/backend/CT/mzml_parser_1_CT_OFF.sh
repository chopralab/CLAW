#!/bin/bash

# SLURM configuration
#SBATCH --account=gchopra
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/mzml_parser_CT_%j_output.txt
#SBATCH --error=logs/mzml_parser_CT_%j_err.txt

# Configuration
INPUT_DIR="Projects/CT/mzml/OFF"
OUTPUT_DIR="Projects/CT/mzml_parsed/OFF"
OUTPUT_PREFIX="OFF"

# Environment setup
module load anaconda/2024.02-py311
source activate /home/iyer95/.conda/envs/CLAW

# Create log directory
mkdir -p logs

# Combine stdout and stderr
exec &> >(tee -a "logs/mzml_parser_${SLURM_JOB_ID}.log")

# Log start time and working directory
echo "Starting job at $(date)"
echo "Working directory: $(pwd)"

# Run parser with arguments
SECONDS=0
python core/python/CT/mzml_parser_1_CT_OFF.py \
    "$INPUT_DIR" \
    "$OUTPUT_DIR" \
    --prefix "$OUTPUT_PREFIX"

# Log completion time
echo "Job completed in $SECONDS seconds"