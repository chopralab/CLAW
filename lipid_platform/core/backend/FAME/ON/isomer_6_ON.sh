#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --job-name=isomer_6_FAME_ON_%A_%a
#SBATCH --output=logs/FAME/ON/isomer/%A_%a_output.txt
#SBATCH --error=logs/FAME/ON/isomer/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --array=0

# ============================
# Environment Setup
# ============================

# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# ============================
# Debugging Information
# ============================

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "SLURM Array Task ID: ${SLURM_ARRAY_TASK_ID}"

# ============================
# Variable Definitions
# ============================

# Absolute Paths for Directories
INPUT_DIR="Projects/FAME/analysis/ON/"
OFF_POSSIBLE_DIR="Projects/FAME/analysis/OFF/"
OUTPUT_DIR="Projects/FAME/isomer_filter/"

# Parameters
RETENTION_TIME_TOLERANCE=0.15

# Python Script
PYTHON_SCRIPT="core/python/FAME/ON/isomer_filter_6_ON.py"

# ============================
# Debugging: List Files Before Execution
# ============================

echo "Listing files in OzOFF directory: ${OFF_POSSIBLE_DIR}"
ls -l "${OFF_POSSIBLE_DIR}"

echo "Listing files in OzON directory: ${INPUT_DIR}"
ls -l "${INPUT_DIR}"

# ============================
# Execute Python Script
# ============================

python -u "${PYTHON_SCRIPT}" \
    --input_dir "${INPUT_DIR}" \
    --off_possible_dir "${OFF_POSSIBLE_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --retention_time_tolerance "${RETENTION_TIME_TOLERANCE}"

# ============================
# Job Completion
# ============================

echo "Job finished at: $(date)"
