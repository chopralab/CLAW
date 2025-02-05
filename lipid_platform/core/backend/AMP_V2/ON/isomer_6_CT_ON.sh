#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --job-name=isomer_6_CT_ON_%A_%a
#SBATCH --output=logs/AMP_V2/ON/isomer/%A_%a_output.txt
#SBATCH --error=logs/AMP_V2/ON/isomer/%A_%a_err.txt
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

# Directories
INPUT_DIR="/scratch/negishi/iyer95/Projects/AMP_V2/analysis/ON/"
OFF_POSSIBLE_DIR="/scratch/negishi/iyer95/Projects/AMP_V2/analysis/OFF/off_possible/"
OUTPUT_DIR="/scratch/negishi/iyer95/Projects/AMP_V2/isomer_filter_6/"

# Parameters
RETENTION_TIME_TOLERANCE=0.15
ISOMER_FILTER_OUTPUT="isomer_filter_output"

# Python Script
PYTHON_SCRIPT="core/python/AMP_V2/ON/isomer_filter_6_CT_ON.py"

# Optional: Specify specific files (Uncomment and modify if needed)
# SPECIFIC_OFF_FILES=("file1.parquet" "file2.parquet")
# SPECIFIC_ON_FILES=("fileA.parquet" "fileB.parquet")

# ============================
# Execute Python Script
# ============================

python -u "${PYTHON_SCRIPT}" \
    --input_dir "${INPUT_DIR}" \
    --off_possible_dir "${OFF_POSSIBLE_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --retention_time_tolerance "${RETENTION_TIME_TOLERANCE}" \
    --isomer_filter_output "${ISOMER_FILTER_OUTPUT}" \
    # Uncomment the following lines to include specific files
    # --specific_off_files "${SPECIFIC_OFF_FILES[@]}" \
    # --specific_on_files "${SPECIFIC_ON_FILES[@]}" \
    #> core/backend/logs/isomer/${SLURM_ARRAY_TASK_ID}_output.log 2>&1

# ============================
# Job Completion
# ============================

echo "Job finished at: $(date)"
