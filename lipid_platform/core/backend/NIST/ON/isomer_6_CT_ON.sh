#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --job-name=isomer_6_NIST_ON_%A_%a
#SBATCH --output=logs/NIST/ON/isomer/%A_%a_output.txt
#SBATCH --error=logs/NIST/ON/isomer/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --array=0

# ============================
# Environment Setup
# ============================
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
INPUT_DIR="Projects/NIST/analysis/ON/"
OFF_POSSIBLE_DIR="Projects/NIST/analysis/OFF/"
OUTPUT_DIR="Projects/NIST/isomer_filter/"
RETENTION_TIME_TOLERANCE=0.25

# Python Script
PYTHON_SCRIPT="core/python/NIST/ON/isomer_filter_6_CT_ON.py"

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
    --retention_time_tolerance "${RETENTION_TIME_TOLERANCE}"
    # Uncomment and adjust if you want specific files only:
    # --specific_off_files "${SPECIFIC_OFF_FILES[@]}" \
    # --specific_on_files "${SPECIFIC_ON_FILES[@]}" \

# ============================
# Job Completion
# ============================
echo "Job finished at: $(date)"
