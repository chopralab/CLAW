#!/bin/bash

# =============================================================================
# SLURM Job Script for Grouping Lipid Data
# =============================================================================

#SBATCH --account=gchopra
#SBATCH --job-name=group_4_NP_%A_%a
#SBATCH --output=logs/CT/NP/group/group_4_CT_%A_%a_output.txt
#SBATCH --error=logs/CT/NP/group/group_4_CT_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --array=0-40  # Adjust this to the number of files you have

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=group_4_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# =============================================================================
# Configuration Variables
# =============================================================================

# Define paths and script
PYTHON_SCRIPT="core/python/CT/NP/group_4_AMP_notpossible.py"

# These can be set as environment variables or passed as arguments to the script
INPUT_DIR="${INPUT_DIR:-Projects/AMP/match/OFF/notpossible/}"
OUTPUT_DIR="${OUTPUT_DIR:-Projects/AMP/group/OFF/notpossible/}"

# STD_Only flag can be set to 'STD' or 'Sample'
STD_ONLY="${STD_ONLY:-Sample}"

# =============================================================================
# Job Execution
# =============================================================================

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Get the specific file for this array task
input_file_path=$(ls "$INPUT_DIR"*.parquet* | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

# Check if input_file_path is not empty
if [ -z "$input_file_path" ]; then
    echo "No input file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path with STD_Only: $STD_ONLY" >&2

# Run the Python script with the input file path, STD_Only flag, and output directory
python "$PYTHON_SCRIPT" --input_file "$input_file_path" --std_only "$STD_ONLY" --output_dir "$OUTPUT_DIR"

# Print the current working directory to the error log after running the Python script
pwd >&2
