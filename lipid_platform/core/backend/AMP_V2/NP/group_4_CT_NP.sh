#!/bin/bash

#SBATCH --account=gchopra
#SBATCH --job-name=group_4_NP_%A_%a
#SBATCH --output=logs/CT/NP/group/group_4_CT_%A_%a_output.txt
#SBATCH --error=logs/CT/NP/group/group_4_CT_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --array=0-3  # Adjust this to the number of files you have

module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Define hardcoded paths
INPUT_DIR="Projects/CT/match/NP/"
OUTPUT_DIR="Projects/CT/group/NP/"
STD_ONLY="Sample"
PYTHON_SCRIPT="core/python/CT/NP/group_4_CT_NP.py"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Get the specific file for this array task
input_file_path=$(ls "$INPUT_DIR"*.parquet* | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

if [ -z "$input_file_path" ]; then
    echo "No input file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "Processing file: $input_file_path with STD_Only: $STD_ONLY" >&2

python "$PYTHON_SCRIPT" --input_file "$input_file_path" --std_only "$STD_ONLY" --output_dir "$OUTPUT_DIR"
