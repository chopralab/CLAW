#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=logs/NIST/OFF/analysis/%A_%a_output.txt
#SBATCH --error=logs/NIST/OFF/analysis/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-3  # Adjust this based on the number of files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=analysis_5_task_${current_date_time}_%a

# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Define input and output directories
INPUT_DIR="Projects/NIST/group/OFF"
OUTPUT_DIR="Projects/NIST/analysis/OFF"

# Remove trailing slashes if they exist
INPUT_DIR=$(echo "$INPUT_DIR" | sed 's:/*$::')
OUTPUT_DIR=$(echo "$OUTPUT_DIR" | sed 's:/*$::')

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# List all input files
input_files=("$INPUT_DIR"/*.parquet)

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Check if the input file exists (to avoid errors if the task index is out of range)
if [ -z "$input_file_path" ]; then
  echo "No file found for task ID $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

# Set the flags and parameters
IGNORE_COLUMNS_FLAG="ignore"  # Options: "ignore" or "keep"
HEIGHT=1000
WIDTH=2
REL_HEIGHT=0.5
OUTPUT_DIR_SCRIPT="$OUTPUT_DIR"

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path with ignore_columns_flag: $IGNORE_COLUMNS_FLAG" >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Run the Python script with the input file path and parameters
python core/python/NIST/OFF/analysis_5_CT_OFF.py \
  "$input_file_path" \
  "$HEIGHT" \
  "$WIDTH" \
  "$REL_HEIGHT" \
  "$IGNORE_COLUMNS_FLAG" \
  --output_dir "$OUTPUT_DIR_SCRIPT"

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the current working directory to the error log after running the Python script
pwd >&2

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2
