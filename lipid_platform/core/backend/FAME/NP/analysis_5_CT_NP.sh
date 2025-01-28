#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=logs/CT/NP/analysis/%A_%a_output.txt
#SBATCH --error=logs/CT/NP/analysis/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-3  # Adjust this based on the number of files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=analysis_5_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW


# Define configurable variables
PYTHON_SCRIPT="core/python/CT/NP/analysis_5_CT_NP.py"

INPUT_DIR="Projects/CT/group/NP/"
OUTPUT_DIR="Projects/CT/analysis/NP/"
HEIGHT=500
WIDTH=2
REL_HEIGHT=0.5
IGNORE_COLUMNS_FLAG=true  # Set to true or false as needed
MAX_PEAKS_FLAG=false      # Set to true or false as needed

# Remove trailing slashes if they exist
INPUT_DIR=$(echo "$INPUT_DIR" | sed 's:/*$::')
OUTPUT_DIR=$(echo "$OUTPUT_DIR" | sed 's:/*$::')

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# List all input files
input_files=($(ls "$INPUT_DIR"/*.parquet))

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Check if the input file exists (to avoid errors if the task index is out of range)
if [ -z "$input_file_path" ]; then
  echo "No file found for task ID $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path" >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Determine the flags for ignore_columns and max_peaks
IGNORE_COLUMNS_ARG=""
if [ "$IGNORE_COLUMNS_FLAG" = true ]; then
  IGNORE_COLUMNS_ARG="--ignore_columns"
fi

MAX_PEAKS_ARG=""
if [ "$MAX_PEAKS_FLAG" = true ]; then
  MAX_PEAKS_ARG="--max_peaks"
fi

# Run the Python script with the specified flags
python "$PYTHON_SCRIPT" \
  --input_file "$input_file_path" \
  --output_dir "$OUTPUT_DIR" \
  --height "$HEIGHT" \
  --width "$WIDTH" \
  --rel_height "$REL_HEIGHT" \
  $IGNORE_COLUMNS_ARG \
  $MAX_PEAKS_ARG

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the current working directory to the error log after running the Python script
pwd >&2

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2
