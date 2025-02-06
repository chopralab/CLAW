#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=logs/FAME/OFF/off_possible/%A_%a_output.txt
#SBATCH --error=logs/FAME/OFF/off_possible/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-1  # Adjust this based on the number of files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=off_possible_6_task_${current_date_time}_%a

# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Define input and output directories
# The input directory holds the base Parquet files and the off_possible output will be stored in a subdirectory.
INPUT_DIR="Projects/FAME/analysis/OFF"
OUTPUT_DIR="Projects/FAME/analysis/OFF/off_possible"

# Define the Python script to run
PYTHON_SCRIPT="core/python/FAME/OFF/off_possible_6.py"

# Define additional flags/parameters to pass to the Python script
HOW_MANY=2
THRESHOLD=10000

# Remove trailing slashes if they exist
INPUT_DIR=$(echo "$INPUT_DIR" | sed 's:/*$::')
OUTPUT_DIR=$(echo "$OUTPUT_DIR" | sed 's:/*$::')

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# List all input files (assumed to be .parquet files) in the INPUT_DIR
input_files=("$INPUT_DIR"/*.parquet)

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Check if the input file exists (to avoid errors if the task index is out of range)
if [ -z "$input_file_path" ]; then
  echo "No file found for task ID $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

# Construct the output file path based on the input file's basename
output_file_path="$OUTPUT_DIR/$(basename "$input_file_path")"

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path" >&2
echo "Output will be saved to: $output_file_path" >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Run the off_possible_6.py Python script with all flags and parameters
python "$PYTHON_SCRIPT" \
    "$input_file_path" \
    "$output_file_path" \
    --how_many "$HOW_MANY" \
    --threshold "$THRESHOLD"

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the current working directory and execution time to the log
pwd >&2
echo "Script execution time: ${elapsed_time} seconds" >&2
