#!/bin/bash

#SBATCH --account=gchopra
#SBATCH --job-name=analysis_5_CT_ON_%A_%a
#SBATCH --output=logs/AMP_V2/ON/analysis/%A_%a_output.txt
#SBATCH --error=logs/AMP_V2/ON/analysis/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-2  # Adjust the array size based on the number of files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=analysis_5_task_${current_date_time}_%a

# Load Anaconda module
# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# ============================
# Configuration Variables
# ============================

# Define the Python script path
PYTHON_SCRIPT="core/python/AMP_V2/ON/analysis_5_CT_ON.py"

# Define input and output directories
INPUT_DIR="Projects/AMP_V2/group/ON/"
OUTPUT_DIR="Projects/AMP_V2/analysis/ON/"

# Define peak detection parameters
HEIGHT=500          # Example value; adjust as needed
WIDTH=2             # Example value; adjust as needed
REL_HEIGHT=0.5      # Example value; adjust as needed

# Flag to create max peaks DataFrame
MAX_PEAKS=true      # Set to true or false

# ============================
# End of Configuration
# ============================

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Verify that the Python script exists
if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "Python script $PYTHON_SCRIPT not found!" >&2
    exit 1
fi

# List all input files in the input directory
input_files=("$INPUT_DIR"/*.parquet)

# Check if there are any input files
if [[ ${#input_files[@]} -eq 0 ]]; then
    echo "No Parquet files found in $INPUT_DIR" >&2
    exit 1
fi

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Verify that the input file exists
if [[ ! -f "$input_file_path" ]]; then
    echo "Input file $input_file_path does not exist!" >&2
    exit 1
fi

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path" >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Determine the max_peaks flag
if $MAX_PEAKS; then
    MAX_PEAKS_FLAG="--max_peaks"
else
    MAX_PEAKS_FLAG=""
fi

# Run the Python script with the specified parameters
python "$PYTHON_SCRIPT" \
    --input_file "$input_file_path" \
    --output_dir "$OUTPUT_DIR" \
    --height "$HEIGHT" \
    --width "$WIDTH" \
    --rel_height "$REL_HEIGHT" \
    $MAX_PEAKS_FLAG

# Capture the exit status of the Python script
exit_status=$?

# End timing
end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the current working directory to the error log after running the Python script
pwd >&2

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2

# Exit with the Python script's exit status
exit $exit_status
