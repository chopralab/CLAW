#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=logs/CT/NP/FAME_filter/%A_%a_output.txt
#SBATCH --error=logs/CT/NP/FAME_filter/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-99  # Adjust this range based on the number of files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=fame_filter_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW


# Define variables for input, output, and fame_std paths
INPUT_DIR="Projects/AMP/notpossible/rt_adjustment_6/"
OUTPUT_DIR="Projects/AMP/fame_filter_7/notpossible_nov6/"
FAME_STD="Projects/STD/off_possible/FAME_off_possible_top2.parquet"
FAME_RT_WINDOW=0.5  # Adjust as needed

# Remove trailing slashes if they exist
INPUT_DIR="${INPUT_DIR%/}"
OUTPUT_DIR="${OUTPUT_DIR%/}"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# List all input parquet files
mapfile -t input_files < <(ls "$INPUT_DIR"/*.parquet 2>/dev/null)

# Check if there are any parquet files
if [ ${#input_files[@]} -eq 0 ]; then
  echo "No parquet files found in input directory: $INPUT_DIR" >&2
  exit 1
fi

# Ensure the array index is within the range of input files
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#input_files[@]}" ]; then
  echo "Array task ID $SLURM_ARRAY_TASK_ID is out of range. Only ${#input_files[@]} files available." >&2
  exit 1
fi

# Get the specific file for this array task
INPUT_FILE_PATH="${input_files[$SLURM_ARRAY_TASK_ID]}"

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $INPUT_FILE_PATH with fame_std: $FAME_STD" >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Run the Python script with the input parameters
python core/python/FAME_filter_CT_NP.py \
  --input_dir "$(dirname "$INPUT_FILE_PATH")" \
  --fame_std "$FAME_STD" \
  --output_dir "$OUTPUT_DIR" \
  --fame_rt_window "$FAME_RT_WINDOW"

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the current working directory to the error log after running the Python script
pwd >&2

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2
