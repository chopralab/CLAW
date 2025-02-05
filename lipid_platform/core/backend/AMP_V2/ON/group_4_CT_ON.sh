#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --job-name=group_4_CT_ON_%A_%a
#SBATCH --output=logs/AMP_V2/ON/group/%A_%a_output.txt
#SBATCH --error=logs/AMP_V2/ON/group/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --array=0-2  # Ensure this matches the number of input files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=group_4_task_${current_date_time}_%a


# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Define variables
PYTHON_SCRIPT="core/python/AMP_V2/ON/group_4_CT_ON.py"
INPUT_DIR="Projects/AMP_V2/match/ON"
OUTPUT_DIR="Projects/AMP_V2/group/ON"

# Remove trailing slashes
INPUT_DIR="${INPUT_DIR%/}"
OUTPUT_DIR="${OUTPUT_DIR%/}"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# List all input files and store in an array
mapfile -t input_files < <(ls "$INPUT_DIR"/*.parquet)

# Check if SLURM_ARRAY_TASK_ID is within the range
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#input_files[@]}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID is out of bounds." >&2
    exit 1
fi

# Get the specific file for this array task
input_file_path="${input_files[$SLURM_ARRAY_TASK_ID]}"

# Log the processing details
echo "Processing file: $input_file_path" >&2
echo "Output directory: $OUTPUT_DIR" >&2
echo "Python script: $PYTHON_SCRIPT" >&2

# Run the Python script with the input file and output directory
python "$PYTHON_SCRIPT" --input_file "$input_file_path" --output_dir "$OUTPUT_DIR" --log_level INFO

# Log completion
echo "Completed processing file: $input_file_path" >&2
