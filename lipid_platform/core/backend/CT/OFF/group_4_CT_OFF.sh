#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=logs/CT/OFF/group/%A_%a_output.txt
#SBATCH --error=logs/CT/OFF/group/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --array=0-3  # Adjust this to the number of files you have

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=group_4_task_${current_date_time}_%a


# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW


# Define input and output directories
input_dir="Projects/CT/match/OFF/"
output_dir="Projects/CT/group/OFF/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Get the specific file for this array task
input_file_path=$(ls "$input_dir"*.parquet* | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

# Define STD_Only flag
STD_Only="Sample"  # You can modify this as needed or make it dynamic

# Check if input_file_path is empty
if [ -z "$input_file_path" ]; then
    echo "Error: No input file found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path with STD_Only: $STD_Only" >&2

# Run the Python script with the input file path, output directory, and STD_Only flag
python core/python/CT/OFF/group_4_CT_OFF.py \
    --input_file "$input_file_path" \
    --output_dir "$output_dir" \
    --STD_Only "$STD_Only"

# Print the current working directory to the error log after running the Python script
pwd >&2
