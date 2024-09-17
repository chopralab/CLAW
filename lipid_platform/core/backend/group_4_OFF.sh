#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/group/group_4_%A_%a_output.txt
#SBATCH --error=core/backend/logs/group/group_4_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --array=0-1  # Array indices for the first 5 files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=group_4_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input and output directories
input_dir="Projects/STD/match/OFF/"
output_dir="Projects/STD/group/OFF/"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# List the first 5 input files
input_files=($(ls $input_dir*.parquet | head -n 5))

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Check if STD_Only is set
STD_Only="STD"  # You can pass this dynamically or keep it hardcoded

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path with STD_Only: $STD_Only" >&2

# Run the Python script with the input file path and STD_Only flag
python core/python/group_4_OFF.py "$input_file_path" "$STD_Only"

# Print the current working directory to the error log after running the Python script
pwd >&2
