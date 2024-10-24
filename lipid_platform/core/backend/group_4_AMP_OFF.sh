#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/group/group_4_AMP_%A_%a_output.txt
#SBATCH --error=core/backend/logs/group/group_4_AMP_%A_%a_err.txt
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

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input and output directories
input_dir="Projects/AMP/match/OFF/"
output_dir="Projects/AMP/group/OFF/"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Get the specific file for this array task
input_file_path=$(ls $input_dir*.parquet* | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

# Check if STD_Only is set 
STD_Only="Sample"  # You can pass this dynamically or keep it hardcoded

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path with STD_Only: $STD_Only" >&2

# Run the Python script with the input file path and STD_Only flag
python core/python/group_4_AMP_OFF.py "$input_file_path" "$STD_Only"

# Print the current working directory to the error log after running the Python script
pwd >&2
