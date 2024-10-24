#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/group/group_4_%A_%a_output.txt
#SBATCH --error=core/backend/logs/group/group_4_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --array=0-41  # Array indices for the first 5 files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=group_4_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input and output directories
input_dir="Projects/AMP/match/"
output_dir="Projects/AMP/group/"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# List the first 5 input files
input_files=($(ls $input_dir/*.parquet | head -n 5))

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path" >&2

# Run the Python script with the input file path
python core/python/group_4.py "$input_file_path"

# Print the current working directory to the error log after running the Python script
pwd >&2
