#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/group/group_4_AMPon_%A_%a_output.txt
#SBATCH --error=core/backend/logs/group/group_4_AMPon_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --array=0-40 # make sure its the same as the match_3_AMP.sh array

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=group_4_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input and output directories
input_dir="Projects/AMP/match/ON/"
output_dir="Projects/AMP/group/ON/"

# Remove trailing slashes if they exist
input_dir=$(echo $input_dir | sed 's:/*$::')
output_dir=$(echo $output_dir | sed 's:/*$::')

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# List all input files in the input directory
input_files=($(ls $input_dir/*.parquet))

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path" >&2

# Run the Python script with the input file path
python core/python/group_4_AMP.py "$input_file_path"

# Print the current working directory to the error log after running the Python script
pwd >&2
