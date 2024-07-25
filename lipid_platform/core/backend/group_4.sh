#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/group_4_%j_ouptut.txt
#SBATCH --error=core/backend/logs/group_4_%j_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=group_4_task_${current_date_time}

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input variables
input_file_path="df_match_3.parquet"

# Print the current working directory to the error log
pwd >&2

# Run the Python script with the input variables
python core/python/group_4.py "$input_file_path"

# Print the current working directory to the error log after running the Python script
pwd >&2
