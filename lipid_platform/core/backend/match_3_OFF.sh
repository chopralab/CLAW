#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/match/match_3_%A_%a_output.txt
#SBATCH --error=core/backend/logs/match/match_3_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-1  # Array indices for 100 jobs, with a check to only process 41 files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=match_3_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input variables
OzOFF_database_path="lipid_database/OzOFF_database/Database_OzOFF.parquet"
input_dir="Projects/STD/samples/OFF/"
output_dir="Projects/STD/match/OFF/"

# List all files in the input directory and get the file corresponding to the SLURM_ARRAY_TASK_ID
files=($(ls $input_dir*.parquet))
num_files=${#files[@]}

# Check if the current task ID is within the range of available files
if [ $SLURM_ARRAY_TASK_ID -lt ${#files[@]} ]; then
    sample_path=${files[$SLURM_ARRAY_TASK_ID]}

    # Print the current working directory to the error log
    pwd >&2

    # Run the Python script with the input variables
    python core/python/match_3_OFF.py "$OzOFF_database_path" "$sample_path" "$output_dir"

    # Print the current working directory to the error log after running the Python script
    pwd >&2
else
    echo "Task ID $SLURM_ARRAY_TASK_ID is out of range for available files." >&2
fi
