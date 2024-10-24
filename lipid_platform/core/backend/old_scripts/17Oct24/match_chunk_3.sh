#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/match_chunk_3_%A_%a_output.txt
#SBATCH --error=core/backend/logs/match_chunk_3_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-13  # Create 14 jobs in the array

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=match_3_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input variables
OzOFF_database_path="Projects/AMP/results/csv_data/OzOFF_Possible_species.parquet"
OzON_database_path="lipid_database/OzON_databases/OzON_Possible_Database_0.parquet"
df_sample_2_path="df_sample_2.parquet"

# Calculate the chunk size and indices for this job
CHUNK_SIZE=100
START_INDEX=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
END_INDEX=$((START_INDEX + CHUNK_SIZE))

# Print the current working directory to the error log
pwd >&2

# Run the Python script with the input variables and indices for this chunk
python core/python/match_chunk_3.py "$OzOFF_database_path" "$OzON_database_path" "$df_sample_2_path" $START_INDEX $END_INDEX

# Print the current working directory to the error log after running the Python script
pwd >&2
