#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/analysis_5_%j_output.txt
#SBATCH --error=core/backend/logs/analysis_5_%j_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=analysis_5_task_${current_date_time}

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input variables
input_file_path="df_group_OFF_4.parquet"
output_file_path="df_analysis_OFF_5.parquet"
height=1000
width=2
rel_height=0.5

# Print the current working directory to the error log
pwd >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Run the Python script with the input variables
python core/python/off/analysis_OFF_5.py "$input_file_path" "$output_file_path" "$height" "$width" "$rel_height"

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the current working directory to the error log after running the Python script
pwd >&2

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2
