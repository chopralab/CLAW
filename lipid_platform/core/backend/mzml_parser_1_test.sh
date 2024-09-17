#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/mzml_parser/mzml_parser_1_%j_output.txt
#SBATCH --error=core/backend/logs/mzml_parser/mzml_parser_1_%j_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=cpu_task_${current_date_time}

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Print the current working directory to the error log
pwd >&2

# Define the path to the mzML data
mzml_data='/home/iyer95/CLAW/CLAW/lipid_platform/Projects/STD/mzml/ON/'
mzml_data_OLD='Projects/STD/mzml/ON/'

# Print the current working directory to the error log before running the Python script
pwd >&2

# Record the start time
start_time=$(date +%s)

# Run the Python script with the mzML data path as an argument
python core/python/mzml_parser_1_test.py $mzml_data

# Record the end time
end_time=$(date +%s)

# Calculate the elapsed time in seconds
elapsed_time=$((end_time - start_time))

# Print the elapsed time to the error log
echo "Script execution time: ${elapsed_time} seconds" >&2

# Print the current working directory to the error log after running the Python script
pwd >&2
