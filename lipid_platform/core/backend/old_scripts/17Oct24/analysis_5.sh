#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/analysis/analysis_5_%A_%a_output.txt
#SBATCH --error=core/backend/logs/analysis/analysis_5_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-41  # Adjust the array size based on the number of files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=analysis_5_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input and output directories
input_dir="Projects/AMP/group/"
output_dir="Projects/AMP/analysis/"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# List the first 5 input files
input_files=($(ls $input_dir/*.parquet | head -n 5))

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path" >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Run the Python script with the input file path
python core/python/analysis_5.py "$input_file_path" 1000 2 0.5

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the current working directory to the error log after running the Python script
pwd >&2

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2
