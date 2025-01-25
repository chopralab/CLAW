#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/fame_filter/notpossible/fame_filter_6_notpossible_%A_%a_output.txt
#SBATCH --error=core/backend/logs/fame_filter/notpossible/fame_filter_6_notpossible_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --array=0  # Adjust this based on the number of files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=fame_filter_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define input and output directories
input_dir="Projects/AMP/notpossible/rt_adjustment_6/"
output_dir="Projects/AMP/fame_filter_7/notpossible_nov6/"

# Remove trailing slashes if they exist
input_dir=$(echo $input_dir | sed 's:/*$::')
output_dir=$(echo $output_dir | sed 's:/*$::')

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# List all input files
input_files=($(ls $input_dir/*.parquet))

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Check if the input file exists (to avoid errors if the task index is out of range)
if [ -z "$input_file_path" ]; then
  echo "No file found for task ID $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

# Load the fame_std parquet file
fame_std="Projects/STD/off_possible/FAME_off_possible_top2.parquet"

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path with fame_std: $fame_std" >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Run the Python script with the input file path and fame_std as arguments
python core/python/not_possible/FAME_filter_6_notpossible.py "$input_file_path" "$fame_std" "$output_dir"

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the current working directory to the error log after running the Python script
pwd >&2

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2
