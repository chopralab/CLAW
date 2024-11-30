#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/analysis/analysis_5_STDon_%A_%a_output.txt
#SBATCH --error=core/backend/logs/analysis/analysis_5_STDon_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --array=0 # Adjust the array size based on the number of files

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=analysis_5_task_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /scratch/negishi/iyer95/conda/CLAW

# Define input and output directories
input_dir="Projects/STD/group/ON/"
output_dir="Projects/STD/analysis/ON/"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# List all input files in the input directory
input_files=($(ls $input_dir/*.parquet))

# Get the specific file for this array task
input_file_path=${input_files[$SLURM_ARRAY_TASK_ID]}

# Print the current working directory and input file to the error log
pwd >&2
echo "Processing file: $input_file_path" >&2

# Define parameters for the Python script
height=500
width=2
rel_height=0.5
noise_start=25  # Adjust the noise_start as needed
noise_end=26   # Adjust the noise_end as needed

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Run the Python script with the input file path and additional parameters
python core/python/STD/analysis_5_ON.py "$input_file_path" $height $width $rel_height $noise_start $noise_end

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the current working directory to the error log after running the Python script
pwd >&2

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2
