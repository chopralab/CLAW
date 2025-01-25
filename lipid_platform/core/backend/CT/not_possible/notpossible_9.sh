#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/analysis/notpossible/analysis_AMP_notpossible_%A_%a_output.txt
#SBATCH --error=core/backend/logs/analysis/notpossible/analysis_AMP_notpossible_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-40  # Adjust this based on the number of input files or tasks

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=lipid_matching_${current_date_time}_%a

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define directories
ozon_dir="Projects/AMP/isomer_filter_6/"
ozoff_dir="Projects/AMP/analysis/OFF/notpossible/"
output_dir="Projects/AMP/notpossible_9/"

# Remove trailing slashes if they exist
ozon_dir=$(echo $ozon_dir | sed 's:/*$::')
ozoff_dir=$(echo $ozoff_dir | sed 's:/*$::')
output_dir=$(echo $output_dir | sed 's:/*$::')

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Accept retention time window as input (default is 0.1 if not provided)
rt_window=${1:-0.1}

# Print the current working directory and parameters to the error log
pwd >&2
echo "Running lipid matching with parameters:" >&2
echo "  OzON Directory: $ozon_dir" >&2
echo "  OzOFF Directory: $ozoff_dir" >&2
echo "  Output Directory: $output_dir" >&2
echo "  Retention Time Window: $rt_window" >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Run the Python script
python core/python/not_possible/notpossible_9_AMP.py --ozon_dir "$ozon_dir" --ozoff_dir "$ozoff_dir" --output_dir "$output_dir" --rt_window $rt_window

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2
