#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=logs/CT/NP/notpossible/%A_%a_output.txt
#SBATCH --error=logs/CT/NP/notpossible/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-3  # Adjust this based on the number of input files or tasks

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=lipid_matching_${current_date_time}_%a


# Load Anaconda module
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW


# Define default directories and parameters
DEFAULT_OZON_DIR="Projects/CT/isomer_filter_6/"
DEFAULT_OZOFF_DIR="Projects/CT/analysis/OFF/notpossible/"
DEFAULT_OUTPUT_DIR="Projects/CT/notpossible_9/"
DEFAULT_RT_WINDOW=0.05

# Accept input parameters with defaults
ozon_dir=${1:-$DEFAULT_OZON_DIR}
ozoff_dir=${2:-$DEFAULT_OZOFF_DIR}
output_dir=${3:-$DEFAULT_OUTPUT_DIR}
rt_window=${4:-$DEFAULT_RT_WINDOW}

# Remove trailing slashes if they exist
ozon_dir=$(echo "$ozon_dir" | sed 's:/*$::')
ozoff_dir=$(echo "$ozoff_dir" | sed 's:/*$::')
output_dir=$(echo "$output_dir" | sed 's:/*$::')

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Print the current working directory and parameters to the error log
pwd >&2
echo "Running lipid matching with parameters:" >&2
echo "  OzON Directory: $ozon_dir" >&2
echo "  OzOFF Directory: $ozoff_dir" >&2
echo "  Output Directory: $output_dir" >&2
echo "  Retention Time Window: $rt_window" >&2

# Measure and print the time taken by the Python script
start_time=$(date +%s)

# Run the Python script with all flags
python core/python/CT/NP/notpossible_9_CT_NP.py \
    --ozon_dir "$ozon_dir" \
    --ozoff_dir "$ozoff_dir" \
    --output_dir "$output_dir" \
    --rt_window "$rt_window"

end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))

# Print the execution time to the log
echo "Script execution time: ${elapsed_time} seconds" >&2
