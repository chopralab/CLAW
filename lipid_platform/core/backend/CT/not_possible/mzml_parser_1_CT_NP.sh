#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --job-name=mzml_1_CT_NP_%A_%a
#SBATCH --output=logs/CT/NP/mzml/%A_%a_output.txt
#SBATCH --error=logs/CT/NP/mzml/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00

# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Print the current working directory to the standard output
echo "Current working directory: $(pwd)"

# Define variables for each flag
PYTHON_SCRIPT="core/python/CT/NP/mzml_parser_1_CT_NP.py"
INPUT_FOLDER="Projects/CT/mzml/OFF"  # <-- Update this path as needed
OUTPUT_TRANSITION="Projects/CT/mzml_parsed/NP"  # <-- Define your desired output prefix
OUTPUT_OZESI="Projects/CT/mzml_parsed/NP"  # <-- Define your desired output prefix

# Ensure the output directories exist
mkdir -p "$(dirname "$OUTPUT_TRANSITION")"
mkdir -p "$(dirname "$OUTPUT_OZESI")"

# Record the start time
start_time=$(date +%s)

# Run the Python script with the defined flags
python "$PYTHON_SCRIPT" \
    --input_folder "$INPUT_FOLDER" \
    --output_transition "$OUTPUT_TRANSITION" \
    --output_ozesi "$OUTPUT_OZESI"

# Record the end time
end_time=$(date +%s)

# Calculate the elapsed time in seconds
elapsed_time=$((end_time - start_time))

# Print the elapsed time to the standard output
echo "Script execution time: ${elapsed_time} seconds"

# Print the current working directory after running the Python script
echo "Current working directory after execution: $(pwd)"
