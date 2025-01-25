#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --job-name=mzml_1_CT_ON_%j
#SBATCH --output=logs/CT/ON/mzml/%j_output.txt
#SBATCH --error=logs/CT/ON/mzml/%j_error.txt

# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Define input and output paths
INPUT_FOLDER="Projects/CT/mzml/ON/"
OUTPUT_FOLDER="Projects/CT/mzml_parsed/ON/"
TRANSITION_SUBDIR="transitions/"
OUTPUT_TRANSITION_FILE="mzml_transition_summed_1_CT_ON"
OUTPUT_OZESI_FILE="mzml_parser_1_CT_ON"

# Define the Python script path
PYTHON_SCRIPT="core/python/CT/ON/mzml_parser_1_CT_ON.py"

# Run the Python script
python "${PYTHON_SCRIPT}" \
    --input_folder "${INPUT_FOLDER}" \
    --output_folder "${OUTPUT_FOLDER}" \
    --transition_subdir "${TRANSITION_SUBDIR}" \
    --output_transition_file "${OUTPUT_TRANSITION_FILE}" \
    --output_ozesi_file "${OUTPUT_OZESI_FILE}"

# Print script completion message
echo "Script execution completed." >&2
