#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=cpu_task_${current_date_time}

# Define log directories
OUTPUT_LOG_DIR="core/backend/logs/mzml_parser/mzml_parser_1_CT_ON"
mkdir -p "${OUTPUT_LOG_DIR}/output" "${OUTPUT_LOG_DIR}/error"

# SBATCH output and error paths
#SBATCH --output=${OUTPUT_LOG_DIR}/output/mzml_parser_1_%A_%a_output.txt
#SBATCH --error=${OUTPUT_LOG_DIR}/error/mzml_parser_1_%A_%a_err.txt

# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Print the current working directory to the error log
pwd >&2

# Define input and output paths
INPUT_FOLDER="Projects/CT/mzml/ON/"
OUTPUT_TRANSITION="Projects/CT/mzml_parsed/ON/sum/df_transition_summed_1_CT_ON"
OUTPUT_OZESI="Projects/CT/mzml_parsed/ON/df_mzml_parser_1_CT_ON"

# Define the Python script path
PYTHON_SCRIPT="core/python/CT/ON/mzml_parser_1_CT_ON.py"

# Record the start time
start_time=$(date +%s)

# Run the Python script with the defined arguments
python "${PYTHON_SCRIPT}" \
    --input_folder "${INPUT_FOLDER}" \
    --output_transition "${OUTPUT_TRANSITION}" \
    --output_ozesi "${OUTPUT_OZESI}"

# Record the end time
end_time=$(date +%s)

# Calculate the elapsed time in seconds
elapsed_time=$((end_time - start_time))

# Print the elapsed time to the error log
echo "Script execution time: ${elapsed_time} seconds" >&2

# Print the current working directory to the error log after running the Python script
pwd >&2
