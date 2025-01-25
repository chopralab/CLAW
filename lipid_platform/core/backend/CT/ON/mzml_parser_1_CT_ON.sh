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

# Redirect all output to both stdout and stderr
exec > >(tee -a core/backend/logs/mzml_parser/mzml_parser_1_${SLURM_JOB_ID}_output.txt)
exec 2> >(tee -a core/backend/logs/mzml_parser/mzml_parser_1_${SLURM_JOB_ID}_err.txt >&2)

# Print the current working directory to the error log
pwd >&2

# Define variables for each flag
PYTHON_SCRIPT="core/python/mzml_parser_1_CT_NP.py"
INPUT_FOLDER="/path/to/your/mzml/data"  # Update this path as needed
OUTPUT_TRANSITION="core/output/transition_summed"  # Define your desired output prefix
OUTPUT_OZESI="core/output/ozesi_data"  # Define your desired output prefix

# Optionally, you can make the paths dynamic or pass them as environment variables

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

# Print the elapsed time to the error log
echo "Script execution time: ${elapsed_time} seconds" >&2

# Print the current working directory to the error log after running the Python script
pwd >&2
