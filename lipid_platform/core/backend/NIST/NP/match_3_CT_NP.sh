#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --job-name=match_NP_%A_%a
#SBATCH --output=logs/NIST/NP/match/match_3_CT_%A_%a_output.txt
#SBATCH --error=logs/NIST/NP/match/match_3_CT_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --array=0-5  # Array indices for 41 jobs

# Load Anaconda module
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Define input variables
MRM_DATABASE_PATH="lipid_database/OzOFF_database/Database_OzOFF.parquet"
SAMPLE_INPUT_DIR="Projects/NIST/samples/OFF/"
OUTPUT_DIR="Projects/NIST/match/NP/"
PYTHON_SCRIPT="core/python/NIST/NP/match_3_CT_NP.py"  # Updated to match the Python script name

# List all match files in the input directory
sample_files=($(ls "${SAMPLE_INPUT_DIR}"*.parquet))
total_files=${#sample_files[@]}

# Check if the current task ID is within the range of available files
if [ ${SLURM_ARRAY_TASK_ID} -lt ${total_files} ]; then
    # Get the match file corresponding to the current task ID
    sample_path="${sample_files[${SLURM_ARRAY_TASK_ID}]}"
    sample_name=$(basename "${sample_path}" .parquet)

    # Log the sample being processed
    echo "Processing sample: ${sample_name}" >&2

    # Run the Python script with the specified flags
    python "${PYTHON_SCRIPT}" \
        --mrm_database "${MRM_DATABASE_PATH}" \
        --sample_file "${sample_path}" \
        --output_dir "${OUTPUT_DIR}" \
        --tolerance 0.3  # You can modify or parameterize this as needed

    # Log completion of the sample processing
    echo "Completed processing sample: ${sample_name}" >&2
else
    echo "Task ID ${SLURM_ARRAY_TASK_ID} is out of range for available files (${total_files} files)." >&2
fi
