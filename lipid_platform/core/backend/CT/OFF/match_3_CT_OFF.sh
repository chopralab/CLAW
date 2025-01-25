#!/bin/bash

#SBATCH --account=gchopra
#SBATCH --output=logs/CT/OFF/match/%A_%a_output.txt
#SBATCH --error=logs/CT/OFF/match/%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-40
#SBATCH --job-name=match_3_$(date +"%Y%m%d_%H%M%S")_%a

# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Define input variables (modifiable flags)
PYTHON_SCRIPT="core/python/CT/OFF/match_3_CT_OFF.py"  # Path to the Python script
DATABASE="lipid_database/OzOFF_database/Database_OzOFF.parquet"
INPUT_DIR="Projects/CT/samples/OFF/"
OUTPUT_DIR="Projects/CT/match/OFF/"
TOLERANCE=0.3
LOG_LEVEL="INFO"

# Ensure the Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Python script $PYTHON_SCRIPT not found." >&2
    exit 1
fi

# List all Parquet files in the input directory
files=($(ls "${INPUT_DIR}"*.parquet))
num_files=${#files[@]}

# Check if the current task ID is within the range of available files
if [ "$SLURM_ARRAY_TASK_ID" -lt "$num_files" ]; then
    sample_path="${files[$SLURM_ARRAY_TASK_ID]}"
    
    # Log the current working directory
    pwd >&2
    
    # Run the Python script with the defined flags
    python "$PYTHON_SCRIPT" \
        --database "$DATABASE" \
        --sample "$sample_path" \
        --output "$OUTPUT_DIR" \
        --tolerance "$TOLERANCE" \
        --log-level "$LOG_LEVEL"
    
    # Log the current working directory after script execution
    pwd >&2
else
    echo "Task ID $SLURM_ARRAY_TASK_ID exceeds file count (${num_files})." >&2
    exit 1
fi
