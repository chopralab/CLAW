#!/bin/bash

#SBATCH --account=gchopra
#SBATCH --job-name=match_3_CT_ON_%j
#SBATCH --output=logs/FAME/ON/match/%j_output.txt
#SBATCH --error=logs/FAME/ON/match/%j_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --array=0-2  # Adjust this range based on the number of input files

### ============================ ###
###       Environment Setup      ###
### ============================ ###

# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

### ============================ ###
###        Variable Setup        ###
### ============================ ###

# Define the Python script path
PYTHON_SCRIPT="core/python/FAME/ON/match_3_CT_ON.py"  # Ensure this path is correct

# Define directories and file paths
OZOFF_DIR="Projects/FAME/analysis/OFF/"  # Directory containing OzOFF parquet files
OZON_DATABASE="lipid_database/OzON_databases/OzON_Possible_Database_0.parquet"  # Path to OzON database
SAMPLE_DIR="Projects/FAME/samples/ON/"                # Directory containing sample parquet files
OUTPUT_DIR="Projects/FAME/match/ON/"                  # Directory to save output files

# Parameters for lipid matching (can be adjusted as needed)
TOLERANCE=0.3
RETENTION_TIME_WINDOW=0.5
LOG_LEVEL="INFO"

# Additional processing parameters (if applicable)
HEIGHT=500          # Example value; adjust as needed
WIDTH=2             # Example value; adjust as needed
REL_HEIGHT=0.5      # Example value; adjust as needed

# Flag to create max peaks DataFrame
MAX_PEAKS=true      # Set to true or false

### ============================ ###
###      Prepare Input Files     ###
### ============================ ###

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Gather all .parquet files in the sample directory
mapfile -t files < <(ls "${SAMPLE_DIR}"*.parquet 2>/dev/null)

# Total number of files found
num_files=${#files[@]}

# Log the number of files found
echo "Found $num_files .parquet files in input directory: $SAMPLE_DIR"

### ============================ ###
###          Processing          ###
### ============================ ###

# Check if the SLURM_ARRAY_TASK_ID is within the range of available files
if [ "$SLURM_ARRAY_TASK_ID" -lt "$num_files" ]; then
    # Select the sample file based on the task ID
    sample_path="${files[$SLURM_ARRAY_TASK_ID]}"
    sample_filename=$(basename "$sample_path")
    
    echo "Processing file [Task ID: $SLURM_ARRAY_TASK_ID]: $sample_filename"

    # Determine the max_peaks flag
    if $MAX_PEAKS; then
        MAX_PEAKS_FLAG="--max_peaks"
    else
        MAX_PEAKS_FLAG=""
    fi

    # Execute the Python script with appropriate flags
    python "$PYTHON_SCRIPT" \
        --ozoff_dir "$OZOFF_DIR" \
        --ozon_database "$OZON_DATABASE" \
        --sample_file "$sample_path" \
        --output_dir "$OUTPUT_DIR" \
        --tolerance "$TOLERANCE" \
        --retention_time_window "$RETENTION_TIME_WINDOW" \
        --log_level "$LOG_LEVEL" \
        $MAX_PEAKS_FLAG \
        --height "$HEIGHT" \
        --width "$WIDTH" \
        --rel_height "$REL_HEIGHT"

    # Check if the Python script executed successfully
    if [ $? -eq 0 ]; then
        echo "Successfully processed: $sample_filename"
    else
        echo "Error processing: $sample_filename" >&2
        exit 1
    fi
else
    echo "Warning: Task ID $SLURM_ARRAY_TASK_ID exceeds the number of files ($num_files). Skipping." >&2
fi
