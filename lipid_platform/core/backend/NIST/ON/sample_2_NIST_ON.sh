#!/bin/bash

#SBATCH --account=gchopra
#SBATCH --job-name=sample_2_%j
#SBATCH --output=logs/NIST/ON/sample/%j_output.txt
#SBATCH --error=logs/NIST/ON/sample/%j_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00

# Load Anaconda module and activate the environment
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Define variables
STD="no"  # Change to 'yes' if using STD
PYTHON_SCRIPT="core/python/NIST/ON/sample_2_NIST_ON.py"
INPUT_PARQUET="Projects/NIST/mzml_parsed/ON/mzml_parser_1_CT_ON.parquet"
OUTPUT_DIR="Projects/NIST/samples/ON/"

# Define ion parameters
PARENT_ION=425.40
PRODUCT_ION=183
TOLERANCE=0.3

# Execute the Python script with arguments
python $PYTHON_SCRIPT \
    --std $STD \
    --input_parquet $INPUT_PARQUET \
    --output_dir $OUTPUT_DIR \
    --parent_ion $PARENT_ION \
    --product_ion $PRODUCT_ION \
    --tolerance $TOLERANCE
