#!/bin/bash

#SBATCH --account=gchopra
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --job-name=sample_2_CT_OFF%j
#SBATCH --output=logs/FAME/OFF/sample/sample_2_CT_%j_output.txt
#SBATCH --error=logs/FAME/OFF/sample/sample_2_CT_%j_err.txt

# Load Anaconda module
module load anaconda/2024.02-py311
source activate /scratch/negishi/iyer95/conda/CLAW

# Default values
INPUT_FILE="Projects/FAME/mzml_parsed/OFF/mzml_parser_1_OFF_intensity.parquet"
OUTPUT_DIR="Projects/FAME/samples/OFF/"
COLUMNS_CONFIG="core/config/columns_config.json"  # Corrected path
MAX_WORKERS=64
STD="d2-16:0"
PARENT_ION=425.40
PRODUCT_ION=183
TOLERANCE=0.3

# Run the Python script with parameters
python core/python/FAME/OFF/sample_2_CT_OFF.py \
    --input_file "$INPUT_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --columns_config "$COLUMNS_CONFIG" \
    --max_workers "$MAX_WORKERS" \
    --std "$STD" \
    --parent_ion "$PARENT_ION" \
    --product_ion "$PRODUCT_ION" \
    --tolerance "$TOLERANCE"
