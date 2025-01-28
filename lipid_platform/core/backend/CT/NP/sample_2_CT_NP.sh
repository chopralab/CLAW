#!/bin/bash

#SBATCH --account=gchopra
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --job-name=sample_extract
#SBATCH --output=logs/CT/OFF/sample/sample_2_CT_%j_output.txt
#SBATCH --error=logs/CT/OFF/sample/sample_2_CT_%j_err.txt

# Load Anaconda module
module load anaconda/2024.02-py311
source activate  /scratch/negishi/iyer95/conda/CLAW

# Default values
INPUT_FILE="Projects/CT/mzml_parsed/OFF/df_mzml_parser_1_OFF.parquet"
OUTPUT_DIR="Projects/CT/samples/OFF/"
COLUMNS_CONFIG="config/columns_config.json"
MAX_WORKERS=64
STD="d2-16:0"
PARENT_ION=425.40
PRODUCT_ION=183
TOLERANCE=0.3

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_file) INPUT_FILE="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        --columns_config) COLUMNS_CONFIG="$2"; shift ;;
        --max_workers) MAX_WORKERS="$2"; shift ;;
        --std) STD="$2"; shift ;;
        --parent_ion) PARENT_ION="$2"; shift ;;
        --product_ion) PRODUCT_ION="$2"; shift ;;
        --tolerance) TOLERANCE="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Run the Python script with parameters
python core/python/CT/NP/sample_2_CT_NP.py \
    --input_file "$INPUT_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --columns_config "$COLUMNS_CONFIG" \
    --max_workers "$MAX_WORKERS" \
    --std "$STD" \
    --parent_ion "$PARENT_ION" \
    --product_ion "$PRODUCT_ION" \
    --tolerance "$TOLERANCE"