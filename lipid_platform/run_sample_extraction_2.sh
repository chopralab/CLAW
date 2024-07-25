#!/bin/bash

#SBATCH --account=gchopra
#SBATCH --output=logs/sample_extract_%j_output.txt
#SBATCH --error=logs/sample_extract_%j_error.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=sample_extract

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Run the Python script
python sample_id_extract_2.py
