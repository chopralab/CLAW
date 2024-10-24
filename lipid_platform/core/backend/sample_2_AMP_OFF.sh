#!/bin/bash

#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/sample/sample_2_AMP_%j_output.txt
#SBATCH --error=core/backend/logs/sample/sample_2_AMP_%j_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --job-name=sample_extract

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Run the Python script
python core/python/sample_2_AMP_OFF.py
