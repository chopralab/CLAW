#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/match_combine_3_%j_output.txt
#SBATCH --error=core/backend/logs/match_combine_3_%j_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Run the Python script
python core/python/match_combine_3.py
