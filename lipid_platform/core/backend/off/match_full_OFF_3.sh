#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/match_combine_3_%j_output.txt
#SBATCH --error=core/backend/logs/match_combine_3_%j_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Run the Python script
python core/python/match_full_OFF_3.py lipid_database/AMP_Database_OFF.parquet df_sample_OFF_2.parquet
