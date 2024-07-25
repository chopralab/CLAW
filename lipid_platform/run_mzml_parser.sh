#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=logs/cpu_task_output_%j.txt
#SBATCH --error=logs/cpu_task_error_%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00

# Generate a timestamp-based job name
current_date_time=$(date +"%Y%m%d_%H%M%S")
#SBATCH --job-name=cpu_task_${current_date_time}

# Load Anaconda module
module load anaconda/2024.02-py311

# Activate the conda environment
source activate /home/iyer95/.conda/envs/CLAW

# Define the path to the mzML data
mzml_data='./Projects/AMP/mzml/AMP_ON/'

# Run the Python script with the mzML data path as an argument
python parse_mzml_files.py $mzml_data
