#!/bin/bash
#SBATCH --account=gchopra
#SBATCH --output=core/backend/logs/isomer/isomer_filter_6_%A_%a_output.txt
#SBATCH --error=core/backend/logs/isomer/isomer_filter_6_%A_%a_err.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --array=0

module load anaconda/2024.02-py311
source activate /home/iyer95/.conda/envs/CLAW

# Set the working directory
cd /scratch/negishi/iyer95/iyer95/

# Absolute Paths
input_dir="/scratch/negishi/iyer95/iyer95/Projects/AMP/analysis/ON/"
off_possible_dir="/scratch/negishi/iyer95/iyer95/Projects/AMP/analysis/OFF/off_possible/"
output_dir="/scratch/negishi/iyer95/iyer95/Projects/AMP/isomer_filter_6/"

# Create output directory
mkdir -p "$output_dir"

# Run the Python script
python core/python/isomer_filter_6_AMP.py \
    --input_dir "$input_dir" \
    --off_possible_dir "$off_possible_dir" \
    --output_dir "$output_dir" \
    --retention_time_tolerance 0.15 \
    > core/backend/logs/isomer/${SLURM_ARRAY_TASK_ID}_output.log 2>&1
