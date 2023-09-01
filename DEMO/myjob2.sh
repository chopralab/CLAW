#!/bin/bash

# Activate the R environment
source /home/sanjay/anaconda3/etc/profile.d/conda.sh
conda activate R

# Run the R script
R --vanilla --no-save < edgeR.r
