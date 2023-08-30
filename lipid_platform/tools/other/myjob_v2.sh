#!/bin/bash
# FILENAME:  myjob.sub


module load r

# --vanilla:
# --no-save: do not save datasets at the end of an R session
Rscript --vanilla --no-save edgeR_v2.r 
Rscript --vanilla --no-save edgeR_No_replicates.R 