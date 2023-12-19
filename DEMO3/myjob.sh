#!/bin/bash
# FILENAME:  myjob.sub

module load r

# --vanilla:
# --no-save: do not save datasets at the end of an R session
R --vanilla --no-save < edgeR.r