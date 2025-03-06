#!/bin/bash
# FILENAME:  run_LIGER

module load conda
conda activate /home/cbeveri/.conda/envs/cent7/2020.11-py38/CLAW

python prep_LIGER.py
python Run_LIGER2.py
python Run_LIGER1.py


conda deactivate
conda activate /home/cbeveri/.conda/envs/cent7/2020.11-py38/data_analysis_llm


python LIGER_AI_lit_search.py
