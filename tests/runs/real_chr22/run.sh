#!/bin/bash
# Real-data check of the DSL2 pipeline: chr22, 16 threshold combinations, on SLURM, from the v5 matrices
# (the v5 WGS VCFs no longer exist); the results are compared with those of the v5 run.
set -euo pipefail
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/ParetoGT
export PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH NXF_ANSI_LOG=false
nextflow run main.nf -work-dir tests/runs/real_chr22/work \
  --matrixDir /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/tuning.concordance.v5/results/01.prepare_format_matrix \
  --chromosomes 22 --dp 1,8 --gq 10,20 --laf 0.0,0.2 --haf 1.0,0.8 --combosPerTask 8 \
  --outDir tests/runs/real_chr22/results
