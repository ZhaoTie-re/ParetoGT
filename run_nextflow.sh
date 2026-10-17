#!/bin/bash
#SBATCH -p gr10478b
#SBATCH -t 7-00:00:00
#SBATCH --rsc p=1:t=1:c=4:m=8G
#SBATCH -o nextflow_head.log
#SBATCH -e nextflow_head.err

# ==============================================================================
# ParetoGT launcher: runs the Nextflow head as a SLURM job; the tasks are SLURM jobs of their own.
#
# Usage:
#   sbatch run_nextflow.sh [clean] [Nextflow options ...]
#   - clean: start afresh (without -resume)
#   - anything else is passed to Nextflow, e.g.
#       sbatch run_nextflow.sh --wgsDir <dir> --trueVcf <vcf>
#       sbatch run_nextflow.sh --matrixDir results/01.prepare_format_matrix --dp 6-12
# Progress: ./monitor_status.sh (reads results/pipeline_info/trace.txt and nextflow_head.log).
# Author: Zhao Tie
# ==============================================================================

# Nextflow 25.04 or later (DSL2). The environment's bin goes first on PATH: at this site an activated
# conda environment is not reliably first.
export PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH

# Java heap of the head process: 8G job, 7G heap (many tasks to schedule)
export NXF_OPTS='-Xms2g -Xmx7g'

RESUME_FLAG="-resume"
NF_ARGS=()
for arg in "$@"; do
    if [[ "$arg" == "clean" ]]; then
        echo "Starting afresh (no -resume)..."
        RESUME_FLAG=""
    else
        NF_ARGS+=("$arg")
    fi
done
[ -n "$RESUME_FLAG" ] && echo "Resuming from the last run..."

# trace, report and timeline: results/pipeline_info/ (nextflow.config)
nextflow run main.nf $RESUME_FLAG "${NF_ARGS[@]}"
