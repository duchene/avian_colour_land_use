#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_a2d_col
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-10
#SBATCH --output=logs/%x_%A_%a.out
#
# A2d refitted once per Delhey plumage colour (one colour per array task).
# Needs fits/A2d_model_setup.RData. Summarise afterwards with
# scripts/cooney/A2d_colours_03_summary.R (hpc/16_a2d_colours_summary.sh).
set -euo pipefail
source hpc/env.sh
COL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<< $'black\nwhite\ngrey\nbrown\nrufous\nyellow\nred\ngreen\nblue\npurple')
echo "==== A2d colour $COL  $(date +%T)"
export ONLY=$COL   # apptainer passes the host environment through
$RSCRIPT scripts/cooney/A2d_colours_02_fit.R
