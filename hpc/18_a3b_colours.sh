#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_a3b_col
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-9
#SBATCH --output=logs/%x_%A_%a.out
#
# Assemblage-level colour (A3b) refitted once per Delhey plumage colour
# (one colour per array task). Needs fits/A3b_model_setup.RData and the
# reported A3b colourfulness fit. Summarise afterwards with
#   A3C_COL=summary $RSCRIPT scripts/cooney/A3b_colours.R
set -euo pipefail
source hpc/env.sh
COL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<< $'black\nwhite\ngrey\nbrown\nrufous\nyellow\nred\ngreen\nblue')
echo "==== A3b colour $COL  $(date +%T)"
export A3C_COL=$COL
$RSCRIPT scripts/cooney/A3b_colours.R
# Presence check for the zero-heavy colours (rufous, red, green, blue):
#   sbatch --array=5,7,8,9 --export=ALL,A3C_FAM=bernoulli hpc/18_a3b_colours.sh
