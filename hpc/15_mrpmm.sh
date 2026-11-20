#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_mrpmm
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#
# A4 colour MR-PMM (MCMCglmm, 4 chains in parallel). Needs
# fits/A1b_model_setup.RData.
set -euo pipefail
source hpc/env.sh
export MRPMM_NITT=${MRPMM_NITT:-205000} MRPMM_BURN=${MRPMM_BURN:-5000} MRPMM_THIN=${MRPMM_THIN:-100}
$RSCRIPT scripts/cooney/A4_colour_mrpmm.R
