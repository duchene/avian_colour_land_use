#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_a2d_collect
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
# Gather the two A2d fits and write the A2d diagnostics and summaries.
set -euo pipefail
source hpc/env.sh
$RSCRIPT scripts/cooney/A2d_02_fit_models.R
$RSCRIPT scripts/cooney/A2d_03_diagnostics_summary.R
