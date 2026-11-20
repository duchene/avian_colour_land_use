#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_smoke
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#
# Smoke test before the long refits: run the four setup scripts, then
# fit each new metric in each model type for a few iterations on a data
# subset. Catches data, family and compile errors in minutes, not hours.
# Usage, from the repo root:  sbatch hpc/01_smoke.sh
set -euo pipefail
source hpc/env.sh
for s in A1b A2d A3a A3b; do
  echo "==== setup $s"; $RSCRIPT scripts/cooney/${s}_01_setup.R
done
$RSCRIPT hpc/smoke_fits.R
