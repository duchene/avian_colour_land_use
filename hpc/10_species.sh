#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_species
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#
# Species-level analyses: A1b (colour ~ land-use association, phylogeny)
# for the primary metrics, its diagnostics and the sample-bias test (C1).
# Usage, from the repo root:  sbatch hpc/10_species.sh
set -euo pipefail
source hpc/env.sh
st=$(date +%s)
for s in cooney/A1b_01_setup cooney/A1b_02_fit_models cooney/A1b_03_diagnostics_summary \
         cooney/C1_sample_bias; do
  echo "==== $s  $(date +%T)"; $RSCRIPT scripts/$s.R
done
echo "done in $(( ($(date +%s) - st) / 60 )) min"
