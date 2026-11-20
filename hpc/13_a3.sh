#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_a3
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
# Assemblage models: A3b (reported) and its LOO partner A3a, two metrics each (male colourfulness, sexual dichromatism).
set -euo pipefail
source hpc/env.sh
for s in A3b_01_setup A3b_02_fit_models A3b_03_diagnostics_summary \
         A3a_01_setup A3a_02_fit_models A3a_03_diagnostics_summary; do
  echo "==== $s  $(date +%T)"; $RSCRIPT scripts/cooney/$s.R
done
