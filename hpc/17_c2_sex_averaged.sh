#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_c2
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-3
#SBATCH --output=logs/%x_%A_%a.out
#
# C2 sensitivity: refit A1b, A2d and A3b colourfulness models with the
# sex-averaged colourfulness (scripts/cooney/C2_sex_averaged.R). Then run
# the comparison as a short job:
#   sbatch --dependency=afterok:<this job> --account=DeiC-KU-L04 --mem=32G \
#     --time=02:00:00 --output=logs/%x_%j.out --job-name=colour_c2_summary \
#     --wrap="source hpc/env.sh && C2_MODEL=summary \$RSCRIPT scripts/cooney/C2_sex_averaged.R"
set -euo pipefail
source hpc/env.sh
export C2_MODEL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<< $'A1b\nA2d\nA3b')
echo "==== C2 $C2_MODEL  $(date +%T)"
$RSCRIPT scripts/cooney/C2_sex_averaged.R
