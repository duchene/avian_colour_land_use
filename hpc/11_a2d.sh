#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_a2d
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-2
#SBATCH --output=logs/%x_%A_%a.out
#
# A2d (abundance change vs primary, phylogeny): one colour metric per
# array task. The setup must exist first (hpc/run_all.sh runs it), and the
# collect step (12_a2d_collect.sh) gathers the fits.
set -euo pipefail
source hpc/env.sh
CV=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<< $'z_malecolcooney\nz_dichroplate')
echo "==== A2d $CV  $(date +%T)"
export ONLY=$CV   # apptainer passes the host environment through
$RSCRIPT scripts/cooney/A2d_02_fit_models.R
