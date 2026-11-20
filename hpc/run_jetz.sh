#!/bin/bash
# Tree sensitivity: repeat every phylogenetic analysis on the Jetz et al.
# (2012) MCC tree (scripts/00c_jetz_mcc.R must have written
# data/external/jetz2012_hackett_mcc_passerines.tre). From the repo root:
#   bash hpc/run_jetz.sh [dependency job id]
# Runs inside runs/jetz/, which holds its own fits/, results/, figures/ and
# logs/ and links to the shared scripts/, data/ and hpc/; PHYLO=jetz makes
# 00_config.R read the Jetz tree. A3a/A3b have no phylogeny: their
# Claramunt-run fits are linked into runs/jetz/fits for the post-hoc tables.
# Finishes with scripts/cooney/T1_tree_comparison.R in the repo root.
set -euo pipefail
ROOT=$(pwd)
DEP=${1:+--dependency=afterok:$1}
export PHYLO=jetz
cd runs/jetz
SB="sbatch --parsable --account=DeiC-KU-L04 --output=logs/%x_%j.out"
R='source hpc/env.sh && $RSCRIPT'
j1=$(sbatch --parsable $DEP --job-name=jetz_species hpc/10_species.sh)
j4=$(sbatch --parsable --dependency=afterok:$j1 --job-name=jetz_mrpmm hpc/15_mrpmm.sh)
js=$($SB $DEP --job-name=jetz_a2d_setup --cpus-per-task=2 --mem=32G --time=02:00:00 \
       --wrap="$R scripts/cooney/A2d_01_setup.R")
ja=$(sbatch --parsable --dependency=afterok:$js --job-name=jetz_a2d hpc/11_a2d.sh)
jc=$(sbatch --parsable --dependency=afterok:$ja --job-name=jetz_a2d_collect hpc/12_a2d_collect.sh)
jk=$(sbatch --parsable --dependency=afterok:$js --array=2-10 --job-name=jetz_a2d_colours hpc/14_a2d_colours.sh)
jb=$(sbatch --parsable --dependency=afterok:$js --array=1 --job-name=jetz_a2d_black \
       --export=ALL,A2D_ITER=7500,A2D_WARMUP=2500 hpc/14_a2d_colours.sh)
jks=$($SB --dependency=afterok:$jk:$jb --job-name=jetz_a2d_colours_summary --cpus-per-task=2 \
       --mem=32G --time=04:00:00 --wrap="$R scripts/cooney/A2d_colours_03_summary.R")
jf=$($SB --dependency=afterok:$jc:$j4:$jks --job-name=jetz_posthoc --cpus-per-task=2 --mem=32G \
       --time=06:00:00 --wrap="$R scripts/cooney/posthoc_contrasts.R && $R scripts/cooney/posthoc_interactions.R \
         && $R scripts/cooney/pub_figures.R")
cd "$ROOT"
jt=$($SB --dependency=afterok:$jf --job-name=colour_tree_comparison --cpus-per-task=1 --mem=8G \
       --time=01:00:00 --wrap="$R scripts/cooney/T1_tree_comparison.R")
echo "jetz: species=$j1 mrpmm=$j4 a2d_setup=$js a2d=$ja collect=$jc colours=$jk black=$jb" \
     "colour_summary=$jks posthoc=$jf comparison=$jt"
