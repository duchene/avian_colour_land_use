#!/bin/bash
# Submit the 2026-11 revision analyses with dependencies. From the repo root:
#   bash hpc/run_all.sh
# Every fit script skips fits that already exist in fits/, so a rerun only
# fills gaps. Order of dependencies:
#   species (A1b + C1 bias test) -> A4 MR-PMM
#   A2d setup -> A2d (2 metrics) -> collect
#             -> A2d per colour (10 fits) -> per-colour summary
#   A3 (A3b, A3a)
#   A1b + A2d + A3b -> C2 sex-averaged colourfulness check -> summary
#   everything -> publication and supplementary figures, post-hoc tables
set -euo pipefail
source hpc/env.sh
SB="sbatch --parsable --account=DeiC-KU-L04 --output=logs/%x_%j.out"
R='source hpc/env.sh && $RSCRIPT'

j1=$(sbatch --parsable hpc/10_species.sh)
j4=$(sbatch --parsable --dependency=afterok:$j1 hpc/15_mrpmm.sh)
js=$($SB --job-name=colour_a2d_setup --cpus-per-task=2 --mem=32G --time=02:00:00 \
       --wrap="$R scripts/cooney/A2d_01_setup.R")
ja=$(sbatch --parsable --dependency=afterok:$js hpc/11_a2d.sh)
jc=$(sbatch --parsable --dependency=afterok:$ja hpc/12_a2d_collect.sh)
# colours 2-10 at the A2d chain length; black (task 1) needs longer chains
# to reach Rhat < 1.01 (see ANALYSIS_NOTES.md)
jk=$(sbatch --parsable --dependency=afterok:$js --array=2-10 hpc/14_a2d_colours.sh)
jb=$(sbatch --parsable --dependency=afterok:$js --array=1 \
       --export=ALL,A2D_ITER=7500,A2D_WARMUP=2500 hpc/14_a2d_colours.sh)
jks=$($SB --dependency=afterok:$jk:$jb --job-name=colour_a2d_colours_summary --cpus-per-task=2 \
       --mem=32G --time=04:00:00 --wrap="$R scripts/cooney/A2d_colours_03_summary.R")
j3=$(sbatch --parsable hpc/13_a3.sh)
j5=$(sbatch --parsable --dependency=afterok:$j1:$jc:$j3 hpc/17_c2_sex_averaged.sh)
j6=$($SB --dependency=afterok:$j5 --job-name=colour_c2_summary --cpus-per-task=2 --mem=32G \
       --time=02:00:00 --wrap="source hpc/env.sh && C2_MODEL=summary \$RSCRIPT scripts/cooney/C2_sex_averaged.R")
jf=$($SB --dependency=afterok:$j4:$jks:$j6 --job-name=colour_figures --cpus-per-task=2 --mem=32G \
       --time=06:00:00 --wrap="$R scripts/cooney/posthoc_contrasts.R && $R scripts/cooney/posthoc_interactions.R \
         && $R scripts/cooney/pub_figures.R && $R scripts/cooney/supp_fig_S1_tree.R \
         && $R scripts/cooney/supp_fig_S3_S4.R && $R scripts/cooney/supp_tables.R")
echo "species=$j1 mrpmm=$j4 a2d_setup=$js a2d=$ja a2d_collect=$jc a2d_colours=$jk black=$jb" \
     "colour_summary=$jks a3=$j3 c2=$j5 c2_summary=$j6 figures=$jf"
