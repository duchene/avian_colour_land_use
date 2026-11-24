# Avian plumage colour and land-use change

How land-use change relates to passerine plumage colouration, and whether colour
predicts abundance across land-use types, accounting for biome, trophic niche, body
mass, the PREDICTS study hierarchy, and phylogeny. Built on the PREDICTS biodiversity
database merged with Cooney et al. (2022) colour data (colourfulness measured from
multispectral photographs of museum specimens under an ultraviolet-sensitive visual
model, LociUVS, and plate-scored sexual dichromatism), the Delhey et al. (2023) plumage
colour classes, and AVONET morphological traits.

**Status: refitted.** The analytical dataset was rebuilt on 2026-08-11 after a filter
bug was corrected and Temperate Open was dropped, and all ten models have since been
fitted against it. `results/` and `figures/` hold the current outputs. The
pre-correction outputs are preserved in the git tag `prefilter-results-2026-08` and
should not be quoted. `ANALYSIS_NOTES.md` records the design decisions and the open
questions.

## 2026-11 revision (current)

After co-author review the reported analyses were revised. The decisions are recorded
in `ANALYSIS_NOTES.md` under "2026-11 revision"; in short:

- **Two colour metrics** in every reported model: male colourfulness (`malecolcooney`,
  male LociUVS) and sexual dichromatism (`dichroplate`, the Cooney et al. 2022 plate
  score), set by `RESP_MAIN` in `scripts/00_config.R`. The plate score replaces the
  male-minus-female LociUVS difference.
- **Specific colours** come from Delhey et al. (2023) as the proportion of the male
  plumage in each of ten colour classes (`data/external/`, built by
  `scripts/00b_external_colour.R`; helpers in `scripts/cooney/delhey_colours.R`). Each
  colour is a separate variable; no ratios or composite scores are formed.
- **Reported models**, refitted: A1b, A2d and A3b (with A3a for its LOO comparison).
  A2c, A1 and the Dale (B) family were not refitted and their results predate the
  revision.
- **New analyses**:

| Analysis | Script | Question |
|---|---|---|
| A2d per colour | `A2d_colours_02_fit.R`, `A2d_colours_03_summary.R` | A2d refitted once per colour class, colour proportion as the predictor |
| A4 MR-PMM | `A4_colour_mrpmm.R` (MCMCglmm) | Eight colours and four land-use scores as joint responses (Halliwell et al. 2025); correlations, phylogenetic signal, conditional colour on land-use slopes |
| A3b per colour | `A3b_colours.R` | Assemblage colour refitted once per colour (logit, z across species), with a presence (Bernoulli) check for rufous, red, green and blue |
| C1 sample bias | `C1_sample_bias.R` | Do the PREDICTS passerines differ in colour from all passerines? PGLS and KS tests |
| C2 sex-averaged | `C2_sex_averaged.R` | Do A1b, A2d and A3b reach the same conclusions with mean male and female colourfulness? |

- **Supplementary outputs** in `figures/pub/cooney/`: Figure S1 (tree with colour
  rings, `supp_fig_S1_tree.R`), S2 (sample bias, `C1_sample_bias.R`), S3 and S4
  (MR-PMM and per-colour effects, `supp_fig_S3_S4.R`), Tables S1 to S3
  (`supp_tables.R`) and S4 (assemblage colours, from `results/cooney/A3b_colours_*.csv`). Legends are in `figures/pub/cooney/Figure_and_Table_legends.txt`.
- **Two trees.** All phylogenetic analyses are repeated on the Jetz et al. (2012) MCC
  tree with identical data (`hpc/run_jetz.sh`, outputs in `runs/jetz/`), and compared
  with the reported Claramunt et al. (2025) runs by `T1_tree_comparison.R`
  (Table S3). 509 of 550 quantities reach the same conclusion under both trees. Two
  results reported in the main text are not credible with the Jetz tree: more black
  plumage in pasture-associated species (A4) and smaller losses of white-plumaged
  species in croplands (A2d per colour).
- The revision ran on GenomeDK through Slurm; see "Running on GenomeDK" below.

The rest of this README describes the pre-revision pipeline, which the scripts still
support.

## Analyses

Five models, each fitted twice: once with Cooney UVS colour (the **A** family,
primary) and once with independent Dale et al. plumage scores (the **B** family, a
colour-source robustness check). The two families are identical in data, priors,
sampler, and random structure, and differ only in the colour columns.

| Model | Question | Response | Phylogeny |
|---|---|---|---|
| **A1** / B1 | Species-level colour ~ land-use association | colour | yes |
| **A1b** / B1b | A1 with association scores from detections only (**reported**) | colour | yes |
| **A2c** / B2c | Relative abundance ~ colour x biome x land use | abundance | yes |
| **A2d** / B2d | Abundance change against each species' own primary-vegetation baseline | difference | yes |
| **A3b** / B3b | Community colour ~ biome x land use | colour | no |
| **A3a** / B3a | A3b plus land use x trophic niche and land use x body mass | colour | no |

A3b is the primary community model and biome x land use is reported from it. A3a adds
the guild and size interactions and is a nested superset, so leave-one-out
cross-validation adjudicates between them. A3a and A3b carry no phylogenetic term
because colour is constant within a species, which makes a per-species phylogenetic
intercept unidentifiable at the record level. A1 answers that question at the species
level instead.

**A1b / B1b are the reported species-level models, not A1 / B1.** A1 builds each
species' land-use association score by counting rows in the raw extract, and 84% of
those rows are structural zeros from the near-complete site by species matrix, so the
score largely reflects how a study spread its sites rather than where the bird was
found. A1b rebuilds the scores from detections only and changes nothing else. The two
predictors correlate at 0.58 to 0.79 depending on land use, and the difference decides
a result: pasture below primary vegetation is credible under A1 at -0.147
[-0.273, -0.023] and is not under A1b at -0.085 [-0.193, +0.021]. B1b shows the same
contrast is not credible under Dale scores either, at +0.024 [-0.012, +0.062], so the
two colour systems agree once both are scored from detections. A1 and B1 are kept
so the check stays auditable. They appear in `results/*/posthoc_contrasts.csv` but are
excluded from the published `figures/pub/*/Table_1_posthoc_contrasts.csv`, which
carries one model per analysis.

**A2e / A2f** are sensitivity checks that add the full three-way colour x biome x land
use interaction to A2c / A2d. They exist to justify reporting the two-way models.

Before the revision, colour was summarised by four metrics: mean colourfulness, male
colourfulness, the dichromatism ratio, and the male-minus-female difference. The
reported models now use male colourfulness and the plate-scored dichromatism (above).

## Structure

```
scripts/
  00_config.R              Every shared decision. Read this first.
  00_data_preparation.R    Raw data -> data/present.csv
  00_descriptive_figures.R Descriptive plots
  00b_external_colour.R    Cooney 2022 and Delhey 2023 tables -> data/external/
  cooney/                  A1 A1b A2c A2d A2e A2f A3a A3b, pub_figures.R,
                           2026-11: A2d_colours_*, A4_colour_mrpmm, C1, C2,
                           supp_fig_S1_tree, supp_fig_S3_S4, supp_tables
  dale/                    B1 B1b B2c B2d B3a B3b, exploratory_figures_dale.R
data/                      present.csv, raw .csv.gz, BBtree2.tre phylogeny
data/external/             Derived Cooney 2022 and Delhey 2023 species tables
hpc/                       Slurm job scripts for GenomeDK (run_all.sh submits everything)
results/cooney|dale/       Summary tables (CSV): convergence, fixed effects, variance, R2, LOO
figures/pub/cooney|dale/   Publication figures 3, 5, 6, Tables 1 to 2, and (cooney) the
                           supplementary Figures S1 to S4 and Tables S1 to S2
                           (1 and 2 are introductory and live outside this repo)
figures/exploratory/       Earlier exploratory figure set, not for publication
figures/cooney|dale/       Per-analysis diagnostics, with desc_* descriptives at the root
fits/                      brms model objects (gitignored, flat, prefixed)
```

`scripts/00_config.R` is the single source of truth. It owns the biome map and the
biome drop, the land-use reference level, the tree-matching synonyms, the priors, the
metric-to-family mapping, and the shared data loaders. Nothing else may reimplement
them. If you change a modelling decision, change it there and it propagates to all ten
analyses at once.

## Running

All scripts assume the working directory is the repository root. Each `_02_fit_models.R`
skips any model whose `fits/*.RData` already exists, so an interrupted run resumes from
the last saved fit. That also means **you must clear `fits/` before refitting on new
data**, or nothing will refit.

```r
source("scripts/00_data_preparation.R")   # check the printed tables

# Cooney (A), in publication order. A3b before A3a, because the A3a LOO loads the A3b fits.
source("scripts/cooney/A1_01_setup.R");  source("scripts/cooney/A1_02_fit_models.R");  source("scripts/cooney/A1_03_diagnostics_summary.R")
source("scripts/cooney/A1b_01_setup.R"); source("scripts/cooney/A1b_02_fit_models.R"); source("scripts/cooney/A1b_03_diagnostics_summary.R")
source("scripts/cooney/A2c_01_setup.R"); source("scripts/cooney/A2c_02_fit_models.R"); source("scripts/cooney/A2c_03_diagnostics_summary.R")
source("scripts/cooney/A2d_01_setup.R"); source("scripts/cooney/A2d_02_fit_models.R"); source("scripts/cooney/A2d_03_diagnostics_summary.R")
source("scripts/cooney/A3b_01_setup.R"); source("scripts/cooney/A3b_02_fit_models.R"); source("scripts/cooney/A3b_03_diagnostics_summary.R")
source("scripts/cooney/A3a_01_setup.R"); source("scripts/cooney/A3a_02_fit_models.R"); source("scripts/cooney/A3a_03_diagnostics_summary.R")
source("scripts/cooney/A2e_3way_sensitivity.R"); source("scripts/cooney/A2f_3way_sensitivity.R")
source("scripts/cooney/pub_figures.R")

# Dale (B), same order. B3b before B3a.
source("scripts/dale/B1_01_setup.R");  source("scripts/dale/B1_02_fit_models.R");  source("scripts/dale/B1_03_diagnostics_summary.R")
source("scripts/dale/B1b_01_setup.R"); source("scripts/dale/B1b_02_fit_models.R"); source("scripts/dale/B1b_03_diagnostics_summary.R")
source("scripts/dale/B2c_01_setup.R"); source("scripts/dale/B2c_02_fit_models.R"); source("scripts/dale/B2c_03_diagnostics_summary.R")
source("scripts/dale/B2d_01_setup.R"); source("scripts/dale/B2d_02_fit_models.R"); source("scripts/dale/B2d_03_diagnostics_summary.R")
source("scripts/dale/B3b_01_setup.R"); source("scripts/dale/B3b_02_fit_models.R"); source("scripts/dale/B3b_03_diagnostics_summary.R")
source("scripts/dale/B3a_01_setup.R"); source("scripts/dale/B3a_02_fit_models.R"); source("scripts/dale/B3a_03_diagnostics_summary.R")
source("scripts/dale/pub_figures_dale.R")
```

A3 and B3 models take minutes. Each phylogenetic fit (A2c, A2d, A2e, A2f and their B
counterparts) carries a dense species covariance matrix of about 1,150 by 1,150 and
takes a few hours. Run the scripts sequentially rather than in parallel shells, since
each model wants 4 chains times `nthreads` cores, and budget roughly 20 GB of disk for
a full `fits/`. Check `results/*/[AB]*_convergence.csv` after each family, and raise
`adapt_delta` before changing any model.

## Running on GenomeDK

The 2026-11 revision runs through Slurm, charged to the project account. From the
repository root:

```bash
bash hpc/run_all.sh     # submits every step with dependencies; see the header
```

`hpc/env.sh` defines `$RSCRIPT`, an R 4.4.1 container with the project libraries
(`tools/Rlib_colour` before `tools/Rlib`). Fit scripts skip fits that exist, so a rerun
only fills gaps. Do not edit an R script while a job is running it, since R reads the
script as it executes. Any script that loads a saved brmsfit must attach `rstan` first.

## Requirements

R packages: `brms`, `cmdstanr`, `tidyverse`, `posterior`, `loo`, `ape`, `phytools`,
`patchwork`, `scales`, plus `svglite` and `plotrix` for the figures. The revision
analyses also use `MCMCglmm` (A4), `phylolm` (C1) and `igraph` (built from source into
`tools/Rlib_colour` for phytools in the container).
