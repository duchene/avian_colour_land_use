# Avian plumage colour and land-use change

How land-use change relates to passerine plumage colouration, and whether colour
predicts abundance across land-use types, accounting for biome, trophic niche, body
mass, the PREDICTS study hierarchy, and phylogeny. Built on the PREDICTS biodiversity
database merged with Cooney spectrophotometric colour data (ultraviolet-sensitive
visual model, LociUVS) and AVONET morphological traits.

**Status: awaiting refit.** The analytical dataset was rebuilt on 2026-08-11 after a
filter bug was corrected and Temperate Open was dropped. No model has been fitted
against it yet, so `results/` and `figures/` are empty by design. The pre-correction
outputs are preserved in the git tag `prefilter-results-2026-08` and should not be
quoted. `ANALYSIS_NOTES.md` records the design decisions and the open questions.

## Analyses

Five models, each fitted twice: once with Cooney UVS colour (the **A** family,
primary) and once with independent Dale et al. plumage scores (the **B** family, a
colour-source robustness check). The two families are identical in data, priors,
sampler, and random structure, and differ only in the colour columns.

| Model | Question | Response | Phylogeny |
|---|---|---|---|
| **A1b** / B1b | Community colour ~ biome x land use | colour | no |
| **A1a** / B1a | A1b plus land use x trophic niche and land use x body mass | colour | no |
| **A2c** / B2c | Relative abundance ~ colour x biome x land use | abundance | yes |
| **A2d** / B2d | Abundance change against each species' own primary-vegetation baseline | difference | yes |
| **A3** / B3 | Species-level colour ~ land-use association | colour | yes |

A1b is the primary community model and biome x land use is reported from it. A1a adds
the guild and size interactions and is a nested superset, so leave-one-out
cross-validation adjudicates between them. A1a and A1b carry no phylogenetic term
because colour is constant within a species, which makes a per-species phylogenetic
intercept unidentifiable at the record level. A3 answers that question at the species
level instead.

**A2e / A2f** are sensitivity checks that add the full three-way colour x biome x land
use interaction to A2c / A2d. They exist to justify reporting the two-way models.

Colour is summarised by four metrics: mean colourfulness, male conspicuousness, the
sexual dichromatism ratio, and the sexual dichromatism difference.

## Structure

```
scripts/
  00_config.R              Every shared decision. Read this first.
  00_data_preparation.R    Raw data -> data/present.csv
  00_descriptive_figures.R Descriptive plots
  cooney/                  A1a A1b A2c A2d A2e A2f A3, pub_figures.R
  dale/                    B1a B1b B2c B2d B3, pub_figures_dale.R
data/                      present.csv, raw .csv.gz, BBtree2.tre phylogeny
results/cooney|dale/       Summary tables (CSV): convergence, fixed effects, variance, R2, LOO
figures/pub/cooney|dale/   Publication figures 1 to 7
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

# Cooney (A). A1b before A1a, because the A1a LOO loads the A1b fits.
source("scripts/cooney/A1b_01_setup.R"); source("scripts/cooney/A1b_02_fit_models.R"); source("scripts/cooney/A1b_03_diagnostics_summary.R")
source("scripts/cooney/A1a_01_setup.R"); source("scripts/cooney/A1a_02_fit_models.R"); source("scripts/cooney/A1a_03_diagnostics_summary.R")
source("scripts/cooney/A2c_01_setup.R"); source("scripts/cooney/A2c_02_fit_models.R"); source("scripts/cooney/A2c_03_diagnostics_summary.R")
source("scripts/cooney/A2d_01_setup.R"); source("scripts/cooney/A2d_02_fit_models.R"); source("scripts/cooney/A2d_03_diagnostics_summary.R")
source("scripts/cooney/A3_01_setup.R");  source("scripts/cooney/A3_02_fit_models.R");  source("scripts/cooney/A3_03_diagnostics_summary.R")
source("scripts/cooney/A2e_3way_sensitivity.R"); source("scripts/cooney/A2f_3way_sensitivity.R")
source("scripts/cooney/pub_figures.R")

# Dale (B), same order. B1b before B1a.
source("scripts/dale/B1b_01_setup.R"); source("scripts/dale/B1b_02_fit_models.R"); source("scripts/dale/B1b_03_diagnostics_summary.R")
source("scripts/dale/B1a_01_setup.R"); source("scripts/dale/B1a_02_fit_models.R"); source("scripts/dale/B1a_03_diagnostics_summary.R")
source("scripts/dale/B2c_01_setup.R"); source("scripts/dale/B2c_02_fit_models.R"); source("scripts/dale/B2c_03_diagnostics_summary.R")
source("scripts/dale/B2d_01_setup.R"); source("scripts/dale/B2d_02_fit_models.R"); source("scripts/dale/B2d_03_diagnostics_summary.R")
source("scripts/dale/B3_01_setup.R");  source("scripts/dale/B3_02_fit_models.R");  source("scripts/dale/B3_03_diagnostics_summary.R")
source("scripts/dale/pub_figures_dale.R")
```

A1 and B1 models take minutes. Each phylogenetic fit (A2c, A2d, A2e, A2f and their B
counterparts) carries a dense species covariance matrix of about 1,150 by 1,150 and
takes a few hours. Run the scripts sequentially rather than in parallel shells, since
each model wants 4 chains times `nthreads` cores, and budget roughly 20 GB of disk for
a full `fits/`. Check `results/*/[AB]*_convergence.csv` after each family, and raise
`adapt_delta` before changing any model.

## Requirements

R packages: `brms`, `cmdstanr`, `tidyverse`, `posterior`, `loo`, `ape`, `phytools`,
`patchwork`, `scales`, plus `svglite` and `plotrix` for the figures.
