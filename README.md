# Avian plumage colour and land-use change

How land-use change relates to passerine plumage colouration, and whether colour
predicts abundance across land-use types — accounting for biome, trophic niche, body
mass, the spatial study/block/site structure, and phylogeny. Built on the PREDICTS
biodiversity database merged with Cooney spectrophotometric colour data (ultraviolet-
sensitive visual model, LociUVS) and AVONET morphological traits.

## Analyses

The manuscript focuses on four models. `ANALYSIS_NOTES.md` has the full methods and
results, including the simpler earlier variants (A1, A2, A2b) that these supersede.

- **A1b** — Community plumage colour ~ biome × land use (+ trophic niche + body mass),
  with the full PREDICTS random structure `(1|SS)+(1|SSB)+(1|SSBS)`. No phylogenetic
  term: colour is a species-level constant, so a per-species phylogenetic intercept is
  not identifiable at the record level (that question is answered by A3).
- **A2c** — Relative abundance ~ colour × biome × land use (all two-way interactions),
  with the PREDICTS random structure and a phylogenetic random effect.
- **A2d** — Paired-difference abundance change (vs each species' primary-vegetation
  baseline) ~ colour × biome × land use, same structure as A2c.
- **A3** — Phylogenetic regression of species-level colour on land-use association.

**A2e / A2f** are sensitivity checks that add the full three-way `colour × biome ×
land use` interaction to A2c / A2d; it is not favoured by LOO, so the two-way models
are reported.

Colour is summarised by four metrics: mean colourfulness, male conspicuousness, the
sexual dichromatism ratio, and the sexual dichromatism difference.

## Structure

```
data/      Analytical dataset (present.csv), raw data (.csv.gz), phylogeny (BBtree2.tre)
scripts/   R scripts: data prep, analyses (A1b / A2c / A2d / A3, + A2e/A2f), figures
results/   Summary tables (CSV): convergence, fixed effects, variance components, R2, LOO
figures/   Plots; figures/pub holds the publication figures
fits/       brms model objects (gitignored — regenerate via the scripts)
```

## Requirements

R packages: `brms`, `cmdstanr`, `tidyverse`, `posterior`, `loo`, `ape`, `phytools`,
`patchwork`, `scales` (plus `svglite` and `plotrix` for the figures).

## Running

All scripts assume the working directory is the repository root. Prepare the data,
then run each focal analysis (setup → fit → diagnostics), the sensitivity checks, and
finally the figures:

```r
source("scripts/00_data_preparation.R")   # requires the raw .csv.gz in data/

# A1b — community colour
source("scripts/A1b_01_setup.R"); source("scripts/A1b_02_fit_models.R"); source("scripts/A1b_03_diagnostics_summary.R")
# A2c — relative abundance (reuses A3 tree-matching)
source("scripts/A2c_01_setup.R"); source("scripts/A2c_02_fit_models.R"); source("scripts/A2c_03_diagnostics_summary.R")
# A2d — paired-difference abundance change
source("scripts/A2d_01_setup.R"); source("scripts/A2d_02_fit_models.R"); source("scripts/A2d_03_diagnostics_summary.R")
# A3 — phylogenetic regression
source("scripts/A3_01_setup.R"); source("scripts/A3_02_fit_models.R"); source("scripts/A3_03_diagnostics_summary.R")

# Sensitivity checks (three-way interaction)
source("scripts/A2e_3way_sensitivity.R"); source("scripts/A2f_3way_sensitivity.R")

# Publication figures
source("scripts/pub_figures.R")
```

The phylogenetic models (A2c/A2d and their sensitivity checks) each take a few hours;
model objects are large and are not tracked in git.
