# Avian plumage colour and land-use change

Analyses of how land-use change affects passerine plumage colouration and whether colour predicts abundance across land-use types. Uses the PREDICTS biodiversity database merged with Cooney spectrophotometric colour data and Avonet morphological traits.

## Structure

```
data/           Analytical dataset (present.csv)
scripts/        R scripts for data preparation and analyses
results/        Summary tables (CSV)
figures/        Plots (PNG)
fits/           brms model objects (not tracked — regenerate via scripts)
```

## Analyses

- **Analysis 1 (A1):** Colour ~ land use + body mass + habitat + biome + trophic niche
- **Analysis 2 (A2):** Abundance ~ colour × land use

See `ANALYSIS_NOTES.md` for detailed methods and current results.

## Requirements

R packages: `brms`, `cmdstanr`, `tidyverse`, `bayesplot`, `posterior`, `loo`

## Running

All scripts assume the working directory is the repository root. Run in order:

```r
source("scripts/00_data_preparation.R")  # requires raw CSV in parent directory
source("scripts/A1_01_setup.R")
source("scripts/A1_02_fit_models.R")
source("scripts/A1_03_model_comparison.R")
source("scripts/A1_04_diagnostics.R")
source("scripts/A1_05_summarize.R")
source("scripts/A2_01_setup.R")
source("scripts/A2_02_fit_models.R")
source("scripts/A2_03_diagnostics_summary.R")
```
