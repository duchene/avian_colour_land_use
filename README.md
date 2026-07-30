# Avian plumage colour and land-use change

How land-use change relates to passerine plumage colouration, and whether colour
predicts abundance across land-use types — accounting for biome, trophic niche, body
mass, the spatial study/block/site structure, and phylogeny. Built on the PREDICTS
biodiversity database merged with Cooney spectrophotometric colour data (ultraviolet-
sensitive visual model, LociUVS) and AVONET morphological traits.

## Analyses

The manuscript focuses on five models across four analyses. `ANALYSIS_NOTES.md` has
the full methods and results, including the simpler earlier variants (A2, A2b) that
these supersede.

- **A1b** — Community plumage colour ~ biome × land use (+ trophic niche + body mass),
  with the full PREDICTS random structure `(1|SS)+(1|SSB)+(1|SSBS)`. The primary
  community-colour model. No phylogenetic term: colour is a species-level constant, so
  a per-species phylogenetic intercept is not identifiable at the record level (that
  question is answered by A3).
- **A1a** — Extends A1b by adding land use × trophic niche and land use × body mass
  interactions, resolving guild- and size-specific colour responses (e.g. frugivores
  in plantation forest). Identical data, priors, sampler and random structure to A1b,
  so it is a nested superset and LOO tests whether those interactions earn their place.
  Biome × land use is reported from A1b.
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

Every focal analysis has a Dale-colour counterpart, prefixed **B**: B1a, B1b, B2c, B2d,
and B3 mirror A1a, A1b, A2c, A2d, and A3, with Dale et al. plumage scores in place of
the Cooney UVS metrics. Same data, priors, sampler, random structure, and metric sets,
differing only in the colour source and the `B`-prefixed outputs. This tests whether the
findings hold under an independent colour-scoring system.

## Structure

```
data/            Analytical dataset (present.csv), raw data (.csv.gz), phylogeny (BBtree2.tre)
scripts/         Shared data prep (00_*.R) at the root, then two parallel families:
  cooney/        Cooney-colour analyses A1a / A1b / A2 / A2b / A2c / A2d / A2e / A2f / A3, pub_figures.R
  dale/          Dale-colour analyses B1a / B1b / B2c / B2d / B3, pub_figures_dale.R
results/
  cooney/        Cooney (A) summary tables (CSV): convergence, fixed effects, variance, R2, LOO
  dale/          Dale (B) summary tables, same set
figures/
  pub/cooney/    Cooney publication figures (Figures 1-8)
  pub/dale/       Dale publication figures (parallel numbering; Figure 4 has no Dale counterpart)
  cooney/, dale/ Per-analysis diagnostic plots; desc_* descriptive plots at the root
fits/            brms model objects (gitignored, kept flat and prefixed — regenerate via the scripts)
```

The Cooney (A) and Dale (B) families are deliberately kept separate on disk: same
pipeline, one folder each, so the two colour systems never share an output path.

## Requirements

R packages: `brms`, `cmdstanr`, `tidyverse`, `posterior`, `loo`, `ape`, `phytools`,
`patchwork`, `scales` (plus `svglite` and `plotrix` for the figures).

## Running

All scripts assume the working directory is the repository root. Prepare the data,
then run each focal analysis (setup → fit → diagnostics), the sensitivity checks, and
finally the figures:

```r
source("scripts/00_data_preparation.R")   # requires the raw .csv.gz in data/; writes both colour systems

# --- Cooney (A) family: scripts/cooney/ ---
# A1b — community colour (primary; run before A1a so its fits exist for the LOO)
source("scripts/cooney/A1b_01_setup.R"); source("scripts/cooney/A1b_02_fit_models.R"); source("scripts/cooney/A1b_03_diagnostics_summary.R")
# A1a — community colour extension (land use x trophic and x mass); LOO vs A1b
source("scripts/cooney/A1a_01_setup.R"); source("scripts/cooney/A1a_02_fit_models.R"); source("scripts/cooney/A1a_03_diagnostics_summary.R")
# A2c — relative abundance (reuses A3 tree-matching)
source("scripts/cooney/A2c_01_setup.R"); source("scripts/cooney/A2c_02_fit_models.R"); source("scripts/cooney/A2c_03_diagnostics_summary.R")
# A2d — paired-difference abundance change
source("scripts/cooney/A2d_01_setup.R"); source("scripts/cooney/A2d_02_fit_models.R"); source("scripts/cooney/A2d_03_diagnostics_summary.R")
# A3 — phylogenetic regression
source("scripts/cooney/A3_01_setup.R"); source("scripts/cooney/A3_02_fit_models.R"); source("scripts/cooney/A3_03_diagnostics_summary.R")
# Sensitivity checks (three-way interaction)
source("scripts/cooney/A2e_3way_sensitivity.R"); source("scripts/cooney/A2f_3way_sensitivity.R")
# Cooney publication figures
source("scripts/cooney/pub_figures.R")

# --- Dale (B) family: scripts/dale/, the same analyses with Dale colour ---
# Run B1b before B1a (B1a's LOO needs the B1b fits), as with A1b/A1a.
source("scripts/dale/B1b_01_setup.R"); source("scripts/dale/B1b_02_fit_models.R"); source("scripts/dale/B1b_03_diagnostics_summary.R")
source("scripts/dale/B1a_01_setup.R"); source("scripts/dale/B1a_02_fit_models.R"); source("scripts/dale/B1a_03_diagnostics_summary.R")
source("scripts/dale/B2c_01_setup.R"); source("scripts/dale/B2c_02_fit_models.R"); source("scripts/dale/B2c_03_diagnostics_summary.R")
source("scripts/dale/B2d_01_setup.R"); source("scripts/dale/B2d_02_fit_models.R"); source("scripts/dale/B2d_03_diagnostics_summary.R")
source("scripts/dale/B3_01_setup.R"); source("scripts/dale/B3_02_fit_models.R"); source("scripts/dale/B3_03_diagnostics_summary.R")
# Dale publication figures
source("scripts/dale/pub_figures_dale.R")
```

The phylogenetic models (A2c/A2d and their sensitivity checks) each take a few hours;
model objects are large and are not tracked in git.
