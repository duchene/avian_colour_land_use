# Rerun note: corrected primary-vegetation filter

Written 2026-08-05. Read before running anything on the new machine.

## What was wrong

`scripts/00_data_preparation.R` built `fulldat_primary` (studies containing primary
vegetation), then wrote `present.csv` from `fulldat` instead. The filter never reached the
saved data, so every analysis in the repo (A1a, A1b, A2, A2b, A2c, A2d, A2e, A2f, A3 and
the whole B family) was fitted on the unfiltered file.

Confirmed three ways: `present.csv` had 34,653 rows, exactly the non-zero-abundance row
count of the raw file, with a matching SSBS multiset; 16 of 54 references (8,561 records,
24.7%) held no primary-vegetation record; and those 16 held no primary-vegetation sites in
the raw file either, including its 182,181 zero rows, so this was not an artefact of
dropping absences.

## The fix

`00_data_preparation.R` now runs three reporting steps: drop zero-abundance rows
(216,834 to 34,653), keep only units with primary vegetation **and** at least one other
land use (34,653 to 26,092), then save and print the land-use and biome-by-land-use
tables. A `stopifnot` asserts the criterion holds, so it cannot regress silently.

Two choices, both commented in the script. The filter runs *after* dropping absences,
because testing the raw file only asks whether primary-vegetation sites were sampled.
`FILTER_UNIT <- "Reference"` matches the original code; `"SS"` gives 26,060 instead of
26,092, so the choice is immaterial. I also dropped the unused `site_study` line, which
was computed on the discarded object and never reached `present.csv`.

Dry-run verified against a scratch path. **It has not been run against
`data/present.csv`**, which still holds the old 34,653-row file, so the existing fits stay
interpretable until you deliberately replace it.

## The dataset you will get

| | Old (wrong) | New |
|---|---|---|
| Records | 34,653 | 26,092 |
| Species | 1,684 | 1,437 |
| References / studies | 54 / 59 | 38 / 41 |
| Blocks / sites | 361 / 3,716 | 205 / 2,377 |
| A1 fitted N (Cooney) | 33,210 | 24,838 |
| B1 fitted N (Dale) | 34,481 | 25,920 |

| Land use | Old | New | Change |
|---|---|---|---|
| Primary vegetation | 10,920 | 10,920 | 0 |
| Secondary | 8,645 | 5,769 | -33% |
| Plantation forest | 7,080 | 6,524 | -8% |
| Cropland | 4,577 | 2,225 | **-51%** |
| Pasture | 3,431 | 654 | **-81%** |

Phylogenetic subsets shrink too. Stage-1 (`jetz_sp`) matching on the corrected data gives
1,145 species / 20,491 records, and with Cooney colour roughly 1,032 species / 19,423
records for A2c, against 1,348 / 26,951 before. The setup scripts recount this themselves.

## Decide before fitting

**Pasture and Cropland lose most of their data.** Cerezo et al. 2011 alone contributed
1,340 Cropland and 1,387 Pasture records and is now excluded, leaving Pasture at 654
records. The pasture results were load-bearing (ANALYSIS_NOTES.md finding 10, the sorting
argument), so treat every pasture and cropland conclusion as unverified until the refits
land.

**The sparse-cell collapse needs revisiting.** The setup scripts hardcode one collapse,
Temperate Open x Plantation forest into Secondary, under a threshold of roughly <50 records
or <15 species. On the corrected data the whole Temperate Open non-primary row is thin:

| Biome4 | Primary | Secondary | Plantation | Cropland | Pasture |
|---|---|---|---|---|---|
| Tropical Forest | 4,733 | 3,304 | 3,690 | 1,585 | 326 |
| Tropical Open | 2,330 | 1,706 | 2,367 | 315 | 96 |
| Temperate Forest | 3,742 | 723 | 425 | 325 | 197 |
| Temperate Open | 115 | **36** | **42** | **0** | **35** |

Species counts for that row: 32, 23, 12, 0, 10. Temperate Open x Cropland is now **empty**,
where it held 1,339 records before, so the empty interaction cell that `A1b_01_setup.R`
comments on moves from Plantation to Cropland, and the existing collapse leaves a Temperate
Open Secondary cell of only 78 records. Options, in my order of preference:

1. Collapse all Temperate Open non-primary land uses into one non-primary level: keeps the
   biome, one honest contrast of 113 records against 115 primary.
2. Drop Temperate Open (228 records, 0.9%), making the interaction fully estimable.
   Simplest to defend.
3. Merge Temperate Open into Temperate Forest for a three-level Biome4.

I have not changed this. It is a modelling decision, and it must be identical across A1a,
A1b, A2c, A2d, A2e, A2f and all five B scripts, so set it in one place.

## Critical: clear the old fits first

`fits/` holds 91 files and 18 GB, and every fit script starts with
`if (file.exists(f)) { ... return(fit) }`, so **with the old fits present nothing refits and
you silently rebuild the wrong results.** `fits/` is gitignored, so a fresh clone is clean,
but copying the Dropbox folder carries it across.

```sh
mv fits fits_PREFILTER_2026-08-05 && mkdir fits
cp -r results results_PREFILTER_2026-08-05
cp -r figures figures_PREFILTER_2026-08-05
```

`results/` and `figures/` are tracked and get overwritten in place, so snapshot and commit
them first: the diff then shows exactly what the filter changed.

## Run order

From the repository root. Sampler settings unchanged.

```r
source("scripts/00_data_preparation.R")   # check the printed tables against this note

# Cooney (A). A1b before A1a (LOO loads A1b fits); A3 before A2c (synonym table).
source("scripts/cooney/A1b_01_setup.R"); source("scripts/cooney/A1b_02_fit_models.R"); source("scripts/cooney/A1b_03_diagnostics_summary.R")
source("scripts/cooney/A1a_01_setup.R"); source("scripts/cooney/A1a_02_fit_models.R"); source("scripts/cooney/A1a_03_diagnostics_summary.R")
source("scripts/cooney/A3_01_setup.R");  source("scripts/cooney/A3_02_fit_models.R");  source("scripts/cooney/A3_03_diagnostics_summary.R")
source("scripts/cooney/A2c_01_setup.R"); source("scripts/cooney/A2c_02_fit_models.R"); source("scripts/cooney/A2c_03_diagnostics_summary.R")
source("scripts/cooney/A2d_01_setup.R"); source("scripts/cooney/A2d_02_fit_models.R"); source("scripts/cooney/A2d_03_diagnostics_summary.R")
source("scripts/cooney/A2e_3way_sensitivity.R"); source("scripts/cooney/A2f_3way_sensitivity.R")
source("scripts/cooney/pub_figures.R")

# Dale (B), same order
source("scripts/dale/B1b_01_setup.R"); source("scripts/dale/B1b_02_fit_models.R"); source("scripts/dale/B1b_03_diagnostics_summary.R")
source("scripts/dale/B1a_01_setup.R"); source("scripts/dale/B1a_02_fit_models.R"); source("scripts/dale/B1a_03_diagnostics_summary.R")
source("scripts/dale/B3_01_setup.R");  source("scripts/dale/B3_02_fit_models.R");  source("scripts/dale/B3_03_diagnostics_summary.R")
source("scripts/dale/B2c_01_setup.R"); source("scripts/dale/B2c_02_fit_models.R"); source("scripts/dale/B2c_03_diagnostics_summary.R")
source("scripts/dale/B2d_01_setup.R"); source("scripts/dale/B2d_02_fit_models.R"); source("scripts/dale/B2d_03_diagnostics_summary.R")
source("scripts/dale/pub_figures_dale.R")
```

The README lists A2c before A3, but `A2c_01_setup.R` reads
`results/cooney/A3_synonym_table.csv`, so run A3 first.

Machine notes. The constraint is the phylogenetic models (A2c, A2d, A2e, A2f and their B
counterparts), each carrying a dense species covariance matrix, previously 1,348 x 1,348
and now around 1,032 x 1,032, so per-chain memory should fall by about a third. Those took
a few hours each before. Each model wants 4 chains x `nthreads` cores, so run the scripts
sequentially rather than in parallel shells, and budget 25 GB of disk for a fresh `fits/`.
The A1 and B1 models are minutes. Run detached with logs captured, then check
`results/*/[AB]*_convergence.csv` after each family: with Pasture at 654 records and
several thin biome cells, divergences and low ESS are likelier than before. Raise
`adapt_delta` before touching the model.

## Update after the refits

- `ANALYSIS_NOTES.md:7` and `:11`, plus every N in the results sections.
- `draft_methods_results.txt:18` onward. It states n = 34,653 for the site-level analyses.
  Two errors there: the filter was not applied, and even on the old data the Cooney models
  fit 33,210 rows. State the fitted N per family.
- `README.md`, if the biome or collapse decision changes. It also says the analyses account
  for body mass, but mass enters only A1a, A1b, B1a and B1b, not A2c, A2d or A3.
- Passerine coverage: 1,437 species now, roughly 24% of BirdTree passerines rather than 28%.
