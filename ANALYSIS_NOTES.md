# Analysis notes

Design decisions and model specifications. This file records what the analyses do and
why. It deliberately holds no fitted numbers. Results live in `results/`, and quoting
them from prose here is what allowed the documentation to drift out of step with the
code in the first place.

Last restructured 2026-08-11.

## Status

The refit began 2026-08-11 16:08. The whole Cooney family (A1b, A1a, A2c, A2d, A3),
both sensitivity checks, and the Cooney figures are complete, all with zero divergent
transitions. In the Dale family B1b and B1a are complete, B2c is fitting, and B2d and B3
follow. Progress is logged to `logs/` (gitignored).

Fitted numbers live in `results/` and are written up in `draft_methods_results.txt`.
They are deliberately absent here, because quoting results in prose is what let this
file drift out of step with the code.

Observed runtimes on this machine (16 cores, 4 chains at 4 threads), per model: A1b 17
minutes, A3 28 minutes, A1a 50 minutes, A2c 80 minutes, and A2d 373 minutes. A2d and its
B counterpart are the constraint, at 53,551 records against A2c's 19,785, 3,500
iterations, and `adapt_delta` 0.95. Together with A2f they account for most of the
wall time. Budget about 60 hours for the full A and B sequence.

## What changed on 2026-08-11

**The primary-vegetation filter now reaches the saved data.** `00_data_preparation.R`
built a filtered object and then wrote `present.csv` from the unfiltered one, so every
analysis in the repository had been fitted on 34,653 records that should have been
filtered. A `stopifnot` now asserts the criterion on the saved file, so this cannot
recur silently.

**Temperate Open was dropped.** On the corrected data it holds 115 primary-vegetation
records against 36, 42, 0, and 35 in the four modified land uses, and Temperate Open by
Cropland is empty. No collapse leaves a defensible contrast. Dropping the biome costs
about 0.9% of records and makes biome by land use fully estimable, so the sparse-cell
collapse that each setup script used to hardcode is gone. **Biome4 now carries three
levels.**

**Shared decisions moved into `scripts/00_config.R`.** The biome map had been copied
into three scripts, the tree-synonym table into six, the prior blocks into twelve, and
the sparse-cell collapse into eight, in two mutually inconsistent variants. Each now
exists once.

**A2 and A2b were removed.** Both were superseded by A2c and A2d, which add biome, the
full PREDICTS hierarchy, and a phylogenetic random effect. Their scripts and results
are recoverable from the archive tag.

## Data

### Raw

PREDICTS site-level bird surveys merged with Cooney et al. spectrophotometric plumage
measurements (ultraviolet-sensitive visual model, LociUVS) and AVONET traits. Stored as
`data/aves_predicts_passerines_dale_cooney_avonet.csv.gz`, holding 216,834 records from
1,703 passerine species.

### Analytical dataset

`data/present.csv`, built by `00_data_preparation.R` in four steps.

1. Drop zero-abundance records. PREDICTS stores a near-complete site by species matrix
   per study, so about 84% of raw rows are structural zeros. This leaves 34,653 rows.
2. Apply the biome policy: assign Biome4 and drop Temperate Open. This leaves 31,259
   rows.
3. Keep only publications holding primary vegetation plus at least one other land use,
   since a unit without that contrast cannot inform a land-use effect. This leaves
   25,864 rows from 37 of 49 publications.
4. Set the land-use reference level and save.

The biome drop runs before the qualifying filter, so a publication is judged on the
records that actually enter the models.

**Resulting dataset:** 25,864 records, 1,420 species, 40 studies, 203 blocks, 2,335
sites. Land use: Primary vegetation 10,805, Plantation forest 6,482, Secondary 5,733,
Cropland 2,225, Pasture 619.

All fifteen biome by land-use cells are estimable, the thinnest being Tropical Open by
Pasture at 96 records and 34 species.

**Pasture is thin.** It fell from 3,431 records to 619, largely because Cerezo et al.
2011 contributed 1,340 Cropland and 1,387 Pasture records and holds no primary
vegetation. Treat pasture effects with corresponding caution.

### Colour metrics

Four Cooney metrics for the A family, computed in `00_data_preparation.R`:
`meancolcooney` (mean of male and female LociUVS), `malecolcooney` (male LociUVS),
`dichrocooney` (male over female, a ratio), and `dichrodiff` (male minus female, which
can be negative). The four Dale counterparts for the B family are `meancoldale`,
`malecoldale`, `dichrodale`, and `dichrodiffdale`. Dale coverage is broader, at 1,679
species with a male score against 1,530 for Cooney.

### Phylogeny

BigBirdTree (BBtree2) from [evolucionario/BigBirdTree](https://github.com/evolucionario/BigBirdTree),
9,072 tips under Jetz-era taxonomy, stored as `data/BBtree2.tre`.

Species are matched in three stages by `match_to_tree()` in the config: the `jetz_sp`
column, then the underscored `Best_guess_binomial`, then a curated table of fifteen
genus reclassifications. Tips claimed by more than one species are dropped. Species
absent from the tree are dropped rather than grafted onto congeners, because grafting
invents branch lengths that distort the covariance matrix.

## Shared design decisions

All of these live in `scripts/00_config.R`.

**Random structure.** The canonical PREDICTS hierarchy `(1|SS) + (1|SSB) + (1|SSBS)`,
meaning study, block within study, and site within block. `SS` is the PREDICTS study
code, not `Reference`, because one publication can hold several studies. Since `SSBS`
is globally unique and nested, no `/` nesting syntax is needed.

**Biome.** Eleven WWF biomes collapse to three: Tropical Forest (reference), Tropical
Open, and Temperate Forest. `assign_biome4()` errors on any unlisted biome instead of
bucketing it silently, which is how the old nested `ifelse` could quietly reassign
records.

**Land use.** Five levels with Primary vegetation as reference: Secondary, Plantation
forest, Cropland, Pasture. In A2d and B2d the primary level is absorbed into the
response, leaving four levels with Secondary as reference.

**Interactions.** Two-way only. The three-way colour by biome by land-use term is
omitted everywhere and tested separately by A2e and A2f.

**Families.** Lognormal for the three positive colour metrics and for relative
abundance. Gaussian for the male-minus-female difference and for the paired difference.

**Phylogeny and colour compete.** Where colour enters as a predictor (A2c, A2d), the
phylogenetic random effect and colour claim the same species-level variance, because
colour is strongly conserved. Colour fixed effects there should be read as effects
beyond phylogeny. Biome and land use vary within a species and are not affected.

## Why A1a and A1b carry no phylogenetic term

A1b was originally specified with `(1 | gr(phylo, cov = A))`. A feasibility fit showed
the specification is not identifiable, and the term was dropped.

Plumage colour is a species-level trait, identical across all of a species' occurrence
records. A phylogenetic term with one level per species therefore behaves as a
per-species random intercept fitted to a response with zero within-species variance. It
absorbs essentially all between-species variation, while the site and study random
effects, the residual, and every biome by land-use fixed effect collapse toward zero.
The diagnostics showed structural degeneracy rather than slow mixing: sd(phylo) was
0.35 while the three site-level standard deviations and sigma were all near zero, with
Rhat between 1.8 and 3.0 and bulk effective sample sizes of 2 to 7 despite no
divergences.

A3 places a phylogenetic term on colour legitimately because it is fitted at the
species level, one row per species, where the covariance matrix constrains the
phylogenetic effect and sigma captures the residual. The abundance models A2c and A2d
keep the term because abundance genuinely varies within species. Community level
without phylogeny and species level with phylogeny is the only identifiable way to ask
both questions.

## Model specifications

### A1b, B1b: community colour by biome and land use

```
colour ~ Biome4 * Predominant_simple [+ Trophic.Niche + z_logMass]
         + (1|SS) + (1|SSB) + (1|SSBS)
```

Two variants per response, base and covariate. The primary community model. Biome by
land use is reported from here. Trophic niche and body mass enter as main effects only,
since their land-use interactions are A1a's job.

### A1a, B1a: adding guild and size interactions

```
colour ~ Biome4 * Predominant_simple + Trophic.Niche + z_logMass
         + Predominant_simple:z_logMass + Predominant_simple:Trophic.Niche
         + (1|SS) + (1|SSB) + (1|SSBS)
```

A nested superset of the A1b covariate model on identical rows, so leave-one-out
cross-validation tests directly whether the guild and size interactions earn their
place. Both scripts call `load_community_data()`, which guarantees the shared rows the
comparison requires.

Trophic niche has an Aquatic predator reference and includes Herbivore terrestrial at
two records, whose land-use interactions are prior-only and flagged in the fixed-effects
table.

### A2c, B2c: relative abundance

```
relabund ~ colour * Biome4 + colour * Predominant_simple + Biome4 * Predominant_simple
           + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))
```

Relative abundance is a species' effort-corrected count over the site total. Colour
enters standardised. 20,850 records across 1,157 tips.

### A2d, B2d: abundance change against a within-species baseline

```
diff_abund ~ colour * Biome4 + colour * Predominant_simple + Biome4 * Predominant_simple
             + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))
```

Each non-primary record is its species' relative abundance minus that species' mean
relative abundance in primary vegetation within the same study. Anchoring each species
to its own baseline admits absences, which the presence-only A2c cannot do.

Built from the raw file rather than `present.csv`, because the design needs species
that disappeared, giving a negative difference, or appeared, giving a positive one.
Differences of exactly zero mean the species was absent from both and carry no
information, so they are dropped. 53,551 records across 1,154 tips from 39 studies,
from 66,225 non-zero differences of which 54,038 are negative.

### A2e, A2f: three-way sensitivity checks

The full `colour * Biome4 * Predominant_simple` interaction fitted for male colour and
compared to the reported two-way model by leave-one-out cross-validation. A2e extends
A2c, A2f extends A2d. These exist to justify reporting the two-way models and are not
themselves reported. The check uses male colour only. Extend it to the dichromatism
difference if an all-metric statement is needed.

### A3, B3: species-level phylogenetic regression

```
colour ~ 0 + prop_Cropland + prop_Pasture + prop_Plantation_forest
             + prop_Primary_vegetation + prop_Secondary
         + (1 | gr(phylo, cov = A))
```

One row per species. Predictors are within-species land-use proportions, the records of
a species in a land use over that species' total records. Because they sum to 1 per
species they span the intercept, which is therefore removed, and each coefficient reads
as the expected colour of a species found exclusively in that land use. 1,391 species
matched to 1,391 tips.

## Open decisions

These need a call before or during the refit. None has been made silently.

**A3 and B3 compute land-use proportions from the full raw file.** They therefore
include studies with no primary vegetation and records in the dropped Temperate Open
biome, while every other analysis works from the filtered `present.csv`. The rationale
is that a species' habitat association is a property of the species and is best
estimated from all available records. The alternative is consistency with the rest of
the pipeline. Flagged in the header of both scripts.

**A3 posterior contrasts: resolved 2026-08-13.** The compositional coefficients are
strongly correlated in the posterior, so overlapping marginal credible intervals do not
settle whether pasture-associated species differ from primary-associated ones. The
contrast is now computed inside `A3_03_diagnostics_summary.R` and `B3_03_...` and
written to `results/*/[AB]3_landuse_contrasts.csv`, so it can no longer be skipped.
The answer: the pasture difference is credible for colourfulness where the marginal
intervals had suggested it was not, which vindicates the earlier claim but only because
it was finally tested the right way.

**Mediterranean Forests, Woodlands and Scrub sits in Temperate Forest.** Its WWF name
spans forest and shrubland, so the assignment is genuinely borderline. Changing it means
editing `BIOME_MAP` in the config and refitting.

**A3 has an E-BFMI problem, and it got worse.** On the corrected data all four models
return an energy Bayesian fraction of missing information below 0.3, on all four chains
for mean colour, male colour and the dichromatism difference, and on one chain for the
ratio. Previously only the dichromatism difference was affected. Rhat (max 1.002) and
bulk ESS (min 863) look fine, which is the point: neither detects this.

It matters because A3's one substantive result, the pasture contrast, is a statement
about a posterior tail. Before that result is submitted, refit A3 with a non-centred
parameterisation of the phylogenetic term, or at higher `adapt_delta` with a longer
warmup, and confirm the contrast is unchanged. B3 will show whether the problem is
specific to the Cooney metrics or inherent to the species-level model.

## Watchlist: results sensitive to the correction

Which findings should be read cautiously, and why. This is a record of where the data
are thin, not a summary of results.

**Anything resting on pasture or cropland.** Pasture fell from 3,431 records to 619 and
cropland from 4,577 to 2,225, because the publications that supplied most of them held
no primary vegetation and are correctly excluded. Any effect whose credibility depends
on those two land uses now rests on far less data, even where the interval still
excludes zero.

**The mass by land-use interactions changed sign on the corrected data**, and they sit
squarely in the cropland and pasture cells above. Treat the direction as unsettled
until the Dale replication is in, since B1a fits the same interaction on a different
colour metric and offers an independent read.

**Effects at the credibility boundary.** Several A1a guild interactions have intervals
that only just exclude zero. Report the interval rather than the verdict, and do not
build an argument on a term whose upper bound is within rounding distance of zero.

**Cross-family comparison is the real test.** The A-versus-B comparison on the previous
data found three of four analyses replicating, with the trophic interactions and part of
the A2c conclusion not carrying over. Redo that comparison from scratch on the corrected
fits rather than assuming the earlier verdict holds.

## Removed analyses

Recoverable with `git checkout prefilter-results-2026-08 -- <path>`.

- **Original A1.** A single `(1|SSBS)` random intercept. Refitted as A1a under the full
  hierarchy.
- **A2.** Relative abundance on colour by land use with a single `(1|SSBS)` intercept
  and no biome or phylogeny. Superseded by A2c.
- **A2b.** Paired differences with `(1|SSBS) + (1|Reference)` and no biome or
  phylogeny. Superseded by A2d.

Publication Figure 4 compared A2 with A2c to show colour effects vanishing under the
fuller model. It was removed with A2, and the remaining figures renumbered 1 to 7 in
both families. Making that point in a sentence of text costs nothing and avoids
maintaining a whole superseded analysis to feed one panel.
