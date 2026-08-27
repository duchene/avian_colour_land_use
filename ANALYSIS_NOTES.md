# Analysis notes

Design decisions and model specifications. This file records what the analyses do and
why. It deliberately holds no fitted numbers. Results live in `results/`, and quoting
them from prose here is what allowed the documentation to drift out of step with the
code in the first place.

Last restructured 2026-08-11.

## Status

**The refit is complete** (2026-08-11 16:08 to 2026-08-14 02:51, 58.7 hours, zero
failed steps). All 33 models across both colour families converged with zero divergent
transitions and a maximum Rhat of 1.012. Results are in `results/`, figures in
`figures/pub/`, and the write-up is `draft_methods_results.txt`.

Fitted numbers are deliberately absent from this file, because quoting results in prose
is what let it drift out of step with the code.

**One diagnostic problem is outstanding**, described under Open decisions: every
species-level model in both families returned a low E-BFMI.

Observed runtimes on this machine (16 cores, 4 chains at 4 threads), per model: A3b 17
minutes, A1 28 minutes, A3a 50 minutes, A2c 80 minutes, A2f about 400 minutes, and A2d
373 minutes. A2d, B2d and A2f account for most of the wall time, since each carries a
dense species covariance matrix over 53,551 records. Budget about 60 hours for a full
A and B sequence.

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

## Why A3a and A3b carry no phylogenetic term

A3b was originally specified with `(1 | gr(phylo, cov = A))`. A feasibility fit showed
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

A1 places a phylogenetic term on colour legitimately because it is fitted at the
species level, one row per species, where the covariance matrix constrains the
phylogenetic effect and sigma captures the residual. The abundance models A2c and A2d
keep the term because abundance genuinely varies within species. Community level
without phylogeny and species level with phylogeny is the only identifiable way to ask
both questions.

## Model specifications

### A3b, B3b: community colour by biome and land use

```
colour ~ Biome4 * Predominant_simple [+ Trophic.Niche + z_logMass]
         + (1|SS) + (1|SSB) + (1|SSBS)
```

Two variants per response, base and covariate. The primary community model. Biome by
land use is reported from here. Trophic niche and body mass enter as main effects only,
since their land-use interactions are A3a's job.

### A3a, B3a: adding guild and size interactions

```
colour ~ Biome4 * Predominant_simple + Trophic.Niche + z_logMass
         + Predominant_simple:z_logMass + Predominant_simple:Trophic.Niche
         + (1|SS) + (1|SSB) + (1|SSBS)
```

A nested superset of the A3b covariate model on identical rows, so leave-one-out
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

The full `colour * Biome4 * Predominant_simple` interaction fitted for male colourfulness and
compared to the reported two-way model by leave-one-out cross-validation. A2e extends
A2c, A2f extends A2d. These exist to justify reporting the two-way models and are not
themselves reported. The check uses male colourfulness only. Extend it to the dichromatism
difference if an all-metric statement is needed.

### A1, B1: species-level phylogenetic regression

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

**A1 and B1 compute land-use proportions from the full raw file.** They therefore
include studies with no primary vegetation and records in the dropped Temperate Open
biome, while every other analysis works from the filtered `present.csv`. The rationale
is that a species' habitat association is a property of the species and is best
estimated from all available records. The alternative is consistency with the rest of
the pipeline. Flagged in the header of both scripts.

**A1 posterior contrasts: resolved 2026-08-13.** The compositional coefficients are
strongly correlated in the posterior, so overlapping marginal credible intervals do not
settle whether pasture-associated species differ from primary-associated ones. The
contrast is now computed inside `A1_03_diagnostics_summary.R` and `B1_03_...` and
written to `results/*/[AB]3_landuse_contrasts.csv`, so it can no longer be skipped.
The answer: the pasture difference is credible for colourfulness where the marginal
intervals had suggested it was not, which vindicates the earlier claim but only because
it was finally tested the right way.

**Mediterranean Forests, Woodlands and Scrub sits in Temperate Forest.** Its WWF name
spans forest and shrubland, so the assignment is genuinely borderline. Changing it means
editing `BIOME_MAP` in the config and refitting.

**The species-level model has an E-BFMI problem. This is the one thing still to fix.**
All four A1 models and all four B1 models return an energy Bayesian fraction of missing
information below 0.3, on all four chains in seven of the eight cases. Because it
appears in both colour families it is a property of the species-level phylogenetic
regression itself, not of either metric. Rhat (max 1.005) and bulk ESS (min 819) look
fine, which is the point: neither detects this.

It matters because A1's one substantive result, the pasture contrast, is a statement
about a posterior tail. Refit A1 and B1 with a non-centred parameterisation of the
phylogenetic term, or at higher `adapt_delta` with a longer warmup, and confirm the
contrast. Note that the pasture contrast does not reproduce under Dale scores, so there
is already independent reason to doubt it.

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
until the Dale replication is in, since B3a fits the same interaction on a different
colour metric and offers an independent read.

**Effects at the credibility boundary.** Several A3a guild interactions have intervals
that only just exclude zero. Report the interval rather than the verdict, and do not
build an argument on a term whose upper bound is within rounding distance of zero.

**Cross-family comparison is the real test.** The A-versus-B comparison on the previous
data found three of four analyses replicating, with the trophic interactions and part of
the A2c conclusion not carrying over. Redo that comparison from scratch on the corrected
fits rather than assuming the earlier verdict holds.

## Removed analyses

Recoverable with `git checkout prefilter-results-2026-08 -- <path>`.

- **Original A3.** A single `(1|SSBS)` random intercept. Refitted as A3a under the full
  hierarchy.
- **A2.** Relative abundance on colour by land use with a single `(1|SSBS)` intercept
  and no biome or phylogeny. Superseded by A2c.
- **A2b.** Paired differences with `(1|SSBS) + (1|Reference)` and no biome or
  phylogeny. Superseded by A2d.

Publication Figure 4 compared A2 with A2c to show colour effects vanishing under the
fuller model. It was removed with A2, and the remaining figures renumbered 1 to 7 in
both families. Making that point in a sentence of text costs nothing and avoids
maintaining a whole superseded analysis to feed one panel.


## A1 land-use association scores count absences (found 2026-08-27)

A1 builds each species' land-use association score in `scripts/cooney/A1_01_setup.R`
with `count()` over the raw extract. That file is 84% structural zeros, because
PREDICTS stores a near-complete site by species matrix, and `count()` counts rows
rather than detections. The score therefore measures how a study distributed its
sites across land uses, not where the species was actually recorded.

The mechanism is visible in a single study. In `VK1_2012__Otto 1` every one of the
20 species has the identical row vector, 116 primary and 27 secondary, while their
detections range from 13/2 to 53/9. 908 species, 53% of the total, occur in exactly
one study, so for them the predictor is a constant that says nothing about the
species at all.

Rebuilding the scores from detections only gives a genuinely different predictor.
The two correlate at 0.58 for primary vegetation, 0.61 for cropland, 0.64 for
plantation, 0.71 for secondary and 0.79 for pasture, and 218 species move their
pasture score by more than 0.1. The species-level cost is small: 1,684 species with
at least one detection against 1,703 in the extract, 1,375 matched to tips against
1,391, and 1,242 entering the model against 1,254.

A1b is that alternative, identical to A1 in every other respect.

### Resolved, 2026-08-27: A1b is reported

A1b converged cleanly, with a maximum Rhat of 1.0032, a minimum bulk ESS of 788 and
no divergent transitions. It carries the same E-BFMI warning as A1, below 0.3 on all
four chains for mean colourfulness, male colourfulness and sexual dichromatism and on
two chains for the dichromatism ratio, so that diagnostic is a property of the
species-level phylogenetic regression rather than of the predictor.

The structural findings do not move. Bayesian R2 is 86.4% [82.1, 90.1] for male
colourfulness against A1's 86.3%, and the phylogenetic standard deviation is still
3.2 times the residual (0.563 against 0.174).

The land-use result does. Pasture below primary vegetation for male colourfulness was
-0.147 [-0.273, -0.023] under A1, the only credible contrast in a reported metric. In
A1b it attenuates to -0.085 [-0.193, +0.021] and the posterior probability of a
negative effect falls from 0.99 to 0.94. Only 1 of 16 A1b contrasts is credible, mean
colourfulness in pasture at -0.113 [-0.209, -0.017], a supplementary metric. The sign
is unchanged throughout, so this is attenuation and not reversal.

A second problem also disappears. Under A1 the coefficients were extrapolations,
because the highest cropland score was 0.50 and the highest pasture score 0.80, so no
species sat at the corner the coefficient describes. Under A1b all five land uses
contain species at a score of 1: 128 primary vegetation, 124 secondary, 25 plantation,
24 cropland and 16 pasture.

Figure 3 and the post-hoc tables now come from A1b. A1 is retained in
`results/cooney/posthoc_contrasts.csv` under model "A1" so the check stays auditable.


### B1b, the Dale replication of the same check (2026-08-27)

B1b converged cleanly, maximum Rhat 1.0046, minimum bulk ESS 799, no divergent
transitions, and carries the same E-BFMI warning on all four chains of all four
models as the Cooney family does.

Dale male colourfulness barely moves between the two constructions: pasture against
primary vegetation is +0.028 [-0.014, +0.072] under B1 and +0.024 [-0.012, +0.062]
under B1b, neither credible. So the earlier "does not replicate" reading was wrong in
its diagnosis. The Cooney estimate was the one that moved, and once both families are
scored from detections they agree that land-use association does not predict
colourfulness.

B1b does turn up something B1 did not. Species associated with pasture are more
sexually dimorphic than species associated with primary vegetation, at +1.86 Dale
units [+0.06, +3.63], probability 0.98, with the dichromatism ratio agreeing at
+0.038 [+0.003, +0.073]. The Cooney counterpart runs the same way at +6.82 LociUVS
[-3.14, +16.84] but only reaches probability 0.90. This is the one land-use effect on
colour that gains support under Dale scores rather than losing it. Credible on one
family alone, with a lower bound at 0.06 and an E-BFMI warning outstanding, it is not
established.
