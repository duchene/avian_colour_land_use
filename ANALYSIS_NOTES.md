# Analysis notes

## Data

### Raw data

Source: PREDICTS database (site-level bird surveys across land-use types worldwide) merged with Cooney spectrophotometric plumage data (UVS visual model) and Avonet traits (body mass, trophic niche, habitat). Stored as `data/aves_predicts_passerines_dale_cooney_avonet.csv.gz` (gzipped, 32MB; 216,834 records, 1,703 species). R reads directly via `read.csv(gzfile(...))`.

### Analytical dataset (A1/A2)

Filtered to passerines with non-zero abundance from studies that include primary vegetation plus at least one other land-use type. N = 34,653 records. Stored as `data/present.csv`.

### Colour variables (computed in `00_data_preparation.R`)

- `meancolcooney`: mean of male + female LociUVS (overall colourfulness)
- `dichrocooney`: male / female LociUVS ratio (sexual dichromatism)
- `malecolcooney`: male LociUVS (male conspicuousness)

All use LociUVS (ultraviolet-sensitive visual model loci counts) from Cooney et al. Dale scores were excluded.

### Land-use levels (`Predominant_simple`)

Cropland (39,482 records), Pasture (32,407), Plantation forest (38,868), Primary vegetation (55,452), Secondary (50,625).

### Phylogenetic tree

BigBirdTree (BBtree2) from [evolucionario/BigBirdTree](https://github.com/evolucionario/BigBirdTree). 9,072 tips (Newick format). Stored as `data/BBtree2.tre`. Uses Jetz-era taxonomy.

## Analysis 1: Colour ~ Land-use

### Question

Does plumage colouration differ across land-use types after accounting for body mass, habitat, biome, and trophic niche?

### Methods

Bayesian GLMMs via brms/CmdStan. For each of the three colour responses:

```
colour ~ Predominant_simple + z_logMass + Habitat + Biome + Trophic.Niche +
         Predominant_simple:z_logMass +
         Predominant_simple:Habitat +
         Predominant_simple:Trophic.Niche +
         (1 | SSBS)
```

- `Predominant_simple`: land-use category (reference = cropland). Levels: Primary vegetation, Secondary, Pasture, Plantation forest.
- `z_logMass`: standardized log body mass.
- `SSBS`: site nested within study nested within block (random intercept).
- All three responses fitted with **lognormal** family. Dichrocooney was initially fitted with Student-t family but had severe non-convergence (Rhat 1.5-1.7, ESS < 10); the ratio data were well-behaved (range 0.43-4.44) and lognormal resolved convergence completely.
- Priors: Normal(0, 2) on intercept, Normal(0, 0.5) on fixed effects, Exponential(2) on SD and sigma.
- 4 chains, 3,000 iterations (1,500 warmup), adapt_delta = 0.90, threading enabled.

A **minimal** model (drops Habitat interaction) was also fitted for LOO-CV comparison.

Scripts: `A1_01_setup.R`, `A1_02_fit_models.R`, `A1_03_model_comparison.R`, `A1_04_diagnostics.R`, `A1_05_summarize.R`.

### Model comparison

Reduced (with Habitat interaction) preferred over Minimal for all three responses by LOO-CV.

### Convergence

All models converged (Rhat < 1.01, ESS > 400, 0 divergences), except dichrocooney minimal which was borderline (Rhat 1.011).

### Key results

**Variance explained:** R² = 17.8% (meancolcooney), 16.0% (malecolcooney), 12.5% (dichrocooney).

**Land-use main effects:** No land-use category has a 95% CI excluding zero for any colour metric. Land use alone does not shift community colour.

**Body mass:** Consistent negative effect across all three metrics (estimate -0.08 to -0.09 for mean/male colour, -0.017 for dichromatism). Larger birds are less colourful.

**Trophic niche:** Strongest predictor.
- Nectarivores and frugivores are significantly more colourful (+0.23 to +0.35).
- Frugivores, granivores, and invertivores are significantly less dichromatic (-0.19 to -0.21).

**Interactions:**
- Plantation forest x Frugivore: frugivores in plantations are significantly less colourful (mean -0.38, male -0.35). Strongest interaction in the dataset.
- Primary vegetation x Nectarivore: unexpectedly lower colour (-0.28 mean, -0.32 male).
- Plantation forest x Grassland/Human Modified habitats: consistently negative for colour.
- Woodland habitats show reduced dichromatism in plantation forest and primary vegetation.

## Analysis 2: Abundance ~ Colour x Land-use

### Question

Do more colourful species have different abundance across land-use types?

### Methods

Bayesian GLMMs via brms/CmdStan. For each standardized colour predictor (z_meancolcooney, z_dichrocooney, z_malecolcooney):

```
abundance ~ z_colour * Predominant_simple + (1 | SSBS)
```

- Family: lognormal.
- Priors: Normal(0, 2) on intercept, Normal(0, 0.5) on fixed effects, Exponential(2) on SD and sigma.
- 4 chains, 3,000 iterations (1,500 warmup), adapt_delta = 0.90, threading enabled.

Scripts: `A2_01_setup.R`, `A2_02_fit_models.R`, `A2_03_diagnostics_summary.R`.

### Convergence

All three models converged (Rhat < 1.01, ESS > 400, 0 divergences).

### Key results

**Variance explained:** R² ~11% for all models (land use drives most abundance variation).

**Main colour effect:** Positive for all three metrics (~0.03 on log scale, 95% CI just excludes zero). More colourful species tend to be slightly more abundant overall.

**Colour x land-use interactions:**
- **Pasture:** positive interaction -- colourful species are disproportionately abundant in pastures (credible for male colour and dichromatism).
- **Plantation forest:** negative interaction -- colourful species are relatively less abundant (credible for all three metrics, estimates -0.05 to -0.06).
- **Primary vegetation:** negative interaction -- similar to plantation forest (credible for all three, estimates -0.03 to -0.05).
- **Secondary:** no credible interaction.

## Analysis 3: Phylogenetic regression -- Colour ~ land-use association

### Question

At the species level, and accounting for shared ancestry, do species that are more associated with particular land-use types differ in colouration? This complements A1 (site-level) by asking whether evolutionary lineages with stronger affinity for degraded or pristine habitats have distinct colour profiles.

### Species-level data (computed in `A3_01_setup.R`)

From the full raw dataset (216,834 records, 1,703 species), for each species and each level of `Predominant_simple`, a **land-use proportion score** is computed as:

```
proportion = (records of species i in land-use j) / (total records in dataset)
```

This gives 5 predictor variables per species: `prop_Cropland`, `prop_Pasture`, `prop_Plantation_forest`, `prop_Primary_vegetation`, `prop_Secondary`. Some species may have 0 for certain land uses.

Colour variables are species-level traits (same value per record of a given species), so the first non-NA value per species is used.

### Taxonomic matching

Tree tip labels use Jetz-era taxonomy. Species were matched in three stages:

1. **jetz_sp column** (primary): 1,367 species matched directly.
2. **Best_guess_binomial** with underscores (fallback): 15 additional species.
3. **Synonym lookup** (genus reclassifications): 15 additional species matched via a curated table of cases where the same biological species appears under a different genus name in the tree due to taxonomic reclassification (e.g., `Trichastoma_bicolor` -> `Pellorneum_bicolor`). The full synonym table is in `results/A3_synonym_table.csv`.

Three species pairs where multiple BirdLife names map to the same Jetz species (Anthus cinnamomeus/richardi, Camaroptera brachyura/brevicaudata, Zosterops kirki/maderaspatanus) were dropped to avoid ambiguity.

**Final dataset:** 1,391 species matched to 1,391 tree tips (81.7% of 1,703). The remaining 306 unmatched species have valid Jetz names that are simply absent from the 9,072-tip BBtree2 tree. No species were grafted onto congeners, as this would introduce arbitrary branch lengths that distort the phylogenetic covariance matrix.

The crosswalk files `birdlife-birdtree_crosswalk.csv` and `clements_jetz_crosswalk.csv` (identical content; 11,291 rows mapping BirdLife to BirdTree names) were used to investigate unmatched species. All 306 remaining unmatched species have their jetz_sp name present in the crosswalk as valid BirdTree names but absent from BBtree2.

### Methods

Bayesian phylogenetic regression via brms/CmdStan. For each of the three colour responses:

```
colour ~ prop_Cropland + prop_Pasture + prop_Plantation_forest +
         prop_Primary_vegetation + prop_Secondary +
         (1 | gr(phylo, cov = A))
```

- `A`: phylogenetic covariance matrix from `ape::vcv.phylo()` (correlation form).
- `(1 | gr(phylo, cov = A))`: phylogenetic random effect accounting for shared ancestry.
- Tree trimmed to species present in both data and tree using `ape::drop.tip()`.
- Family: lognormal (all three responses).
- Priors: Normal(0, 2) on intercept, Normal(0, 1) on fixed effects (wider than A1/A2 because predictors are proportions on a different scale), Exponential(2) on SD and sigma.
- 4 chains, 4,000 iterations (2,000 warmup), adapt_delta = 0.95, max_treedepth = 12, threading enabled (4 threads per chain).
- N = 1,254 species per model (after excluding species with missing colour data).

Scripts: `A3_01_setup.R`, `A3_02_fit_models.R`, `A3_03_diagnostics_summary.R`.

### Convergence

All three models converged: 0 divergences, max Rhat <= 1.003, min bulk ESS >= 1,305, min tail ESS >= 2,524. The malecolcooney model produced an E-BFMI < 0.3 warning (common with large phylogenetic covariance matrices) but Rhat and ESS values were satisfactory.

### Key results

**Variance explained:** R² = 89.2% (meancolcooney, 95% CI: 85.3-92.4%), 86.7% (malecolcooney, 82.5-90.3%), 46.7% (dichrocooney, 37.2-56.7%). The high R² is driven almost entirely by the phylogenetic random effect (shared ancestry), not the land-use predictors.

**Land-use proportion effects:** No land-use proportion score is a credible predictor of any colour metric. All estimates are near zero (range -0.07 to +0.07) with very wide 95% credible intervals spanning approximately +/-2. The phylogenetic signal dominates -- closely related species have similar colouration regardless of which land-use types they occupy.

**Phylogenetic signal:** The phylogenetic SD is large relative to the residual: sd(phylo) = 0.52 (meancolcooney), 0.56 (malecolcooney), 0.26 (dichrocooney), compared to sigma = 0.14, 0.17, 0.17 respectively. This confirms strong phylogenetic conservatism in plumage colouration, especially for overall colourfulness. Dichromatism shows weaker (but still substantial) phylogenetic signal.

## Summary interpretation (A1 + A2 + A3)

Land use does not directly shift community colour in a simple way. Effects are mediated by trophic ecology and habitat context, and there is no evolutionary signal linking land-use affinity to colour:

1. **Trophic niche is a much stronger predictor of colour than land use** (A1). Nectarivores and frugivores are the most colourful guilds.
2. **Plantation forests filter against colourful frugivores** (A1). The plantation forest x frugivore interaction is the strongest in the dataset.
3. **Colourful species are slightly more abundant overall, but this advantage is land-use dependent** (A2). The positive colour-abundance relationship disappears in forested land uses (plantation + primary) and is amplified in pastures.
4. **Larger-bodied birds are consistently less colourful** (A1).
5. **Plumage colour is strongly phylogenetically conserved** (A3). After accounting for shared ancestry, species-level land-use association has no detectable effect on colouration. The community-level patterns in A1/A2 reflect ecological filtering (which species occur where) rather than evolutionary divergence in colour linked to land-use affinity.

## Next steps

- Investigate the plantation forest x frugivore interaction more deeply (e.g., which frugivore lineages are most affected).
- Explore spatial patterns (biome-specific effects).
- Sensitivity analyses: effect of abundance threshold, alternative colour metrics (VolumeUVS).
- Consider whether the A3 proportion score could be refined (e.g., relative to total records per species rather than total dataset size).
