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

- `Predominant_simple`: land-use category (reference = Primary vegetation). Levels: Cropland, Secondary, Pasture, Plantation forest.
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

**Land-use main effects (relative to Primary vegetation):** No land-use category has a 95% CI excluding zero for any colour metric. Land use alone does not shift community colour.

**Body mass:** Consistent negative effect across all three metrics (estimate -0.112 mean colour, -0.123 male colour, -0.024 dichromatism). Larger birds are less colourful. The negative effect is weaker in cropland and pasture (positive interactions with body mass: Pasture +0.055 mean, +0.061 male, credible; Cropland +0.029 mean, +0.031 male, credible).

**Trophic niche (in primary vegetation):**
- Frugivores are significantly more colourful (+0.18 mean colour), and significantly less dichromatic (-0.25).
- Invertivores are significantly less dichromatic (-0.37), and males are less colourful (-0.18).
- Omnivores are significantly less dichromatic (-0.27).

**Trophic interactions (relative to Primary vegetation):**
- Pasture x Frugivore: frugivores in pasture are significantly more colourful than in primary vegetation (+0.30 mean, +0.34 male, credible).
- Pasture x Nectarivore: nectarivores in pasture are significantly more colourful (+0.36 mean, +0.34 male, credible).
- Plantation forest x Frugivore: frugivores in plantations trend less colourful (-0.25 mean, borderline -- 95% CI just includes zero).
- Plantation forest x Nectarivore: nectarivores in plantations are significantly more dichromatic (+0.33, credible).
- Plantation forest x Invertivore/Omnivore: both guilds are significantly more dichromatic in plantations (+0.21 and +0.25 respectively, credible).

**Habitat interactions (relative to Primary vegetation):**
- Cropland x Human Modified: species from human-modified habitats are significantly more colourful in croplands (+0.22 mean, +0.25 male, credible).
- Plantation forest x Grassland: grassland species are significantly less colourful in plantations (-0.25 mean, -0.25 male, credible).
- Plantation forest x Wetland: wetland species are significantly less colourful in plantations (-0.24 mean, -0.32 male, credible).
- Cropland x Grassland: grassland species are significantly more dichromatic in cropland (+0.20, credible).
- Plantation forest x Riverine: riverine species are significantly more dichromatic in plantations (+0.35, credible).
- Cropland x Woodland: woodland species are significantly more dichromatic in cropland (+0.14, credible).

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

All three models converged (Rhat < 1.01, ESS > 400, 0 divergences), except malecolcooney which was borderline (Rhat 1.010).

### Key results

**Variance explained:** R² ~11% for all models (land use drives most abundance variation).

**Main colour effect (in Primary vegetation):** Near zero and not credible for any metric (mean: -0.003, dichro: -0.018, male: -0.010). In primary vegetation, colour does not predict abundance.

**Land-use main effects:** All modified land uses have credibly higher abundance than primary vegetation (Cropland +0.16, Pasture +0.13, Plantation forest +0.35, Secondary +0.17, all credible).

**Colour x land-use interactions (relative to Primary vegetation):**
- **Cropland:** positive interaction -- colourful species are disproportionately more abundant in croplands relative to primary vegetation (credible for all three metrics: mean +0.034, dichro +0.049, male +0.039).
- **Pasture:** strongest positive interaction -- colourful species are disproportionately more abundant in pastures (credible for all three: mean +0.063, dichro +0.090, male +0.084).
- **Plantation forest:** near-zero interaction -- no colour-abundance advantage in plantations (mean -0.023 borderline, dichro +0.000 not credible, male -0.020 borderline).
- **Secondary:** positive interaction -- colourful species relatively more abundant (credible for dichro +0.027 and male +0.028; mean +0.024 borderline).

## Analysis 3: Phylogenetic regression -- Colour ~ land-use association

### Question

At the species level, and accounting for shared ancestry, do species that are more associated with particular land-use types differ in colouration? This complements A1 (site-level) by asking whether evolutionary lineages with stronger affinity for degraded or pristine habitats have distinct colour profiles.

### Species-level data (computed in `A3_01_setup.R`)

From the full raw dataset (216,834 records, 1,703 species), for each species and each level of `Predominant_simple`, a **within-species land-use proportion** is computed as:

```
proportion = (records of species i in land-use j) / (total records of species i)
```

These are relative proportions that sum to 1 per species, removing the confound with overall species prevalence (dividing by total dataset size would conflate "how common is this species" with "which land uses does it prefer"). Because the 5 proportions sum to 1, Primary vegetation is dropped as the reference level (matching A1) to avoid perfect multicollinearity. This gives 4 predictor variables per species: `prop_Cropland`, `prop_Plantation_forest`, `prop_Pasture`, `prop_Secondary`, each interpretable as the fraction of that species' records occurring in that land-use type relative to Primary vegetation. VIFs for these predictors are 3.4-4.1 (acceptable).

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
colour ~ prop_Cropland + prop_Plantation_forest +
         prop_Pasture + prop_Secondary +
         (1 | gr(phylo, cov = A))
```

- Primary vegetation is the reference level (implicit intercept).
- `A`: phylogenetic covariance matrix from `ape::vcv.phylo()` (correlation form).
- `(1 | gr(phylo, cov = A))`: phylogenetic random effect accounting for shared ancestry.
- Tree trimmed to species present in both data and tree using `ape::drop.tip()`.
- Family: lognormal (all three responses).
- Priors: Normal(0, 2) on intercept, Normal(0, 1) on fixed effects (wider than A1/A2 because predictors are proportions on a different scale), Exponential(2) on SD and sigma.
- 4 chains, 4,000 iterations (2,000 warmup), adapt_delta = 0.95, max_treedepth = 12, threading enabled (4 threads per chain).
- N = 1,254 species per model (after excluding species with missing colour data).

Scripts: `A3_01_setup.R`, `A3_02_fit_models.R`, `A3_03_diagnostics_summary.R`.

### Convergence

All three models converged: 0 divergences, max Rhat <= 1.003, min bulk ESS >= 1,418, min tail ESS >= 2,669.

### Key results

**Variance explained:** R² = 88.8% (meancolcooney, 95% CI: 84.9-92.1%), 86.4% (malecolcooney, 82.1-90.2%), 47.1% (dichrocooney, 37.5-57.3%). The high R² is driven almost entirely by the phylogenetic random effect (shared ancestry).

**Land-use proportion effects (relative to Primary vegetation):**
- **Pasture:** credible negative effect on mean colour (-0.15, 95% CI: -0.26 to -0.04) and male colour (-0.15, 95% CI: -0.27 to -0.03). Species more associated with pastures are less colourful than those associated with primary vegetation, after accounting for phylogeny.
- **Cropland:** positive trend for mean colour (+0.10, CI: -0.04 to 0.23) and male colour (+0.09, CI: -0.07 to 0.24), but 95% CIs include zero.
- **Plantation forest:** near zero (+0.02 mean, +0.02 male), CIs include zero.
- **Secondary:** near zero (+0.04 mean, +0.03 male), CIs include zero.
- **Dichromatism:** no land-use proportion is a credible predictor (all estimates within -0.05 to +0.01, CIs include zero).

**Phylogenetic signal:** The phylogenetic SD is large relative to the residual: sd(phylo) = 0.52 (meancolcooney), 0.56 (malecolcooney), 0.26 (dichrocooney), compared to sigma = 0.14, 0.17, 0.17 respectively. This confirms strong phylogenetic conservatism in plumage colouration, especially for overall colourfulness. Dichromatism shows weaker (but still substantial) phylogenetic signal.

## Summary interpretation (A1 + A2 + A3)

All analyses use Primary vegetation as the reference level. Land use does not directly shift community colour in a simple way. Effects are mediated by trophic ecology and habitat context:

1. **Trophic niche is a much stronger predictor of colour than land use** (A1). Frugivores are the most colourful guild in primary vegetation, while invertivores and omnivores are the least dichromatic.
2. **Pastures favour colourful frugivores and nectarivores** (A1). The Pasture x Frugivore and Pasture x Nectarivore interactions are the strongest in the dataset: these guilds are significantly more colourful in pasture than in primary vegetation.
3. **Plantation forests disadvantage grassland and wetland species** (A1). Grassland and wetland species are significantly less colourful in plantations relative to primary vegetation. Plantation forest x Frugivore also trends negative (borderline).
4. **Colour does not predict abundance in primary vegetation, but does in modified land uses** (A2). In primary vegetation, the colour-abundance relationship is near zero. In croplands, pastures, and secondary vegetation, colourful species are disproportionately more abundant. The advantage is strongest in pasture. Plantation forests show no colour-abundance advantage.
5. **Larger-bodied birds are consistently less colourful** (A1), with the negative body mass effect weakened in croplands and pastures.
6. **Plumage colour is strongly phylogenetically conserved** (A3). R² from phylogeny alone is 87-89% for colourfulness and 47% for dichromatism.
7. **Pasture-associated species are less colourful after accounting for phylogeny** (A3). This is the only land-use effect that survives phylogenetic correction relative to primary vegetation, and it applies to overall and male colourfulness but not to dichromatism. Together with A2's finding that colourful species are disproportionately abundant in pastures at the community level, this suggests a sorting process: while colourful individuals gain an abundance advantage in pastures (A2), the lineages evolutionarily associated with open/pastoral habitats are inherently less colourful (A3).

## Next steps

- Investigate the plantation forest x frugivore interaction more deeply (e.g., which frugivore lineages are most affected).
- Explore spatial patterns (biome-specific effects).
- Sensitivity analyses: effect of abundance threshold, alternative colour metrics (VolumeUVS).
