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

## Analysis 1: Colour ~ Land-use

### Question

Does plumage colouration differ across land-use types after accounting for body mass, habitat, biome, and trophic niche?

### Models

Bayesian GLMMs via brms/CmdStan. For each colour response:

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
- All three responses fitted with **lognormal** family.
- Priors: Normal(0, 2) on intercept, Normal(0, 0.5) on fixed effects, Exponential(2) on SD and sigma.
- 4 chains, 3000 iterations (1500 warmup), adapt_delta = 0.90.

A **minimal** model (drops Habitat interaction) was also fitted for LOO-CV comparison.

### Model comparison

Reduced (with Habitat interaction) preferred over Minimal for all three responses by LOO-CV.

### Convergence

All models converged (Rhat < 1.01, ESS > 400, 0 divergences), except dichrocooney minimal which was borderline (Rhat 1.011). Note: dichrocooney was initially fitted with Student-t family but had severe non-convergence (Rhat 1.5–1.7, ESS < 10); switching to lognormal resolved this.

### Key results

**Variance explained:** R² = 17.8% (meancolcooney), 16.0% (malecolcooney), 12.5% (dichrocooney).

**Land-use main effects:** No land-use category has a 95% CI excluding zero for any colour metric. Land use alone does not shift community colour.

**Body mass:** Consistent negative effect across all three metrics (estimate -0.08 to -0.09 for mean/male colour, -0.017 for dichromatism). Larger birds are less colourful.

**Trophic niche:** Strongest predictor.
- Nectarivores and frugivores are significantly more colourful (+0.23 to +0.35).
- Frugivores, granivores, and invertivores are significantly less dichromatic (-0.19 to -0.21).

**Interactions:**
- Plantation forest × Frugivore: frugivores in plantations are significantly less colourful (mean -0.38, male -0.35). Strongest interaction in the dataset.
- Primary vegetation × Nectarivore: unexpectedly lower colour (-0.28 mean, -0.32 male).
- Plantation forest × Grassland/Human Modified habitats: consistently negative for colour.
- Woodland habitats show reduced dichromatism in plantation forest and primary vegetation.

## Analysis 2: Abundance ~ Colour × Land-use

### Question

Do more colourful species have different abundance across land-use types?

### Models

For each standardized colour predictor (z_meancolcooney, z_dichrocooney, z_malecolcooney):

```
abundance ~ z_colour * Predominant_simple + (1 | SSBS)
```

- Family: lognormal.
- Same priors and sampling settings as A1.

### Convergence

All three models converged.

### Key results

**Variance explained:** R² ~11% for all models (land use drives most abundance variation).

**Main colour effect:** Positive for all three metrics (~0.03 on log scale, 95% CI just excludes zero). More colourful species tend to be slightly more abundant overall.

**Colour × land-use interactions:**
- **Pasture:** positive interaction — colourful species are disproportionately abundant in pastures (credible for male colour and dichromatism).
- **Plantation forest:** negative interaction — colourful species are relatively less abundant (credible for all three metrics, estimates -0.05 to -0.06).
- **Primary vegetation:** negative interaction — similar to plantation forest (credible for all three, estimates -0.03 to -0.05).
- **Secondary:** no credible interaction.

## Analysis 3: Phylogenetic regression — Colour ~ land-use association

### Question

At the species level, and accounting for shared ancestry, do species that are more associated with particular land-use types differ in colouration? This complements A1 (site-level) by asking whether evolutionary lineages with stronger affinity for degraded or pristine habitats have distinct colour profiles.

### Species-level data (computed in `A3_01_setup.R`)

From the full raw dataset (216,834 records, 1,703 species), for each species and each level of `Predominant_simple`, a **land-use proportion score** is computed as:

```
proportion = (records of species i in land-use j) / (total records in dataset)
```

This gives 5 predictor variables per species: `prop_Cropland`, `prop_Pasture`, `prop_Plantation_forest`, `prop_Primary_vegetation`, `prop_Secondary`. Some species may have 0 for certain land uses.

Colour variables are species-level traits (same value per record of a given species), so the first non-NA value per species is used.

### Models

Bayesian phylogenetic regression via brms/CmdStan. For each colour response:

```
colour ~ prop_Cropland + prop_Pasture + prop_Plantation_forest +
         prop_Primary_vegetation + prop_Secondary +
         (1 | gr(phylo, cov = A))
```

- `A`: phylogenetic covariance matrix from `ape::vcv.phylo()` (correlation form).
- Tree trimmed to species present in both data and tree using `ape::drop.tip()`.
- Family: lognormal (all three responses).
- Priors: Normal(0, 2) on intercept, Normal(0, 1) on fixed effects (wider than A1/A2 because predictors are proportions on a different scale), Exponential(2) on SD and sigma.
- 4 chains, 4000 iterations (2000 warmup), adapt_delta = 0.95 (phylogenetic models need more careful sampling).

### Status

Scripts written (`A3_01_setup.R`, `A3_02_fit_models.R`, `A3_03_diagnostics_summary.R`). Awaiting phylogenetic tree file to be placed in `data/` — update `tree_file` variable in `A3_01_setup.R` with the filename.

## Summary interpretation (A1 + A2)

Land use does not directly shift community colour in a simple way. Effects are mediated by trophic ecology and habitat context:

1. Trophic niche is a much stronger predictor of colour than land use.
2. Plantation forests appear to filter against colourful frugivores specifically.
3. Colourful species are slightly more abundant overall, but this advantage disappears in forested land uses (plantation + primary) and is amplified in pastures.
4. Larger-bodied birds are consistently less colourful.

## Next steps

- Run A3 once phylogenetic tree is provided.
- Investigate the plantation forest × frugivore interaction more deeply.
- Explore spatial patterns (biome-specific effects).
- Sensitivity analyses: effect of abundance threshold, alternative colour metrics (VolumeUVS).
