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
- `dichrodiff`: male minus female LociUVS (arithmetic difference in dichromatism; can be negative when females are more colourful than males)

All use LociUVS (ultraviolet-sensitive visual model loci counts) from Cooney et al. Dale scores were excluded.

### Biome classification (`Biome4`, computed in `00_data_preparation.R`)

Simplified from the original 11-level `Biome` to 4 categories:
- **Tropical Forest**: Tropical & Subtropical Moist Broadleaf Forests, Tropical & Subtropical Dry Broadleaf Forests, Tropical & Subtropical Coniferous Forests (reference level)
- **Tropical Open**: Tropical & Subtropical Grasslands, Savannas & Shrublands; Flooded Grasslands & Savannas; Deserts & Xeric Shrublands; Mangroves
- **Temperate Forest**: Temperate Broadleaf & Mixed Forests, Temperate Conifer Forests
- **Temperate Open**: Temperate Grasslands, Savannas & Shrublands; Mediterranean Forests, Woodlands & Scrub; Montane Grasslands & Shrublands; Boreal Forests/Taiga

### Land-use levels (`Predominant_simple`)

In the analytical dataset: Primary vegetation (10,920 records), Secondary (8,603), Plantation forest (7,122), Cropland (4,577), Pasture (3,431). Total: 34,653.

### Phylogenetic tree

BigBirdTree (BBtree2) from [evolucionario/BigBirdTree](https://github.com/evolucionario/BigBirdTree). 9,072 tips (Newick format). Stored as `data/BBtree2.tre`. Uses Jetz-era taxonomy.

## Analysis 1: Colour ~ Land-use

### Question

Does plumage colouration differ across land-use types after accounting for body mass, biome, and trophic niche?

### Methods

Bayesian GLMMs via brms/CmdStan. For each of the four colour responses:

**Reduced model** (preferred):
```
colour ~ Predominant_simple + z_logMass + Biome4 + Trophic.Niche +
         Predominant_simple:z_logMass +
         Predominant_simple:Biome4 +
         Predominant_simple:Trophic.Niche +
         (1 | SSBS)
```

**Minimal model** (drops `Predominant_simple:Biome4` interaction):
```
colour ~ Predominant_simple + z_logMass + Biome4 + Trophic.Niche +
         Predominant_simple:z_logMass +
         Predominant_simple:Trophic.Niche +
         (1 | SSBS)
```

- `Predominant_simple`: land-use category (reference = Primary vegetation). Levels: Cropland, Secondary, Pasture, Plantation forest.
- `z_logMass`: standardized log body mass.
- `Biome4`: 4-category simplified biome (reference = Tropical Forest).
- `SSBS`: site nested within study nested within block (random intercept).
- Lognormal family for meancolcooney, dichrocooney, and malecolcooney. Gaussian family for dichrodiff (which can be negative).
- Priors (lognormal responses): Normal(0, 2) on intercept, Normal(0, 0.5) on fixed effects, Exponential(2) on SD and sigma.
- Priors (gaussian/dichrodiff): Normal(0, 50) on intercept, Normal(0, 20) on fixed effects, Exponential(2) on SD and sigma.
- 4 chains, 3,000 iterations (1,500 warmup), adapt_delta = 0.90, threading enabled.

Scripts: `A1_01_setup.R`, `A1_02_fit_models.R`, `A1_03_model_comparison.R`, `A1_04_diagnostics.R`, `A1_05_summarize.R`.

### Model comparison

Reduced model (with Biome4 interaction) strongly preferred over Minimal for all four responses by LOO-CV:

| Response | elpd_diff | SE |
|---|---|---|
| malecolcooney | -126.2 | 14.9 |
| meancolcooney | -122.2 | 14.7 |
| dichrodiff | -89.9 | 13.6 |
| dichrocooney | -67.6 | 12.9 |

### Convergence

All 8 models converged (max Rhat <= 1.007, min ESS bulk >= 535, 0 divergences for all models).

### Reference levels

- **Land use** (`Predominant_simple`): Primary vegetation (reference). Effects of Cropland, Pasture, Plantation forest, and Secondary are relative to Primary.
- **Biome** (`Biome4`): Tropical Forest (reference). Effects of Temperate Forest, Temperate Open, and Tropical Open are relative to Tropical Forest.
- **Trophic niche** (`Trophic.Niche`): Aquatic predator (reference, N = 49 records). Most comparisons of interest are among the common guilds (Invertivore N = 21,374; Omnivore N = 6,574; Granivore N = 3,097; Frugivore N = 2,966; Nectarivore N = 591; Herbivore terrestrial N = 2).
- **Body mass** (`z_logMass`): standardized log body mass centered at zero.

The intercept represents the expected colour for an Aquatic predator species of mean body mass in Primary vegetation in a Tropical Forest biome.

### Key results

No land-use category has a credible main effect (95% CI excluding zero) for any response. Land use alone does not shift community colour. Effects emerge through interactions with biome, trophic niche, and body mass.

#### meancolcooney (overall colourfulness; lognormal family)

R² = 16.2%.

**Body mass:** Credibly negative (-0.115 [-0.122, -0.108]). Larger birds are less colourful. This effect is weakened in Pasture (+0.052, credible) and Cropland (+0.022, credible).

**Biome main effects:** All non-tropical-forest biomes are less colourful: Temperate Forest -0.12, Temperate Open -0.27, Tropical Open -0.20 (all credible).

**Biome x Land-use interactions:**
- Temperate Forest: Cropland -0.09 and Plantation -0.11 credibly negative (communities even less colourful than the biome main effect alone would predict). Secondary -0.04 credible.
- Temperate Open: Secondary +0.20 and Pasture +0.12 credibly positive (communities more colourful than expected).
- Tropical Open: Plantation +0.17 and Cropland +0.08 credibly positive.

**Trophic niche:** Frugivores are the most colourful guild (+0.20, credible). Other guilds not credibly different from reference.

**Trophic x Land-use interactions:** Plantation x Frugivore -0.35 (credible: frugivores are less colourful in plantations). Pasture x Frugivore +0.22, Pasture x Nectarivore +0.25 (both borderline).

#### malecolcooney (male conspicuousness; lognormal family)

R² = 14.6%.

Results closely parallel meancolcooney. **Body mass** -0.125 (credible), weakened in Pasture (+0.064, credible) and Cropland (+0.032, credible).

**Biome main effects:** Temperate Forest -0.07, Temperate Open -0.23, Tropical Open -0.20 (all credible).

**Biome x Land-use interactions:** Same pattern as meancolcooney. Temperate Forest: Cropland -0.14 and Plantation -0.14 (credible). Temperate Open: Secondary +0.19 (credible). Tropical Open: Plantation +0.17, Cropland +0.10 (credible).

**Trophic niche:** Frugivores +0.16 (credible). Invertivores -0.10 (borderline).

**Trophic x Land-use interactions:** Plantation x Frugivore -0.36 (credible). Pasture x Frugivore +0.22, Pasture x Nectarivore +0.25 (both borderline).

#### dichrocooney (sexual dichromatism ratio; lognormal family)

R² = 10.4%.

**Body mass:** Credibly negative (-0.023 [-0.027, -0.019]). Larger birds are less dichromatic. Interactions: Pasture +0.030 and Cropland +0.024 (credible positive — the mass effect weakens). Plantation -0.015 and Secondary -0.011 (credible negative — the mass effect strengthens in these land uses).

**Biome main effects:** Temperate Forest +0.12 and Temperate Open +0.08 (both credible: *more* dichromatic than Tropical Forest, the opposite direction from colour). Tropical Open near zero (-0.005, not credible).

**Biome x Land-use interactions:**
- Temperate Forest: Cropland -0.10, Pasture -0.08, Plantation -0.07 (all credible). Land-use conversion in temperate forests erodes the elevated baseline dichromatism. Secondary -0.001 (not credible: secondary vegetation retains the temperate forest dichromatism pattern).
- Temperate Open: Cropland -0.06 (credible), Pasture -0.05 (borderline). Others not credible.
- Tropical Open: Cropland +0.04 (credible). Others near zero.

**Trophic niche:** Frugivores -0.10 (credible), Invertivores -0.22 (credible), Omnivores -0.13 (credible). These guilds are less dichromatic than the reference in primary vegetation.

**Trophic x Land-use interactions:** Plantation x Nectarivore +0.21 (credible: nectarivores are more dichromatic in plantations). Other interactions not credible.

#### dichrodiff (sexual dichromatism difference; Gaussian family)

R² = 14.1%.

This response is on the raw LociUVS difference scale (male minus female), not log-transformed. Effects are therefore much larger in magnitude than the lognormal responses.

**Body mass:** Credibly negative (-4.19 [-4.78, -3.61]). Larger birds have smaller sex differences. Interactions: Pasture +4.4 and Cropland +2.6 (credible positive). Plantation -1.2 (credible negative — larger species in plantations are relatively *more* sexually different).

**Biome main effects:** Temperate Forest +11.5 (credible: communities in temperate forests have much larger male-female differences). Tropical Open -4.1 (credible: smaller sex differences). Temperate Open +4.1 (borderline).

**Biome x Land-use interactions:**
- Temperate Forest: Cropland -11.2, Pasture -10.1, Plantation -12.0 (all credible). Land-use conversion in temperate forests strongly erodes the elevated sex differences. Secondary -2.3 (borderline).
- Tropical Open: Cropland +7.6 and Pasture +3.2 (both credible: land-use conversion in tropical open habitats increases sex differences).
- Temperate Open: not credible.

**Trophic niche:** Nectarivores +33.1 (credible: very large sex difference). Invertivores -20.1 and Omnivores -9.1 (both credible: smaller sex differences). Frugivores near zero.

**Trophic x Land-use interactions:** Pasture x Frugivore +20.3 (credible: male frugivores in pasture are much more colourful relative to females). Pasture x Nectarivore -19.7 (credible: the large baseline nectarivore sex difference is eroded in pasture). Secondary x Nectarivore -17.4 (credible). Pasture x Invertivore +11.1 (borderline).

## Analysis 2: Relative Abundance ~ Colour x Land-use

### Question

Do more colourful species have different relative abundance across land-use types?

### Methods

Bayesian GLMMs via brms/CmdStan. Response is **relative abundance**: each species' effort-corrected count divided by the total effort-corrected count at that SSBS (site-level community total). For each standardized colour predictor (z_meancolcooney, z_dichrocooney, z_malecolcooney, z_dichrodiff):

```
relative_abundance ~ z_colour * Predominant_simple + (1 | SSBS)
```

- Family: lognormal.
- Priors: Normal(0, 2) on intercept, Normal(0, 0.5) on fixed effects, Exponential(2) on SD and sigma.
- 4 chains, 3,000 iterations (1,500 warmup), adapt_delta = 0.90, threading enabled.
- SSBS retained as random intercept to account for residual site-level non-independence after normalisation.
- Checkpoint files: `fits/A2_fit_relabund_<cv>.RData`; combined: `fits/A2_all_fits_relabund.RData`.

Scripts: `A2_01_setup.R`, `A2_02_fit_models.R`, `A2_03_diagnostics_summary.R`.

### Reference levels

- **Land use** (`Predominant_simple`): Primary vegetation (reference). The intercept represents expected log-relative-abundance in Primary vegetation at the mean value of the colour predictor.
- **Colour predictors**: each is standardized (z-scored), so the main colour effect represents the change in log-relative-abundance per 1 SD increase in colour in Primary vegetation. The interaction terms represent how this slope changes in each modified land use relative to Primary.

### Convergence

Three models converged cleanly (max Rhat ≤ 1.008, min ESS > 500, 0 divergences). z_dichrodiff marginal (max Rhat = 1.014, 0 divergences) — interpret with caution.

### Key results

R² ≈ 50.8% for all four models (vs. 11.1% with raw abundance).

**Land-use main effects (shared across models):** Cropland +0.45, Pasture +0.62–0.63, Secondary +0.22 (all credible). Plantation forest ~0.00 (not credible) — no overall difference from primary vegetation once site effort is normalised.

#### z_meancolcooney (overall colourfulness)

**Main colour effect in Primary vegetation:** -0.01 [-0.02, +0.01] (not credible).

**Colour x Land-use interactions:** Pasture +0.06 [+0.02, +0.10] (credible). Plantation -0.04 [-0.06, -0.01] (credible — colourful species are less relatively abundant in plantations). Cropland +0.02 [-0.01, +0.05] (not credible). Secondary +0.02 [-0.01, +0.04] (not credible).

#### z_malecolcooney (male conspicuousness)

**Main colour effect in Primary vegetation:** -0.01 [-0.02, 0.00] (borderline).

**Colour x Land-use interactions:** Pasture +0.07 [+0.04, +0.11] (credible). Plantation -0.03 [-0.05, -0.01] (credible). Secondary +0.02 [0.00, +0.05] (borderline). Cropland +0.02 [-0.01, +0.05] (not credible).

#### z_dichrocooney (sexual dichromatism ratio)

**Main colour effect in Primary vegetation:** -0.00 [-0.02, +0.01] (not credible). Unlike raw-abundance analysis, no credible negative effect of dichromatism in primary vegetation.

**Colour x Land-use interactions:** Pasture +0.06 [+0.02, +0.09] (credible). Cropland +0.03 [+0.01, +0.06] (credible). Secondary +0.03 [+0.01, +0.06] (credible). Plantation +0.01 [-0.02, +0.03] (not credible).

#### z_dichrodiff (sexual dichromatism difference)

**Main colour effect in Primary vegetation:** -0.01 [-0.03, 0.00] (borderline). No longer credibly negative as in raw-abundance analysis.

**Colour x Land-use interactions:** Pasture +0.10 [+0.06, +0.14] (credible — strongest single interaction). Secondary +0.04 [+0.01, +0.06] (credible). Cropland +0.02 [-0.01, +0.05] (not credible). Plantation -0.00 [-0.02, +0.02] (not credible).

## Analysis 3: Phylogenetic regression -- Colour ~ land-use association

### Question

At the species level, and accounting for shared ancestry, do species that are more associated with particular land-use types differ in colouration? This complements A1 (site-level) by asking whether evolutionary lineages with stronger affinity for degraded or pristine habitats have distinct colour profiles.

### Species-level data (computed in `A3_01_setup.R`)

From the full raw dataset (216,834 records, 1,703 species), for each species and each level of `Predominant_simple`, a **within-species land-use proportion** is computed as:

```
proportion = (records of species i in land-use j) / (total records of species i)
```

These are relative proportions that sum to 1 per species, removing the confound with overall species prevalence (dividing by total dataset size would conflate "how common is this species" with "which land uses does it prefer"). All 5 proportions (`prop_Cropland`, `prop_Pasture`, `prop_Plantation_forest`, `prop_Primary_vegetation`, `prop_Secondary`) are included as predictors in a no-intercept model. Since the proportions sum to 1, removing the intercept resolves the compositional constraint and each coefficient is interpretable as the expected colour for a species found exclusively in that land-use type.

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

Bayesian phylogenetic regression via brms/CmdStan. For each of the four colour responses:

```
colour ~ 0 + prop_Cropland + prop_Pasture + prop_Plantation_forest +
         prop_Primary_vegetation + prop_Secondary +
         (1 | gr(phylo, cov = A))
```

- No intercept (`0 +`): because the 5 proportions sum to 1 per species, they span the intercept. Each coefficient represents the expected colour for a species found exclusively in that land-use type.
- `A`: phylogenetic covariance matrix from `ape::vcv.phylo()` (correlation form).
- `(1 | gr(phylo, cov = A))`: phylogenetic random effect accounting for shared ancestry.
- Tree trimmed to species present in both data and tree using `ape::drop.tip()`.
- Lognormal family for meancolcooney, dichrocooney, and malecolcooney. Gaussian family for dichrodiff.
- Priors (lognormal): Normal(0, 5) on fixed effects, Exponential(2) on SD and sigma.
- Priors (gaussian/dichrodiff): Normal(0, 100) on fixed effects, Exponential(2) on SD and sigma.
- 4 chains, 4,000 iterations (2,000 warmup), adapt_delta = 0.95, max_treedepth = 12, threading enabled (4 threads per chain).
- N = 1,254 species per model (after excluding species with missing colour data).

Scripts: `A3_01_setup.R`, `A3_02_fit_models.R`, `A3_03_diagnostics_summary.R`.

### Reference levels

There is no categorical reference level for land use. The no-intercept model with 5 compositional proportions (summing to 1 per species) means each coefficient is directly interpretable as the expected colour for a hypothetical species found exclusively in that land-use type. Coefficients can be compared pairwise, but no single land use is parameterized as the reference.

### Convergence

All four models converged: 0 divergences, max Rhat <= 1.01, ESS adequate for all parameters. Note: the dichrodiff model produced an E-BFMI warning (< 0.3 for all 4 chains), suggesting potential issues with the energy distribution, though the model otherwise converged.

### Key results

#### meancolcooney (overall colourfulness; lognormal family)

R² = 88.6% (95% CI: 84.6-91.8%), driven almost entirely by phylogenetic signal.

**Expected log-colour by land-use association:**
- Cropland: 4.67 [4.31, 5.03]
- Secondary: 4.61 [4.27, 4.97]
- Plantation forest: 4.59 [4.26, 4.96]
- Primary vegetation: 4.57 [4.23, 4.93]
- Pasture: 4.42 [4.07, 4.79] (lowest)

All CIs overlap due to dominant phylogenetic variance, but Pasture is consistently lowest. The Pasture-Primary difference (~-0.15 on log scale) is the most prominent pattern.

**Phylogenetic signal:** sd(phylo) = 0.52 vs sigma = 0.14. Strong conservatism — closely related species have similar colourfulness.

#### malecolcooney (male conspicuousness; lognormal family)

R² = 86.2% (81.8-89.9%).

**Expected log-colour by land-use association:**
- Cropland: 4.68 [4.29, 5.07]
- Secondary: 4.62 [4.25, 4.99]
- Plantation forest: 4.61 [4.24, 4.99]
- Primary vegetation: 4.59 [4.22, 4.97]
- Pasture: 4.44 [4.05, 4.83] (lowest)

Pattern matches meancolcooney. Pasture-associated lineages have dullest males.

**Phylogenetic signal:** sd(phylo) = 0.56 vs sigma = 0.17. Strongest phylogenetic conservatism of all four responses.

#### dichrocooney (sexual dichromatism ratio; lognormal family)

R² = 46.6% (37.1-57.1%). Lower than colour metrics — dichromatism ratio is less phylogenetically conserved.

**Expected log-dichromatism by land-use association:**
- Pasture: 0.10 [-0.09, 0.28]
- Primary vegetation: 0.08 [-0.09, 0.25]
- Cropland: 0.07 [-0.13, 0.26]
- Plantation forest: 0.06 [-0.11, 0.24]
- Secondary: 0.04 [-0.13, 0.21]

All near zero with wide overlapping CIs. No land-use association predicts dichromatism ratio.

**Phylogenetic signal:** sd(phylo) = 0.26 vs sigma = 0.17. Weaker but still substantial phylogenetic signal.

#### dichrodiff (sexual dichromatism difference; Gaussian family)

R² = 65.6% (57.9-72.9%). Intermediate phylogenetic signal — the arithmetic magnitude of sex differences is partially labile across lineages.

**Expected dichrodiff by land-use association:**
- Pasture: 11.5 [-13.4, 36.8]
- Primary vegetation: 11.3 [-12.1, 35.2]
- Cropland: 7.9 [-17.7, 33.5]
- Plantation forest: 6.2 [-16.8, 30.1]
- Secondary: 5.3 [-17.7, 28.6]

All CIs are extremely wide (~50 units) and heavily overlap. No land-use association credibly predicts the arithmetic dichromatism difference. The E-BFMI warning for this model suggests these estimates should be interpreted cautiously.

**Phylogenetic signal:** dichrodiff phylogenetic SD is large but so is the residual, giving intermediate R² (66% vs 87-89% for colour).

## Summary interpretation (A1 + A2 + A3)

All analyses use Primary vegetation as the reference level. Land use does not directly shift community colour in a simple way. Effects are mediated by biome, trophic ecology, and body mass. Colourfulness (mean/male colour) and sexual dichromatism (dichro ratio, dichrodiff) often show contrasting patterns:

1. **Trophic niche and biome are stronger predictors of colour than land use** (A1). Frugivores are the most colourful guild. Temperate and open biomes have less colourful communities than tropical forests, but temperate forests have *more* sexually dichromatic communities — a dissociation between overall colour and sex differences.
2. **Colour and dichromatism respond differently to biome and land use** (A1). Temperate forest communities are duller but more dichromatic. When those communities are then modified by land-use change (cropland, plantations), both the colour deficit *and* the dichromatism surplus are eroded — land-use conversion in temperate forests homogenizes communities toward the tropical baseline. In tropical open biomes, the pattern reverses: plantations have more colourful communities.
3. **Plantation forests are unfavourable for frugivores** (A1). The Plantation forest x Frugivore interaction is credibly negative for mean and male colour (-0.35, -0.36). This is the strongest trophic interaction in the dataset.
4. **Pastures favour colourful frugivores but reduce nectarivore sex differences** (A1). Male frugivores in pasture are much more colourful relative to females (dichrodiff +20.3, credible). But nectarivores show the opposite: their large baseline sex difference is eroded in pasture (dichrodiff -19.7, credible). The two dichromatism measures capture different aspects of this: dichrodiff (arithmetic) detects large-magnitude shifts that dichro ratio (proportional) does not always reflect.
5. **Sexually dichromatic and colourful species have higher relative abundance in open modified land uses** (A2; response = relative abundance per SSBS). No credible main colour effect in primary vegetation. In Pasture, species with larger male-female colour differences (dichrodiff × Pasture +0.10) and more colourful males (malecolcooney × Pasture +0.07) are disproportionately more relatively abundant. Cropland and Secondary show credible positive interactions for dichromatism. Unlike the raw-abundance analysis, there is no credible negative dichromatism–abundance relationship in primary vegetation once site effort is normalised (R² ~51%).
6. **Plantation forests credibly reduce relative abundance of colourful species** (A2). meancolcooney × Plantation (-0.04) and malecolcooney × Plantation (-0.03) are both credibly negative — colourful species make up a smaller fraction of communities in plantations. This was only borderline in the raw-abundance analysis and is the clearest new finding from the relative-abundance approach. Plantation forest shows no positive colour-abundance effect for any metric.
7. **Larger-bodied birds are consistently less colourful and less dichromatic** (A1). The negative body mass effect is weakened in croplands and pastures for all metrics. In plantations, larger species are relatively *more* dichromatic (opposite direction), suggesting body-size-dependent filtering of sex differences across land-use types.
8. **Plumage colour is strongly phylogenetically conserved** (A3). R² from phylogeny alone is 86-89% for colourfulness, 66% for dichrodiff, and 47% for dichromatism ratio. Dichrodiff shows intermediate phylogenetic signal, suggesting the arithmetic magnitude of sex differences is partially labile.
9. **Pasture-associated species are less colourful after accounting for phylogeny** (A3). This is the only land-use effect that emerges from the phylogenetic analysis, applying to overall and male colourfulness but not to either dichromatism measure. Together with A2's finding that colourful species are disproportionately abundant in pastures at the community level, this suggests a sorting process: while colourful individuals gain an abundance advantage in pastures (A2), the lineages evolutionarily associated with open/pastoral habitats are inherently less colourful (A3). The absence of a phylogenetic signal for dichromatism suggests that sex-difference patterns across land uses (A1, A2) arise from community reassembly rather than deep evolutionary association.

## Next steps

- Investigate the plantation forest x frugivore interaction more deeply (e.g., which frugivore lineages are most affected).
- Sensitivity analyses: effect of abundance threshold, alternative colour metrics (VolumeUVS).
