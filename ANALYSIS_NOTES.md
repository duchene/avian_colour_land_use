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

Simplified from the original `Biome` field. The mapping **as actually coded in `00_data_preparation.R`** (verified against the data; only 11 biomes occur):
- **Tropical Forest** (reference): Tropical & Subtropical Moist Broadleaf Forests; Dry Broadleaf Forests; Coniferous Forests
- **Tropical Open**: Tropical & Subtropical Grasslands, Savannas & Shrublands; Mangroves
- **Temperate Forest**: Temperate Broadleaf & Mixed Forests; Temperate Conifer Forests; **Mediterranean Forests, Woodlands & Scrub**
- **Temperate Open** (`else` catch-all): Temperate Grasslands, Savannas & Shrublands; Montane Grasslands & Shrublands; Tundra (plus any unlisted biome — but none of Deserts & Xeric Shrublands, Flooded Grasslands, or Boreal Forests/Taiga occur in the analytical data)

**Note (unresolved coding choice):** Mediterranean Forests, Woodlands & Scrub (690 records, ~2%) is placed in **Temperate Forest** by the code, and every analysis (A1, A1b, A2c, A2d) uses this coding. Mediterranean sclerophyll is genuinely borderline (its WWF name spans "Forests" and "Woodlands & Scrub"). If a reviewer prefers Mediterranean → Temperate Open, `Biome4` must be recomputed in `00_data_preparation.R` and the affected analyses rerun.

### Land-use levels (`Predominant_simple`)

In the analytical dataset: Primary vegetation (10,920 records), Secondary (8,603), Plantation forest (7,122), Cropland (4,577), Pasture (3,431). Total: 34,653.

### Phylogenetic tree

BigBirdTree (BBtree2) from [evolucionario/BigBirdTree](https://github.com/evolucionario/BigBirdTree). 9,072 tips (Newick format). Stored as `data/BBtree2.tre`. Uses Jetz-era taxonomy.

## Analysis 1: Colour ~ Land-use

> **Superseded.** This section documents the original Analysis 1 (a single `(1 | SSBS)`
> random intercept), whose scripts and result files were removed. It is replaced by two
> models with the full PREDICTS random structure, both under "Extended analyses" below:
> **A1b** (primary, biome × land-use only) and **A1a** (the refit of this model, adding
> land-use × trophic niche and land-use × body mass interactions). The A1a refit
> reproduces the findings below (for example plantation × frugivore) under the correct
> random structure. The numbers in this section are from the original random structure.

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

Scripts: removed. Superseded by `A1a_01_setup.R` / `A1a_02_fit_models.R` / `A1a_03_diagnostics_summary.R` (see Extended analyses › A1a below).

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

## Analysis 2b: Paired-difference — Abundance change ~ Colour × Land-use

### Question

Do colourful or sexually dichromatic species gain or lose relative abundance when habitats are converted from primary vegetation to other land uses? This complements A2 by using a paired-difference design that (a) includes the full dataset (presences and absences), (b) explicitly anchors each species' abundance to its own primary-vegetation baseline within each study, and (c) removes uninformative zeros (species absent from both primary and modified sites).

### Motivation

A2 used only the "present" dataset (non-zero abundance records) with a lognormal family, which excluded absences entirely. Because ~82% of records in the full primary-veg-filtered dataset are zeros, A2 could only ask "among species that are present, does colour predict relative dominance?" and could not detect species being filtered out of communities altogether. A paired-difference approach sidesteps the zero-inflation problem: differences of exactly zero (species absent from both primary and modified) are excluded as uninformative, while species that disappeared (negative difference) or appeared (positive difference) are retained as ecologically meaningful signals.

### Data preparation (computed in `A2b_01_setup.R`)

Starting from the full raw dataset (216,834 records, 1,703 species):

1. **Filter to qualifying studies:** studies containing Primary vegetation plus at least one other land-use type → 143,071 records.
2. **Compute relative abundance per SSBS:** each species' effort-corrected measurement divided by the total effort-corrected measurement at that SSBS.
3. **Compute primary-vegetation baseline:** for each species × study combination, the mean relative abundance across all primary-vegetation SSBS within that study (2,336 species × study combinations; median number of primary sites per combination = 48).
4. **Paired differences:** for each non-primary SSBS record, subtract the species' primary-vegetation baseline → 87,605 paired differences.
5. **Exclude zero differences:** remove records where the difference is exactly zero (species absent from both primary and modified sites) → **66,831 non-zero differences** (54,529 negative, 12,302 positive).
6. **Colour data availability:** 62,975 records with non-NA malecolcooney and dichrodiff.

### Sample structure

- 38 studies, 1,441 SSBS, 1,434 species.
- Land-use breakdown: Plantation forest (26,857), Secondary (23,816), Cropland (12,269), Pasture (3,889).
- Response median: -0.004 (slight net loss relative to primary, as expected).
- 82% of non-zero differences are negative (species lost abundance relative to primary).

### Methods

Bayesian GLMMs via brms/CmdStan. For each of two colour predictors (z_malecolcooney, z_dichrodiff):

```
diff_abund ~ z_colour * Predominant_simple + (1 | SSBS) + (1 | Reference)
```

- `diff_abund`: species relative abundance at modified SSBS minus mean relative abundance at primary-vegetation SSBS within the same study.
- `Predominant_simple`: land-use category of the modified SSBS (Cropland, Pasture, Plantation forest, Secondary). Primary vegetation is absent — it is built into the response as the baseline.
- `z_colour`: standardized colour predictor (z_malecolcooney or z_dichrodiff).
- `(1 | SSBS)`: random intercept for site, accounting for within-site non-independence.
- `(1 | Reference)`: random intercept for study, accounting for between-study variation in difference magnitudes.
- Gaussian family (response can be negative, zero was excluded, positive).
- Priors: Normal(0, 0.5) on intercept, Normal(0, 0.2) on fixed effects, Exponential(5) on SD and sigma.
- 4 chains, 2,000 iterations (1,000 warmup), adapt_delta = 0.90, max_treedepth = 10, threading enabled.

Scripts: `A2b_01_setup.R`, `A2b_02_fit_models.R`, `A2b_03_diagnostics_summary.R`.

### Interpretation of parameters

- **Intercept:** mean abundance change (modified minus primary) for the reference land-use level at mean colour.
- **Land-use main effects:** how each land use's mean abundance change differs from the reference, at mean colour. Negative values indicate the land use causes greater abundance loss on average.
- **Colour main effect:** does colour predict the magnitude of abundance change in the reference land use? A positive coefficient means more colourful species lose less (or gain more) abundance.
- **Colour × land-use interactions:** does the colour effect on abundance change differ by land use?

### Convergence

All four models (2 interaction, 2 main-effects) converged: 0 divergences, max Rhat = 1.00, min ESS bulk >= 2466, min ESS tail >= 1846.

### Key results

R² is effectively zero for all four models (0.01–0.05%), indicating that colour explains negligible variance in paired abundance differences.

#### z_malecolcooney (male conspicuousness)

**Main-effects model:** Main colour effect -0.0003 [-0.0007, +0.0002] (not credible). No land-use main effects credible.

**Interaction model:** No credible effects. All coefficients on the order of 0.0001–0.0007 with CIs spanning zero.

**Random effects:** sd(Reference) = 0.0003, sd(SSBS) = 0.0002, sigma = 0.060. Residual variance dominates entirely.

#### z_dichrodiff (sexual dichromatism difference)

**Main-effects model:** Main colour effect **+0.0011 [+0.0007, +0.0016]** (credible). Species with larger male-female colour differences lose slightly less relative abundance across all modified land uses. No land-use main effects credible.

**Interaction model:** Main colour effect +0.0013 [+0.0002, +0.0025] (credible, but wider CI as it applies only to the reference land use). No interactions credible: the dichrodiff effect does not differ meaningfully across land-use types.

**Random effects:** sd(Reference) = 0.0003, sd(SSBS) = 0.0002, sigma = 0.060. Same as malecolcooney.

#### Interpretation

The paired-difference approach, while conceptually sound for including absences, yields extremely low signal-to-noise. Most species × site differences are near-zero, and the colour signal is swamped by residual noise. The one credible finding — more dichromatic species losing slightly less abundance — is consistent in direction with A2's finding that dichromatic species are disproportionately abundant in modified land uses, but the magnitude is negligible (~0.1 percentage point per SD of dichrodiff). The absence of any malecolcooney effect contrasts with A2's credible Pasture and Plantation interactions, suggesting those patterns are detectable only when conditioning on presence.

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

## Extended analyses: phylogenetic models with full PREDICTS random structure (A1a, A1b, A2c, A2d)

**Status:** A1b finalised as a **no-phylogeny** model (a feasibility fit on 2026-06-23 showed the phylogenetic term is not identifiable when colour is the response — see "A1b — note on complexity" below); A1b and A2c fitted 2026-06-25 (results in their subsections below — A2c overturns A2's colour × land-use findings); A2d fitted 2026-07-02, refit at higher resolution 2026-07-03 (results below; converged — Rhat ≤ 1.015, elevated only on the global intercept). Proposed by collaborator. These extend A1/A2/A2b by adding the canonical PREDICTS nested random structure, biome as an interacting factor, and — for the abundance models A2c/A2d only — a phylogenetic random effect. The aim is to test whether the earlier land-use patterns survive control for the full study/block/site hierarchy and (for the abundance models) shared ancestry. A three-way-interaction sensitivity check (**A2e/A2f**, 2026-07-08) confirms the 2-way abundance models are the appropriate reported models — see their subsection. **A1a** (the refit of the original Analysis 1 with A1b's full random structure, adding land-use × trophic niche and land-use × body mass interactions) was fitted 2026-07-20 and is favoured over A1b by LOO for all four responses (see the A1a subsection).

### Design decisions common to A1b, A2c, A2d

- **Random structure (canonical PREDICTS):** `(1 | SS) + (1 | SSB) + (1 | SSBS)` — study / block-within-study / site-within-block. This replaces the single `(1 | SSBS)` of A1/A2 and the `(1|SSBS) + (1|Reference)` of A2b. Confirmed in `present.csv`: 59 studies (SS) → 361 blocks (SSB) → 3,716 sites (SSBS), cleanly nested. **`SS`** (the PREDICTS study code, 59 levels) is used as the study grouping, not `Reference` (the citation, 54 levels — a few citations contain two studies). Because SSBS is globally unique and nested, `(1|SS) + (1|SSB) + (1|SSBS)` is the nested structure directly — no `/` nesting syntax is needed.
- **Phylogenetic random effect (A2c / A2d only — dropped from A1b; see the A1b complexity note):** `(1 | gr(phylo, cov = A))`, with `A` = `vcv.phylo()` correlation matrix from BBtree2, reusing A3's 3-stage taxonomic matching (jetz_sp → Best_guess_binomial underscored → synonym table `results/A3_synonym_table.csv`). Species absent from BBtree2 are dropped (no grafting). In `present.csv` this retains 1,348 of 1,681 species and 26,951 of 34,653 records (77.8%). Coverage is re-checked per dataset at setup.
- **Interactions: 2-way only — no 3-way terms anywhere.** A `colour × biome × land-use` design is reduced to all main effects plus the three pairwise interactions; the three-way `colour:Biome:land-use` term is omitted.
- **Sparse-cell collapse.** Biome × land-use cells too sparse to support an interaction estimate (threshold ≈ <50 records / <15 species) are collapsed. In `present.csv` the only such cell is **Temperate Open × Plantation forest (42 records, 12 species)**, which is reclassified to Temperate Open × Secondary (the nearest non-primary woody land use). *Take clear note:* this is the only collapse applied, it affects 42 of 34,653 records, and it is re-verified for the A2d dataset at setup. Alternative if preferred: drop those 42 records rather than merge.
- **Families** unchanged from parent analyses: lognormal for malecolcooney / meancolcooney / dichrocooney and for relative abundance; Gaussian for dichrodiff and for the A2d paired difference.
- **Colour / dichromatism metrics:** primary runs use **malecolcooney** (colour) and **dichrodiff** (dichromatism), matching A2b; meancolcooney and dichrocooney are extensions if the primary models fit cleanly.
- **Expected phylogeny–trait confound.** Colour, trophic niche, and body mass are all strongly phylogenetically conserved (A3: colour R² 86–89%). Where colour/niche/mass enter as *predictors* (A2c, A2d), the phylogenetic random effect competes with them for the same species-level variance; their fixed-effect estimates are expected to shrink and should be read as effects *beyond phylogeny*. Biome × land-use effects vary *within* species (a species spans biomes/land-uses) and are not subject to this confound.

### A1b — Community colour/dichromatism ~ biome × land-use (no phylogeny)

Extends **A1**. Dataset: `present.csv` (occurrence records of present species) — **all** of them: A1b has no phylogenetic term, so unlike A2c/A2d no species are dropped for tree-matching. Each model is fitted to the 33,210 records with a non-missing response (of 34,653 in present.csv; vs the ~27k phylo-matched subset used by A2c/A2d, which drop species absent from BBtree2). The response is continuous community colour — this is a community-weighted-trait model, **not** an occupancy/presence-absence (0/1) model. Two variants per response:

```
(base)      colour ~ Biome4 * Predominant_simple
                     + (1|SS) + (1|SSB) + (1|SSBS)
(covariate) colour ~ Biome4 * Predominant_simple + Trophic.Niche + z_logMass
                     + (1|SS) + (1|SSB) + (1|SSBS)
```

`Biome4 * Predominant_simple` expands to both main effects plus their 2-way interaction. `Trophic.Niche` and `z_logMass` enter as **main effects only** — A1's `land-use × trophic` and `land-use × mass` interactions are intentionally dropped (decision: keep only the biome × land-use interaction). Consequently A1b will *not* recover A1's Plantation × Frugivore (-0.35/-0.36) or mass × land-use results. Reference levels as in A1 (Primary vegetation; Tropical Forest). Relative to A1, the only changes are the restructured fixed effects (a single biome × land-use interaction) and the full PREDICTS random hierarchy `(1|SS) + (1|SSB) + (1|SSBS)` in place of A1's single `(1|SSBS)`.

#### A1b — note on complexity (why there is no phylogenetic term)

A1b was originally specified *with* a phylogenetic random effect `(1 | gr(phylo, cov = A))`, matching the collaborator's proposal and the A2c/A2d structure. A feasibility fit (malecolcooney, base model; 2 chains × 1,000 iterations; 2026-06-23) showed this specification is **not identifiable**, and the term was dropped.

**Why it fails.** Plumage colour is a *species-level trait* — identical across all of a species' occurrence records (A3 treats it the same way). The phylogenetic term has one level per species (1,242 levels in this dataset), so it behaves as a per-species random intercept fitted to a response with *zero within-species variance*. The phylogenetic effect therefore absorbs essentially all between-species variation, while the site/study random effects (`SS`, `SSB`, `SSBS`), the residual `sigma`, and every biome × land-use fixed effect collapse toward zero. Diagnostics confirmed structural degeneracy rather than slow mixing: sd(phylo) = 0.35 while sd(SS) = sd(SSB) = sd(SSBS) = sigma ≈ 0; Rhat 1.8–3.0 and bulk ESS 2–7 across nearly all parameters (despite 0 divergences); wall time ≈ 69 min for the reduced 2-chain run.

**Contrast.** A1 worked without a phylogenetic term because the species-constant response was absorbed by the fixed effects + site random effects + a substantial residual (R² ≈ 15%). A3 places a phylogenetic term on colour legitimately because it is fitted at the *species level* (one row per species), where the covariance matrix `A` constrains the phylogenetic effect and `sigma` captures the iid residual. The degeneracy is specific to putting a per-species phylogenetic intercept on a within-species-constant response at the *record* level.

**Consequence for the study.** The phylogenetic angle on colour ↔ environment is handled at the species level by **A3** (and could be extended there to include biome if desired). A1b is therefore the community-level (record-weighted) model *without* phylogeny, differing from A3 by being record-level / occurrence-weighted with categorical biome × land-use and site random effects. The phylogenetic random effect is retained only in the abundance models A2c/A2d, where the response (abundance) genuinely varies within species. This split — community level without phylogeny, species level with phylogeny — is the only identifiable way to ask both questions.

#### A1b — results (fitted 2026-06-25)

All 8 models (4 responses × {base, +trophic/mass}) converged: 0 divergences, max Rhat ≤ 1.017, min bulk ESS ≥ 467. Tables: `results/A1b_{convergence,fixed_effects,r2_summary}.csv` (the empty Temperate Open × Plantation cell is flagged in the fixed-effects table).

Bayesian R²: meancolcooney 10.2% (base) / 16.8% (+covariate); malecolcooney 7.7% / 15.0%; dichrocooney 6.2% / 11.0%; dichrodiff 6.8% / 15.6%. Adding trophic niche + body mass roughly **doubles** R² over the biome × land-use-only base, and the +covariate R² closely matches A1 (16.2 / 14.6 / 10.4 / 14.1%) — i.e. A1b reproduces A1's explanatory power despite the restructured fixed effects (single biome × land-use interaction) and the full PREDICTS random hierarchy.

### A1a — Community colour/dichromatism ~ biome × land-use + land-use × trophic + land-use × mass (no phylogeny)

Extends **A1b** and refits the original Analysis 1 under the full PREDICTS structure.
Same data (`present.csv`), sparse-cell collapse, priors, sampler and random structure as
A1b, plus the two interaction families A1b omits:

```
colour ~ Biome4 * Predominant_simple + Trophic.Niche + z_logMass
         + Predominant_simple:z_logMass + Predominant_simple:Trophic.Niche
         + (1|SS) + (1|SSB) + (1|SSBS)
```

A1a is a nested superset of the A1b covariate model, so LOO adjudicates whether the
guild- and size-specific interactions add predictive value. Scripts: `A1a_01_setup.R` /
`A1a_02_fit_models.R` / `A1a_03_diagnostics_summary.R` (the last also runs the LOO
comparison and needs the A1b covariate fits). Tables:
`results/A1a_{convergence,fixed_effects,r2_summary}.csv` and `results/A1a_vs_A1b_loo.csv`.
Publication Figure 3 derives the trophic and mass interactions by filtering
`A1a_fixed_effects.csv`.

#### A1a — results (fitted 2026-07-20)

All four models converged: 0 divergences, max Rhat ≤ 1.01, min bulk ESS ≥ 846. Bayesian
R² rises over A1b in every response: meancolcooney 16.8 → 18.7%, malecolcooney 15.0 →
17.1%, dichrocooney 11.0 → 13.2%, dichrodiff 15.6 → 17.5%.

**LOO (A1a full vs A1b covariate), A1a favoured for all four responses.** elpd_diff
(A1a minus A1b): meancolcooney +294 (SE 33), malecolcooney +322 (SE 33), dichrocooney
+419 (SE 35), dichrodiff +343 (SE 34), each about 9 to 12 SE. The interaction set earns
its place.

**Trophic × land-use (the headline, all survive the refit under the full structure):**
- Plantation × Frugivore: meancolcooney -0.29 [-0.55, -0.03], malecolcooney -0.31
  [-0.58, -0.03] (credible). Frugivores much duller in plantation forest.
- dichrodiff: Pasture × Frugivore +20.7 [5.9, 35.1], Pasture × Nectarivore -18.1
  [-33.4, -3.4], Secondary × Nectarivore -15.7 [-28.6, -3.1] (all credible).
- dichrocooney: Plantation × Nectarivore +0.20 [0.01, 0.38] (credible).

**Mass × land-use:** the negative mass slope (larger birds duller and less dichromatic)
weakens in Pasture and Cropland (malecolcooney +0.059 and +0.026 per SD log mass, both
credible) and strengthens for dichromatism in Plantation and Secondary (dichrocooney
-0.017 and -0.011, dichrodiff Plantation -1.38, all credible).

**Biome × land-use:** A1a's biome × land-use terms differ from A1b's because they now
condition on the trophic and mass interactions, and some flip sign (for example
meancolcooney Temperate Forest × Plantation is +0.08 credible in A1b but not credible in
A1a). Biome × land-use is therefore reported from A1b, and A1a contributes the guild and
size interactions. Herbivore terrestrial (N=2) land-use interactions are prior-only and
flagged (`sparse_trophic`) in the fixed-effects table.

### A2c — Relative abundance ~ colour × biome × land-use (2-way; phylogenetic)

Extends **A2**. Dataset: `present.csv`, response = relative abundance per SSBS (lognormal), as in A2. Per standardized colour predictor:

```
relabund ~ z_colour + Biome4 + Predominant_simple
           + z_colour:Biome4 + z_colour:Predominant_simple + Biome4:Predominant_simple
           + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))
```

All three pairwise interactions; **no** `z_colour:Biome4:Predominant_simple`. Adds biome and the phylogenetic/PREDICTS structure to A2. Phylogeny–colour confound applies (see common decisions).

#### A2c — results (fitted 2026-06-25; z_malecolcooney, z_dichrodiff)

Both models converged: 0 divergences, max Rhat ≤ 1.005, min bulk ESS ≥ 1,003. R² = 51% for both (≈ A2; carried by the random structure, not colour). Variance components (near-identical across both): sd(phylo) = 0.91, sd(SS) = 0.95, sd(SSB) = 0.23, sd(SSBS) = 0.26, sigma = 0.67 — relative abundance carries strong phylogenetic **and** study-level structure. Tables: `results/A2c_{convergence,fixed_effects,variance_components,r2_summary}.csv`.

**Key finding — A2's colour × land-use effects do not survive.** Adding biome, phylogeny, and the full PREDICTS random structure erases every credible A2 colour × land-use interaction. These are well-estimated zeros (bulk ESS ~1,000, tight CIs), not low power:

- malecolcooney × Pasture: A2 +0.07 (credible) → A2c −0.02 [−0.06, 0.03] (ns)
- malecolcooney × Plantation: A2 −0.03 (credible) → A2c 0.00 [−0.03, 0.03] (ns)
- dichrodiff × Pasture: A2 +0.10 (credible; A2's strongest interaction) → A2c −0.01 [−0.06, 0.04] (ns)
- dichrodiff × Secondary: A2 +0.04 (credible) → A2c +0.02 [−0.01, 0.04] (ns)

What survives is a **colour × biome** effect: malecolcooney × Tropical Open = **−0.09 [−0.15, −0.03] (credible)** — more colourful males are *less* relatively abundant in tropical-open communities. (dichrodiff × Temperate Open = −0.09 [−0.18, 0.00] is borderline.) Colour main effects sit at ~0 (the phylogenetic term absorbs the conserved between-species colour variance). The colour × biome effect survives because biome varies *within* species and so is not confounded with the per-species phylogenetic term.

**Interpretation:** the A2 colour–abundance–land-use associations are confounded with biome and phylogeny. The only robust colour–abundance signal is biome-level (colourful species rarer in tropical-open communities), not land-use-level. This **supersedes Synthesis points 5–6** below.

### A2d — Paired-difference abundance change ~ colour × biome × land-use (2-way; phylogenetic)

Extends **A2b**. Dataset: the A2b paired-difference data (Gaussian; Primary vegetation is the baseline, so `Predominant_simple` has 4 levels: Cropland, Pasture, Plantation forest, Secondary). Per standardized colour predictor:

```
diff_abund ~ z_colour + Biome4 + Predominant_simple
             + z_colour:Biome4 + z_colour:Predominant_simple + Biome4:Predominant_simple
             + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))
```

Adds biome, the phylogenetic random effect, and the full PREDICTS hierarchy (A2b used `(1|SSBS) + (1|Reference)`; A2d switches the study grouping to the canonical `SS`/`SSB`/`SSBS`). A2b's R² was ≈ 0, so biome/phylogeny are unlikely to rescue strong signal — this is primarily a robustness check. Same phylogeny–colour confound as A2c.

#### A2d — results (fitted 2026-07-02; z_malecolcooney, z_dichrodiff)

Data: study unit = `SS` throughout (qualifying filter and primary-veg baseline computed per SS × species); 40 studies, 66,787 non-zero paired differences; final set 54,003 records / 1,167 tips (per-model NA subsetting → malecol n = 50,972). Data-driven collapse: Temperate Open × Plantation has 242 records here (>50) so it was **not** collapsed; Temperate Open × Cropland is **empty (0 records)** → its interaction coefficient is prior-only and flagged in the fixed-effects table. Tables: `results/A2d_{convergence,fixed_effects,variance_components,r2_summary}.csv`.

**Convergence** (refit 2026-07-03 at 4 chains × 3,500 iter, warmup 1,500, adapt_delta 0.95, max_treedepth 12): 0 divergences; min bulk ESS ≈ 615–635. Max Rhat = 1.015, elevated **only on the global Intercept** (ESS 615 — a benign intercept ↔ random-intercept trade-off with four random-intercept levels); every colour, biome, and land-use parameter has Rhat ≤ 1.01. Point estimates were essentially identical to the initial run, so the credible interactions below are **confirmed, not provisional**. (Initial marginal run: max Rhat 1.03, ESS ~340.)

**R² = 22.6% (both models)** — a large jump from A2b's ≈ 0%, driven almost entirely by the phylogenetic random effect: sd(phylo) = 0.072 > sigma = 0.050 > sd(SS) = 0.015 ≫ sd(SSB), sd(SSBS). Abundance *change* relative to primary carries phylogenetic (lineage-level) structure that A2b (no phylo) could not capture — but this is a lineage effect, **not** a colour effect.

**Colour effects are tiny.** On the paired-difference scale (sigma ≈ 0.05; response ~ ±0.06), a few colour × biome interactions are credible but minuscule: malecolcooney × Tropical Open +0.009 [0.004, 0.014]; malecolcooney × Temperate Forest +0.009 [0.003, 0.015]; dichrodiff × Tropical Open +0.011 [0.006, 0.015]; malecolcooney × Plantation +0.002 [0.000, 0.003] (barely). All colour main effects and all other colour × land-use interactions are not credible. A2b's one credible finding (dichrodiff main +0.0011) does not clearly replicate in A2d's interaction model (+0.0011 [−0.0017, 0.0038], ns).

**Interpretation:** consistent with A2b, colour barely predicts net abundance change; the credible colour × biome terms are negligible in magnitude. A2d's substantive addition is that abundance change is phylogenetically structured (a lineage-level pattern), independent of colour. Note the sign contrast with A2c on Tropical Open (A2c: colourful males *less* abundant in standing tropical-open communities, −0.09; A2d: colourful species lose slightly *less* abundance relative to baseline there, +0.009) — these are different responses (standing relative abundance vs change from primary baseline), so not directly contradictory.

### A2e / A2f — three-way interaction sensitivity checks (fitted 2026-07-08)

To justify reporting the 2-way models, the **full three-way** interaction `z_colour * Biome4 * Predominant_simple` was fitted for **malecolcooney** and compared to the reported 2-way model by LOO. **A2e** extends A2c (relative abundance, lognormal); **A2f** extends A2d (paired difference, Gaussian). Same data, priors, and sampler settings as the 2-way parents, so the comparison is fair. Both converged (0 divergences; A2e max Rhat 1.006, A2f 1.009). Tables: `results/A2e_* results/A2f_*` (incl. `_loo_compare.csv`); fits gitignored.

**LOO comparison (2-way vs 3-way):**
- **A2e (relative abundance):** the 2-way is preferred — the 3-way is credibly *worse* (elpd_diff = −10.9, SE 4.2; ≈ 2.6 SE).
- **A2f (paired difference):** the two are statistically indistinguishable (elpd_diff = −3.4, SE 4.5; < 1 SE) — the 3-way adds no predictive value, so parsimony favours the 2-way.
- Caveat: a handful of high Pareto-k observations (8–26 of ~27k–51k), so LOO is slightly optimistic, but the ordering is clear.

**Three-way term credibility:**
- **A2e:** *none* of the 12 three-way terms is credible (all CIs span zero; the Temperate Open cells are very wide, e.g. malecol × Temperate Open × Pasture = −0.24 [−0.52, 0.03]). The empty Temperate Open × Plantation cell is prior-only (flagged).
- **A2f:** 8 of 9 three-way terms not credible; a single boundary-credible term (malecol × Temperate Forest × Plantation = +0.008 [0.0003, 0.017]) — negligible on the ±0.06 difference scale, at the CI boundary, expected by chance across 9 tests, and *not* supported by LOO. The empty Temperate Open × Cropland cell is prior-only (flagged).

**Conclusion:** the three-way interaction is not favoured by LOO and yields no robust three-way effect, for either response. The **2-way models (A2c/A2d) are the appropriate reported models**; A2e/A2f stand as the documented sensitivity check. (Check used malecolcooney only; extend to dichrodiff if an all-metric statement is required.)

### Scripts and run status

`A1b_01_setup.R / _02_fit_models.R / _03_diagnostics_summary.R`; likewise `A2c_*` and `A2d_*`. The 3-way sensitivity checks are single scripts: `A2e_3way_sensitivity.R` and `A2f_3way_sensitivity.R` (each reuses its parent's setup RData and adds LOO). A2c/A2d setup scripts reuse A3's tree-matching block (`A3_01_setup.R`); A1b needs no tree matching (no phylogenetic term). **All three analyses have been fitted** (A1b + A2c: 2026-06-25; A2d: 2026-07-02, refit at higher resolution 2026-07-03). Model objects are gitignored (`fits/*.RData`) — regenerate via the scripts. Runtimes (this machine, 16 cores): the 8 no-phylo A1b models run in the low hours total; each phylogenetic production fit (A2c/A2d, 4 chains) takes ≈ 2–4 h.

**Outstanding:** optionally extend A2c/A2d to meancolcooney + dichrocooney (confirm the pattern across all four metrics). The manuscript (`draft_methods_results.txt`) and `README.md` now cover A1a, A1b, A2c, A2d, and A3, with A1a as Figure 3.

## B family: Dale-colour replication (B1a, B1b, B2c, B2d, B3)

**Status:** prepared 2026-07-28, not yet fitted. All 15 scripts written and all five
setups verified against their A counterparts.

The B family repeats every focal analysis with Dale et al. (2015) plumage scores in place
of the Cooney UVS metrics, as an independent-colour-source robustness check. Each B
analysis mirrors its A counterpart exactly (data, sparse-cell collapse, phylogenetic
matching, priors, sampler, and random structure), differing only in the colour columns
and the `B`-prefixed outputs and fits. Mapping: B1a from A1a, B1b from A1b, B2c from A2c,
B2d from A2d, B3 from A3. The B1a-vs-B1b LOO is the Dale analogue of A1a-vs-A1b.

### Dale metrics (added in 00_data_preparation.R)

- `meancoldale` = mean(Male_plumage_score_dale, Female_plumage_score_dale), lognormal
- `malecoldale` = Male_plumage_score_dale, lognormal
- `dichrodale` = Male_plumage_score_dale / Female_plumage_score_dale, lognormal (ratio)
- `dichrodiffdale` = Male_plumage_score_dale - Female_plumage_score_dale, gaussian (difference)

Dale scores are positive (male 34 to 74, female 34 to 72), the difference ranges -7.9 to
34, so the families match the Cooney set (lognormal for the three positive metrics,
Gaussian for the difference). Coverage is slightly better than Cooney: 1,679 species have
a Dale male score vs 1,530 for Cooney, so the B runs drop fewer records to colour-NA.
present.csv now carries the four Dale columns, added without changing any existing column
or row (verified: all 85 original columns byte-identical after regeneration).

### Metric sets (like-for-like with A)

- B1a, B1b, B3: all four Dale metrics.
- B2c, B2d: `z_malecoldale` and `z_dichrodiffdale`. Colour enters standardised, so the
  predictor scale is identical to the Cooney runs.

Priors, families, sampler settings, and the sparse-cell collapse are all carried over
unchanged from the A scripts.

### Setup verification (2026-07-28)

Setups reproduce the A pipelines. B1b and B1a use 34,653 records (SS 59, SSB 361, SSBS
3,716). B2c uses 27,720 records and 1,375 tips. B2d uses 54,003 records and 1,167 tips
across 40 studies, with Temperate Open × Plantation at 242 records and no collapse, as in
A2d. B3 matches 1,391 tips, with Dale colour present for about 1,387 of them, more than
A3's 1,254, reflecting Dale's broader coverage. Fits are pending. Runtimes mirror the A
family: the no-phylogeny B1a and B1b are cheap, each phylogenetic fit (B2c, B2d, B3)
takes a few hours.

## Summary interpretation (A1 + A2 + A2b + A3)

All analyses use Primary vegetation as the reference level. Land use does not directly shift community colour in a simple way. Effects are mediated by biome, trophic ecology, and body mass. Colourfulness (mean/male colour) and sexual dichromatism (dichro ratio, dichrodiff) often show contrasting patterns:

1. **Trophic niche and biome are stronger predictors of colour than land use** (A1). Frugivores are the most colourful guild. Temperate and open biomes have less colourful communities than tropical forests, but temperate forests have *more* sexually dichromatic communities — a dissociation between overall colour and sex differences.
2. **Colour and dichromatism respond differently to biome and land use** (A1). Temperate forest communities are duller but more dichromatic. When those communities are then modified by land-use change (cropland, plantations), both the colour deficit *and* the dichromatism surplus are eroded — land-use conversion in temperate forests homogenizes communities toward the tropical baseline. In tropical open biomes, the pattern reverses: plantations have more colourful communities.
3. **Plantation forests are unfavourable for frugivores** (A1). The Plantation forest x Frugivore interaction is credibly negative for mean and male colour (-0.35, -0.36). This is the strongest trophic interaction in the dataset.
4. **Pastures favour colourful frugivores but reduce nectarivore sex differences** (A1). Male frugivores in pasture are much more colourful relative to females (dichrodiff +20.3, credible). But nectarivores show the opposite: their large baseline sex difference is eroded in pasture (dichrodiff -19.7, credible). The two dichromatism measures capture different aspects of this: dichrodiff (arithmetic) detects large-magnitude shifts that dichro ratio (proportional) does not always reflect.
5. **Sexually dichromatic and colourful species have higher relative abundance in open modified land uses** (A2; response = relative abundance per SSBS). No credible main colour effect in primary vegetation. In Pasture, species with larger male-female colour differences (dichrodiff × Pasture +0.10) and more colourful males (malecolcooney × Pasture +0.07) are disproportionately more relatively abundant. Cropland and Secondary show credible positive interactions for dichromatism. Unlike the raw-abundance analysis, there is no credible negative dichromatism–abundance relationship in primary vegetation once site effort is normalised (R² ~51%). **[Superseded by A2c (2026-06-25): these colour × land-use interactions are confounded with biome and phylogeny and are no longer credible once both are included — see Extended analyses › A2c results.]**
6. **Plantation forests credibly reduce relative abundance of colourful species** (A2). meancolcooney × Plantation (-0.04) and malecolcooney × Plantation (-0.03) are both credibly negative — colourful species make up a smaller fraction of communities in plantations. This was only borderline in the raw-abundance analysis and is the clearest new finding from the relative-abundance approach. Plantation forest shows no positive colour-abundance effect for any metric. **[Superseded by A2c (2026-06-25): the malecolcooney × Plantation effect is not credible under the biome + phylogeny model — see Extended analyses › A2c results.]**
7. **Colour–abundance relationships operate among present species, not through wholesale filtering** (A2b). When absences are included via a paired-difference design (modified minus primary baseline), only a single weak effect survives: more dichromatic species lose marginally less abundance (+0.0011 per SD dichrodiff; credible but negligible R²). No male-colour effect is detectable. This indicates that the A2 patterns (colour predicting relative dominance) reflect reshuffling of community proportions among species that persist, rather than differential extirpation of colourful species.
8. **Larger-bodied birds are consistently less colourful and less dichromatic** (A1). The negative body mass effect is weakened in croplands and pastures for all metrics. In plantations, larger species are relatively *more* dichromatic (opposite direction), suggesting body-size-dependent filtering of sex differences across land-use types.
9. **Plumage colour is strongly phylogenetically conserved** (A3). R² from phylogeny alone is 86-89% for colourfulness, 66% for dichrodiff, and 47% for dichromatism ratio. Dichrodiff shows intermediate phylogenetic signal, suggesting the arithmetic magnitude of sex differences is partially labile.
10. **Pasture-associated species are less colourful after accounting for phylogeny** (A3). This is the only land-use effect that emerges from the phylogenetic analysis, applying to overall and male colourfulness but not to either dichromatism measure. Together with A2's finding that colourful species are disproportionately abundant in pastures at the community level, this suggests a sorting process: while colourful individuals gain an abundance advantage in pastures (A2), the lineages evolutionarily associated with open/pastoral habitats are inherently less colourful (A3). The absence of a phylogenetic signal for dichromatism suggests that sex-difference patterns across land uses (A1, A2) arise from community reassembly rather than deep evolutionary association.

## Next steps

- Investigate the plantation forest x frugivore interaction more deeply (e.g., which frugivore lineages are most affected).
- Sensitivity analyses: effect of abundance threshold, alternative colour metrics (VolumeUVS).
