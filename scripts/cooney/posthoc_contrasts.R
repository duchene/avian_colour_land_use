# ============================================================
# posthoc_contrasts.R
# All pairwise contrasts between factor levels, from the posterior.
#
# The Bayesian analogue of a Tukey post-hoc test. For each factor we
# build an estimated marginal mean (EMM) per level by averaging the
# linear predictor over a balanced reference grid of the OTHER
# factors, then take every pairwise difference of the posterior
# draws. No multiplicity correction is applied or needed: these are
# posterior probabilities of a difference, not repeated frequentist
# tests, so the Tukey/Bonferroni family-wise machinery does not
# apply. The averaging is balanced (equal weight per cell), matching
# the emmeans default, so the contrasts describe the model rather
# than the sampling design's cell sizes.
#
# Random effects are excluded (re_formula = NA), so contrasts are
# population-level. Body mass is held at z_logMass = 0 (the mean,
# ~20.3 g) and the colour predictors in A2c/A2d at their mean (z=0).
#
# Reported per contrast:
#   estimate, lo, hi   difference of EMMs on the modelled (link) scale
#   pd                 posterior probability of the stated direction
#   credible           95% CrI excludes zero
#   pct_change         (exp(diff)-1)*100, lognormal responses only
#   effect_size        difference / SD of the response on the link
#                      scale, a Cohen's-d-like standardised measure
#
# Run: Rscript scripts/cooney/posthoc_contrasts.R
# ============================================================

source("scripts/cooney/posthoc_helpers.R")
suppressMessages({ library(tidyr); library(stringr); library(readr) })

dir.create("results/cooney", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/pub/cooney", showWarnings = FALSE, recursive = TRUE)

# One factor at a time: EMMs over a balanced grid of the others, then
# every pairwise difference. Interactions are opened up in Table 2
# (posthoc_interactions.R), which reports the same contrasts computed
# separately within each biome.
contrasts_for <- function(fit, model, response, specs) {
  sdl <- link_sd(fit)
  lognorm <- fit$family$family == "lognormal"
  map_dfr(names(specs), function(fac) {
    emm <- emm_draws(fit, specs[[fac]], fac)
    pairwise(emm, sdl, lognorm, model = model, response = response, factor = fac)
  })
}

all_contrasts <- list()

# ------------------------------------------------------------
# A3b: community colour (covariate variant, the reported one)
# ------------------------------------------------------------
message("A3b ...")
lu <- LU5; bio <- BIOME3; tro <- TROPHIC6

grid_a1b <- make_grid(Predominant_simple = lu, Biome4 = bio,
                      Trophic.Niche = tro, z_logMass = 0)

for (resp in c("malecolcooney", "dichrodiff")) {
  f <- paste0("fits/A3b_fit_", resp, "_covariate.RData")
  load(f)
  all_contrasts[[paste0("A3b_", resp)]] <- contrasts_for(
    fit, "A3b", resp,
    list(Predominant_simple = grid_a1b, Biome4 = grid_a1b, Trophic.Niche = grid_a1b))
  rm(fit); gc(verbose = FALSE)
  message("  ", resp, " done")
}

# ------------------------------------------------------------
# A2c: relative abundance
# ------------------------------------------------------------
message("A2c ...")
for (cv in c("z_malecolcooney", "z_dichrodiff")) {
  f <- paste0("fits/A2c_fit_", cv, ".RData")
  load(f)
  g <- make_grid(Predominant_simple = lu, Biome4 = bio)
  g[[cv]] <- 0
  all_contrasts[[paste0("A2c_", cv)]] <- contrasts_for(
    fit, "A2c", cv, list(Predominant_simple = g, Biome4 = g))
  rm(fit); gc(verbose = FALSE)
  message("  ", cv, " done")
}

# ------------------------------------------------------------
# A2d: abundance change. Primary vegetation is absorbed into the
# response, so land use carries four levels with Secondary as
# reference and every contrast is a difference of differences.
# ------------------------------------------------------------
message("A2d ...")
lu_d <- LU4
for (cv in c("z_malecolcooney", "z_dichrodiff")) {
  f <- paste0("fits/A2d_fit_", cv, ".RData")
  load(f)
  g <- make_grid(Predominant_simple = lu_d, Biome4 = bio)
  g[[cv]] <- 0
  all_contrasts[[paste0("A2d_", cv)]] <- contrasts_for(
    fit, "A2d", cv, list(Predominant_simple = g, Biome4 = g))
  rm(fit); gc(verbose = FALSE)
  message("  ", cv, " done")
}

# ------------------------------------------------------------
# A1: species-level compositional model. There is no intercept and
# the five proportions sum to 1, so each coefficient IS the marginal
# mean for a species found exclusively in that land use. Contrast
# the coefficient draws directly rather than through a grid.
# ------------------------------------------------------------
message("A1 ...")
load("fits/A1_all_fits.RData")
lu_prop <- c("prop_Primary_vegetation", "prop_Secondary",
             "prop_Plantation_forest", "prop_Cropland", "prop_Pasture")
lu_lab  <- str_replace_all(str_remove(lu_prop, "^prop_"), "_", " ")

for (resp in c("malecolcooney", "dichrodiff")) {
  fit <- fits_A1[[resp]]
  d <- posterior::as_draws_matrix(fit)
  emm <- d[, paste0("b_", lu_prop), drop = FALSE]
  colnames(emm) <- lu_lab
  all_contrasts[[paste0("A1_", resp)]] <- pairwise(
    emm, link_sd(fit), fit$family$family == "lognormal",
    model = "A1", response = resp, factor = "Predominant_simple")
  message("  ", resp, " done")
}
rm(fits_A1); gc(verbose = FALSE)

# ------------------------------------------------------------
# Write out
# ------------------------------------------------------------
tbl <- bind_rows(all_contrasts) %>%
  mutate(response = recode(response,
           malecolcooney = "Male colourfulness", dichrodiff = "Sexual dichromatism",
           z_malecolcooney = "Male colourfulness", z_dichrodiff = "Sexual dichromatism"),
         factor = recode(factor,
           Predominant_simple = "Land use", Biome4 = "Biome",
           Trophic.Niche = "Trophic niche")) %>%
  select(model, response, factor, contrast, estimate, lo, hi, pd,
         credible, pct_change, effect_size)

tbl <- tbl %>% mutate(across(all_of(c("contrast")), relabel_levels))
write_csv(tbl, "results/cooney/posthoc_contrasts.csv")
write_csv(tbl, "figures/pub/cooney/Table_1_posthoc_contrasts.csv")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("CREDIBLE PAIRWISE CONTRASTS (95% CrI excludes zero)\n")
cat(strrep("=", 70), "\n", sep = "")
cred <- tbl %>% filter(credible)
if (nrow(cred)) {
  cred %>%
    mutate(across(c(estimate, lo, hi, effect_size), ~ round(.x, 4)),
           pd = round(pd, 3)) %>%
    as.data.frame() %>% print(row.names = FALSE)
} else cat("none\n")

cat("\nTotal contrasts:", nrow(tbl), "| credible:", sum(tbl$credible), "\n")
cat("Written to results/cooney/posthoc_contrasts.csv",
    "and figures/pub/cooney/Table_1_posthoc_contrasts.csv\n")
