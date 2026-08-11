# ============================================================
# B2d_03_diagnostics_summary.R
# Convergence, fixed effects, variance components, Bayesian R2
# for the B2d models. Mirrors A2d_03. Outputs:
#   results/B2d_{convergence,fixed_effects,variance_components,r2_summary}.csv
# ============================================================

library(brms)
library(posterior)
source("scripts/00_config.R")
library(tidyverse)

load("fits/B2d_all_fits.RData")     # `fits` (named list, one per colour var)
load("fits/B2d_model_setup.RData")  # `dat_final`
dir.create("results/dale", showWarnings = FALSE)

# Prior-only biome x land-use cells, derived from the fitted data.
# ref_lu is Secondary here: primary vegetation is absorbed into the
# paired-difference response, so it is not a factor level.
EMPTY_PARAMS <- empty_cell_params(dat_final, ref_lu = "Secondary")
cat("Empty (prior-only) interaction cells flagged:",
    if (length(EMPTY_PARAMS)) paste(EMPTY_PARAMS, collapse = ", ") else "none", "\n")

conv <- map_dfr(names(fits), function(cv) {
  fit <- fits[[cv]]; s <- summary(fit)
  rh <- c(s$fixed[, "Rhat"], s$spec_pars[, "Rhat"],
          unlist(lapply(s$random, function(x) x[, "Rhat"])))
  eb <- c(s$fixed[, "Bulk_ESS"], s$spec_pars[, "Bulk_ESS"],
          unlist(lapply(s$random, function(x) x[, "Bulk_ESS"])))
  np <- nuts_params(fit)
  tibble(model = cv,
         max_rhat = round(max(rh, na.rm = TRUE), 4),
         min_ess_bulk = round(min(eb, na.rm = TRUE)),
         n_divergent = sum(np$Value[np$Parameter == "divergent__"]))
})
write_csv(conv, "results/dale/B2d_convergence.csv")
cat("Convergence:\n"); print(conv)

fixef_tbl <- map_dfr(names(fits), function(cv) {
  fx <- as.data.frame(fixef(fits[[cv]]))
  fx$parameter <- rownames(fx); fx$model <- cv
  fx$prior_only_empty_cell <- fx$parameter %in% EMPTY_PARAMS
  as_tibble(fx)
}) %>% select(model, parameter, Estimate, Est.Error, Q2.5, Q97.5, prior_only_empty_cell)
write_csv(fixef_tbl, "results/dale/B2d_fixed_effects.csv")

varcomp <- map_dfr(names(fits), function(cv) {
  vc <- VarCorr(fits[[cv]])
  grp <- map_dfr(names(vc), function(g)
    tibble(group = g, sd = vc[[g]]$sd[, "Estimate"]))
  sig <- summary(fits[[cv]])$spec_pars["sigma", "Estimate"]
  bind_rows(grp, tibble(group = "sigma", sd = sig)) %>% mutate(model = cv)
}) %>% select(model, group, sd)
write_csv(varcomp, "results/dale/B2d_variance_components.csv")
cat("\nVariance components:\n"); print(varcomp)

r2 <- map_dfr(names(fits), function(cv) {
  r <- bayes_R2(fits[[cv]])
  tibble(model = cv, R2 = r[, "Estimate"], Q2.5 = r[, "Q2.5"], Q97.5 = r[, "Q97.5"])
})
write_csv(r2, "results/dale/B2d_r2_summary.csv")
cat("\nBayesian R2:\n"); print(r2)

cat("\nB2d diagnostics complete. Tables in results/.\n")
