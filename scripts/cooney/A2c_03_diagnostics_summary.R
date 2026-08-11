# ============================================================
# A2c_03_diagnostics_summary.R
# Convergence, fixed effects, and Bayesian R2 for A2c models.
# Outputs: results/A2c_convergence.csv, A2c_fixed_effects.csv,
#          results/A2c_r2_summary.csv
# ============================================================

library(brms)
library(posterior)
source("scripts/00_config.R")
library(tidyverse)

load("fits/A2c_all_fits.RData")   # `fits` (named list, one per colour var)
dir.create("results/cooney", showWarnings = FALSE)

# Prior-only biome x land-use cells, derived from the fitted data.
load("fits/A2c_model_setup.RData")   # `dat`
EMPTY_CELLS <- empty_cell_params(dat)
cat("Prior-only empty cells:",
    if (length(EMPTY_CELLS)) paste(EMPTY_CELLS, collapse = ", ") else "none", "\n")

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
write_csv(conv, "results/cooney/A2c_convergence.csv")
cat("Convergence:\n"); print(conv)

fixef_tbl <- map_dfr(names(fits), function(cv) {
  fx <- as.data.frame(fixef(fits[[cv]]))
  fx$parameter <- rownames(fx); fx$model <- cv
  fx$prior_only_empty_cell <- fx$parameter %in% EMPTY_CELLS
  as_tibble(fx)
}) %>% select(model, parameter, Estimate, Est.Error, Q2.5, Q97.5, prior_only_empty_cell)
write_csv(fixef_tbl, "results/cooney/A2c_fixed_effects.csv")

# Variance components (group-level SDs) — phylo vs site/study vs residual.
varcomp <- map_dfr(names(fits), function(cv) {
  vc <- VarCorr(fits[[cv]])
  grp <- map_dfr(names(vc), function(g)
    tibble(group = g, sd = vc[[g]]$sd[, "Estimate"]))
  sig <- summary(fits[[cv]])$spec_pars["sigma", "Estimate"]
  bind_rows(grp, tibble(group = "sigma", sd = sig)) %>% mutate(model = cv)
}) %>% select(model, group, sd)
write_csv(varcomp, "results/cooney/A2c_variance_components.csv")
cat("\nVariance components:\n"); print(varcomp)

r2 <- map_dfr(names(fits), function(cv) {
  r <- bayes_R2(fits[[cv]])
  tibble(model = cv, R2 = r[, "Estimate"], Q2.5 = r[, "Q2.5"], Q97.5 = r[, "Q97.5"])
})
write_csv(r2, "results/cooney/A2c_r2_summary.csv")
cat("\nBayesian R2:\n"); print(r2)

cat("\nA2c diagnostics complete. Tables in results/.\n")
