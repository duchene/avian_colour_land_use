# ============================================================
# B2c_03_diagnostics_summary.R
# Convergence, fixed effects, variance components, and Bayesian R2
# for the B2c models. Mirrors A2c_03. Outputs:
#   results/B2c_convergence.csv, B2c_fixed_effects.csv,
#   results/B2c_variance_components.csv, B2c_r2_summary.csv
# ============================================================

library(brms)
library(posterior)
library(tidyverse)

load("fits/B2c_all_fits.RData")   # `fits` (named list, one per colour var)
dir.create("results/dale", showWarnings = FALSE)

EMPTY_CELL <- "Biome4TemperateOpen:Predominant_simplePlantationforest"

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
write_csv(conv, "results/dale/B2c_convergence.csv")
cat("Convergence:\n"); print(conv)

fixef_tbl <- map_dfr(names(fits), function(cv) {
  fx <- as.data.frame(fixef(fits[[cv]]))
  fx$parameter <- rownames(fx); fx$model <- cv
  fx$prior_only_empty_cell <- grepl(EMPTY_CELL, fx$parameter, fixed = TRUE)
  as_tibble(fx)
}) %>% select(model, parameter, Estimate, Est.Error, Q2.5, Q97.5, prior_only_empty_cell)
write_csv(fixef_tbl, "results/dale/B2c_fixed_effects.csv")

varcomp <- map_dfr(names(fits), function(cv) {
  vc <- VarCorr(fits[[cv]])
  grp <- map_dfr(names(vc), function(g)
    tibble(group = g, sd = vc[[g]]$sd[, "Estimate"]))
  sig <- summary(fits[[cv]])$spec_pars["sigma", "Estimate"]
  bind_rows(grp, tibble(group = "sigma", sd = sig)) %>% mutate(model = cv)
}) %>% select(model, group, sd)
write_csv(varcomp, "results/dale/B2c_variance_components.csv")
cat("\nVariance components:\n"); print(varcomp)

r2 <- map_dfr(names(fits), function(cv) {
  r <- bayes_R2(fits[[cv]])
  tibble(model = cv, R2 = r[, "Estimate"], Q2.5 = r[, "Q2.5"], Q97.5 = r[, "Q97.5"])
})
write_csv(r2, "results/dale/B2c_r2_summary.csv")
cat("\nBayesian R2:\n"); print(r2)

cat("\nB2c diagnostics complete. Tables in results/.\n")
