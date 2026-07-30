# ============================================================
# A1b_03_diagnostics_summary.R
# Convergence, fixed effects, and Bayesian R2 for the A1b models.
# Outputs: results/A1b_convergence.csv, A1b_fixed_effects.csv,
#          results/A1b_r2_summary.csv
# ============================================================

library(brms)
library(posterior)
library(tidyverse)

load("fits/A1b_all_fits.RData")   # `fits` (named list)
dir.create("results/cooney", showWarnings = FALSE)

# Empty interaction cell after the Temperate Open x Plantation collapse:
# its coefficient is prior-only (no data) and is flagged, not interpreted.
EMPTY_CELL <- "Biome4TemperateOpen:Predominant_simplePlantationforest"

conv <- map_dfr(names(fits), function(tag) {
  fit <- fits[[tag]]
  s   <- summary(fit)
  rh  <- c(s$fixed[, "Rhat"], s$spec_pars[, "Rhat"],
           unlist(lapply(s$random, function(x) x[, "Rhat"])))
  eb  <- c(s$fixed[, "Bulk_ESS"], s$spec_pars[, "Bulk_ESS"],
           unlist(lapply(s$random, function(x) x[, "Bulk_ESS"])))
  np  <- nuts_params(fit)
  tibble(model = tag,
         max_rhat = round(max(rh, na.rm = TRUE), 4),
         min_ess_bulk = round(min(eb, na.rm = TRUE)),
         n_divergent = sum(np$Value[np$Parameter == "divergent__"]))
})
write_csv(conv, "results/cooney/A1b_convergence.csv")
cat("Convergence:\n"); print(conv)

fixef_tbl <- map_dfr(names(fits), function(tag) {
  fx <- as.data.frame(fixef(fits[[tag]]))
  fx$parameter <- rownames(fx)
  fx$model <- tag
  fx$prior_only_empty_cell <- grepl(EMPTY_CELL, fx$parameter, fixed = TRUE)
  as_tibble(fx)
}) %>%
  select(model, parameter, Estimate, Est.Error, Q2.5, Q97.5, prior_only_empty_cell)
write_csv(fixef_tbl, "results/cooney/A1b_fixed_effects.csv")

r2 <- map_dfr(names(fits), function(tag) {
  r <- bayes_R2(fits[[tag]])
  tibble(model = tag, R2 = r[, "Estimate"], Q2.5 = r[, "Q2.5"], Q97.5 = r[, "Q97.5"])
})
write_csv(r2, "results/cooney/A1b_r2_summary.csv")
cat("\nBayesian R2:\n"); print(r2)

cat("\nA1b diagnostics complete. Tables in results/.\n")
