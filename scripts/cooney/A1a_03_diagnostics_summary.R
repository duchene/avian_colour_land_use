# ============================================================
# A1a_03_diagnostics_summary.R
# Convergence, fixed effects, Bayesian R2 for the A1a models, and the
# LOO comparison against the A1b covariate model (the nested reduced
# model). Output set mirrors A1b exactly:
#   results/A1a_convergence.csv
#   results/A1a_fixed_effects.csv   (biome/trophic/mass interactions
#                                    live here; the figure filters them)
#   results/A1a_r2_summary.csv
#   results/A1a_vs_A1b_loo.csv      (does the guild/size interaction set
#                                    add predictive value beyond A1b?)
# Requires the A1b covariate fits (run A1b first) for the LOO step.
# ============================================================

library(brms)
library(posterior)
source("scripts/00_config.R")
library(tidyverse)
library(loo)

load("fits/A1a_all_fits.RData")   # `fits` (named by response, single variant)
responses <- names(fits)
dir.create("results/cooney", showWarnings = FALSE)

# Prior-only biome x land-use cells, derived from the fitted data.
load("fits/A1a_model_setup.RData")   # `dat`
EMPTY_CELLS <- empty_cell_params(dat)
cat("Prior-only empty cells:",
    if (length(EMPTY_CELLS)) paste(EMPTY_CELLS, collapse = ", ") else "none", "\n")
SPARSE_TROPHIC <- "Trophic.NicheHerbivoreterrestrial"

# ---- convergence ------------------------------------------------------------
conv <- map_dfr(responses, function(tag) {
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
write_csv(conv, "results/cooney/A1a_convergence.csv")
cat("Convergence:\n"); print(conv)

# ---- fixed effects ----------------------------------------------------------
fixef_tbl <- map_dfr(responses, function(tag) {
  fx <- as.data.frame(fixef(fits[[tag]]))
  fx$parameter <- rownames(fx)
  fx$model <- tag
  fx$prior_only_empty_cell <- fx$parameter %in% EMPTY_CELLS
  fx$sparse_trophic        <- grepl(SPARSE_TROPHIC, fx$parameter, fixed = TRUE)
  as_tibble(fx)
}) %>%
  select(model, parameter, Estimate, Est.Error, Q2.5, Q97.5,
         prior_only_empty_cell, sparse_trophic)
write_csv(fixef_tbl, "results/cooney/A1a_fixed_effects.csv")
cat("\nSaved results/A1a_fixed_effects.csv (",
    nrow(fixef_tbl), "rows )\n")

# ---- Bayesian R2 ------------------------------------------------------------
r2 <- map_dfr(responses, function(tag) {
  r <- bayes_R2(fits[[tag]])
  tibble(model = tag, R2 = r[, "Estimate"], Q2.5 = r[, "Q2.5"], Q97.5 = r[, "Q97.5"])
})
write_csv(r2, "results/cooney/A1a_r2_summary.csv")
cat("\nBayesian R2:\n"); print(r2)

# ---- LOO: A1a (full) vs A1b covariate (nested reduced) ----------------------
# Positive elpd_diff_A1a_minus_A1b => the trophic/mass interaction set
# improves out-of-sample prediction over A1b.
loo_tbl <- map_dfr(responses, function(resp) {
  a1b_f <- paste0("fits/A1b_fit_", resp, "_covariate.RData")
  if (!file.exists(a1b_f)) {
    message("  A1b covariate fit missing for ", resp,
            " (run A1b first) — skipping LOO"); return(NULL)
  }
  e <- new.env(); load(a1b_f, envir = e)   # loads `fit` (A1b covariate)
  la <- loo(fits[[resp]]); lb <- loo(e$fit)
  cmp <- loo_compare(la, lb)               # for the standard error of the diff
  ea  <- la$estimates["elpd_loo", "Estimate"]
  eb  <- lb$estimates["elpd_loo", "Estimate"]
  tibble(response = resp,
         elpd_A1a = ea, elpd_A1b = eb,
         elpd_diff_A1a_minus_A1b = ea - eb,
         se_diff = cmp[2, "se_diff"],
         preferred = if (ea >= eb) "A1a (full interactions)"
                     else "A1b (biome x land-use only)")
})
if (nrow(loo_tbl)) {
  write_csv(loo_tbl, "results/cooney/A1a_vs_A1b_loo.csv")
  cat("\nLOO (A1a full vs A1b covariate):\n"); print(loo_tbl)
}

cat("\nA1a diagnostics complete. Tables in results/.\n")
