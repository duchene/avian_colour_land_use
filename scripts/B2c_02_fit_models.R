# ============================================================
# B2c_02_fit_models.R
# Fit B2c models (one per Dale colour predictor). Phylogenetic +
# full PREDICTS random structure. Mirrors A2c_02.
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/B2c_model_setup.RData")   # dat, A, colour_vars, formulas, priors_B2c, nthreads

fit_one <- function(cv) {
  f <- paste0("fits/B2c_fit_", cv, ".RData")
  if (file.exists(f)) { cat("skip (exists):", cv, "\n"); load(f); return(fit) }

  d_cv <- dat %>% filter(!is.na(.data[[cv]]))
  d_cv$phylo <- droplevels(d_cv$phylo)
  tips <- levels(d_cv$phylo)
  A_cv <- A[tips, tips]                          # align cov matrix to rows used

  cat("\n=== fitting B2c", cv, "| n =", nrow(d_cv),
      "| tips =", length(tips), "===\n")
  fit <- brm(formulas[[cv]], data = d_cv, data2 = list(A = A_cv),
             family = lognormal(), prior = priors_B2c,
             chains = 4, iter = 2000, warmup = 1000,
             cores = 4, threads = threading(nthreads),
             control = list(adapt_delta = 0.9, max_treedepth = 10),
             seed = 1, refresh = 200)
  save(fit, file = f)
  cat("saved", f, "\n")
  fit
}

fits <- setNames(lapply(colour_vars, fit_one), colour_vars)
save(fits, file = "fits/B2c_all_fits.RData")
cat("\nAll B2c fits complete. Saved to fits/B2c_all_fits.RData\n")
