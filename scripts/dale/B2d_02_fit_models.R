# ============================================================
# B2d_02_fit_models.R
# Fit B2d models (one per Dale colour predictor). Gaussian family
# (response = signed paired difference). Mirrors A2d_02.
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/B2d_model_setup.RData")   # dat_final, A, colour_vars, formulas, priors_B2d, nthreads

fit_one <- function(cv) {
  f <- paste0("fits/B2d_fit_", cv, ".RData")
  if (file.exists(f)) { cat("skip (exists):", cv, "\n"); load(f); return(fit) }

  d_cv <- dat_final %>% filter(!is.na(.data[[cv]]))
  d_cv$phylo <- droplevels(d_cv$phylo)
  tips <- levels(d_cv$phylo)
  A_cv <- A[tips, tips]

  cat("\n=== fitting B2d", cv, "| n =", nrow(d_cv),
      "| tips =", length(tips), "===\n")
  fit <- brm(formulas[[cv]], data = d_cv, data2 = list(A = A_cv),
             family = gaussian(), prior = priors_B2d,
             chains = 4, iter = 3500, warmup = 1500,
             cores = 4, threads = threading(nthreads),
             control = list(adapt_delta = 0.95, max_treedepth = 12),
             seed = 1, refresh = 200)
  save(fit, file = f)
  cat("saved", f, "\n")
  fit
}

fits <- setNames(lapply(colour_vars, fit_one), colour_vars)
save(fits, file = "fits/B2d_all_fits.RData")
cat("\nAll B2d fits complete. Saved to fits/B2d_all_fits.RData\n")
