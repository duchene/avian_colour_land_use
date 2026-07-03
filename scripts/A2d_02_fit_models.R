# ============================================================
# A2d_02_fit_models.R
# Fit A2d models (one per colour predictor). Gaussian family
# (response = signed paired difference). Per model, subset to
# non-NA colour rows and align the phylo covariance matrix A.
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/A2d_model_setup.RData")   # dat_final, A, colour_vars, formulas, priors_A2d, nthreads

fit_one <- function(cv) {
  f <- paste0("fits/A2d_fit_", cv, ".RData")
  if (file.exists(f)) { cat("skip (exists):", cv, "\n"); load(f); return(fit) }

  d_cv <- dat_final %>% filter(!is.na(.data[[cv]]))
  d_cv$phylo <- droplevels(d_cv$phylo)
  tips <- levels(d_cv$phylo)
  A_cv <- A[tips, tips]

  cat("\n=== fitting A2d", cv, "| n =", nrow(d_cv),
      "| tips =", length(tips), "===\n")
  fit <- brm(formulas[[cv]], data = d_cv, data2 = list(A = A_cv),
             family = gaussian(), prior = priors_A2d,
             chains = 4, iter = 3500, warmup = 1500,
             cores = 4, threads = threading(nthreads),
             control = list(adapt_delta = 0.95, max_treedepth = 12),
             seed = 1, refresh = 200)
  save(fit, file = f)
  cat("saved", f, "\n")
  fit
}

fits <- setNames(lapply(colour_vars, fit_one), colour_vars)
save(fits, file = "fits/A2d_all_fits.RData")
cat("\nAll A2d fits complete. Saved to fits/A2d_all_fits.RData\n")
