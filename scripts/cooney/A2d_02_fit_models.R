# ============================================================
# A2d_02_fit_models.R
# Fit A2d models (one per colour predictor). Gaussian family
# (response = signed paired difference). Per model, subset to
# non-NA colour rows and align the phylo covariance matrix A.
# ============================================================

# Attach rstan before any saved brmsfit is touched: loaded implicitly (by
# methods dispatch on the stanfit slot) it fails in the GenomeDK container.
suppressMessages(library(rstan))
library(tidyverse)
library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/A2d_model_setup.RData")   # dat_final, A, colour_vars, formulas, priors_A2d, nthreads
# Threads per chain from THIS job's allocation (the setup .RData stores the
# value of whatever job ran the setup).
if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK")))
  nthreads <- max(1, floor(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) / 4))

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

# ONLY=z_dichroplate (comma-separated) fits a subset, so the colour
# metrics can run as parallel Slurm jobs; rerun without ONLY to collect.
only <- Sys.getenv("ONLY")
if (nzchar(only)) {
  todo <- intersect(colour_vars, strsplit(only, ",")[[1]])
  for (cv in todo) fit_one(cv)
  quit(save = "no")
}
fits <- setNames(lapply(colour_vars, fit_one), colour_vars)
save(fits, file = "fits/A2d_all_fits.RData")
cat("\nAll A2d fits complete. Saved to fits/A2d_all_fits.RData\n")
