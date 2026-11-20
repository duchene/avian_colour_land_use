# ============================================================
# A1b_02_fit_models.R
# Analysis 1b: Fit phylogenetic regression models
# Colour ~ land-use proportion scores
# ============================================================

# Attach rstan before any saved brmsfit is touched: loaded implicitly (by
# methods dispatch on the stanfit slot) it fails in the GenomeDK container.
suppressMessages(library(rstan))
source("scripts/00_config.R")   # family_prior(), RESP_MAIN
library(tidyverse)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

# Load setup
load("fits/A1b_model_setup.RData")
# Threads per chain from THIS job's allocation (the setup .RData stores the
# value of whatever job ran the setup).
if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK")))
  nthreads <- max(1, floor(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) / 4))
cat("Setup loaded. nthreads =", nthreads, "\n")
cat("Species:", nrow(spdat), "\n")
cat("Predictors:", paste(lu_predictors, collapse = ", "), "\n")
cat("Phylogenetic matrix:", ifelse(is.null(A), "NOT AVAILABLE", "loaded"), "\n")

dir.create("fits", showWarnings = FALSE)

# ============================================================
# MODEL FORMULA
# ============================================================
# Colour ~ 0 + land-use proportions + (1|gr(phylo, cov = A))
# No intercept because the 5 proportions sum to 1 per species,
# so they span the intercept. Each coefficient represents the
# expected colour for a species found exclusively in that land use.
# The phylogenetic random effect accounts for shared ancestry.

predictor_string <- paste(lu_predictors, collapse = " + ")

make_formula <- function(response) {
  if (!is.null(A)) {
    bf(as.formula(paste0(response, " ~ 0 + ", predictor_string,
                         " + (1 | gr(phylo, cov = A))")))
  } else {
    bf(as.formula(paste0(response, " ~ 0 + ", predictor_string)))
  }
}

# ============================================================
# FIT MODELS
# ============================================================

fits_A1b <- list()

for (resp in responses) {
  checkpoint_file <- paste0("fits/A1b_fit_", resp, ".RData")

  if (file.exists(checkpoint_file)) {
    cat("\nLoading existing fit for:", resp, "\n")
    load(checkpoint_file)
    fits_A1b[[resp]] <- fit
    next
  }

  d <- spdat %>% filter(!is.na(.data[[resp]]))
  f <- make_formula(resp)

  # Select family and priors per response
  fam_name <- response_families[[resp]]
  fp  <- family_prior(resp, "species")   # 2026-11: one source for family + priors
  fam <- fp$family
  pr  <- fp$prior

  cat("\nFitting:", resp, "| Family:", fam_name, "| N =", nrow(d), "\n")

  if (!is.null(A)) {
    fit <- brm(
      formula  = f,
      data     = d,
      data2    = list(A = A),
      family   = fam,
      prior    = pr,
      chains   = 4,
      cores    = 4,
      iter     = 4000,
      warmup   = 2000,
      threads  = threading(nthreads),
      control  = list(adapt_delta = 0.95, max_treedepth = 12),
      seed     = 42,
      refresh  = 200
    )
  } else {
    fit <- brm(
      formula  = f,
      data     = d,
      family   = fam,
      prior    = pr,
      chains   = 4,
      cores    = 4,
      iter     = 3000,
      warmup   = 1500,
      threads  = threading(nthreads),
      control  = list(adapt_delta = 0.90, max_treedepth = 10),
      seed     = 42,
      refresh  = 200
    )
  }

  fits_A1b[[resp]] <- fit
  save(fit, file = checkpoint_file)
  cat(">>> Saved:", checkpoint_file, "<<<\n")
  print(summary(fit))
}

# ============================================================
# SAVE ALL FITS
# ============================================================

save(fits_A1b, file = "fits/A1b_all_fits.RData")
cat("\nAll A1 models saved to fits/A1b_all_fits.RData\n")
