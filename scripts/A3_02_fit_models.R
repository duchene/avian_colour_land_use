# ============================================================
# A3_02_fit_models.R
# Analysis 3: Fit phylogenetic regression models
# Colour ~ land-use proportion scores
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

# Load setup
load("fits/A3_model_setup.RData")
cat("Setup loaded. nthreads =", nthreads, "\n")
cat("Species:", nrow(spdat), "\n")
cat("Predictors:", paste(lu_predictors, collapse = ", "), "\n")
cat("Phylogenetic matrix:", ifelse(is.null(A), "NOT AVAILABLE", "loaded"), "\n")

dir.create("fits", showWarnings = FALSE)

# ============================================================
# MODEL FORMULA
# ============================================================
# Colour ~ land-use proportions + (1|gr(phylo, cov = A))
# The phylogenetic random effect accounts for shared ancestry.

predictor_string <- paste(lu_predictors, collapse = " + ")

make_formula <- function(response) {
  if (!is.null(A)) {
    bf(as.formula(paste0(response, " ~ ", predictor_string,
                         " + (1 | gr(phylo, cov = A))")))
  } else {
    bf(as.formula(paste0(response, " ~ ", predictor_string)))
  }
}

# ============================================================
# FIT MODELS
# ============================================================

fits_A3 <- list()

for (resp in responses) {
  checkpoint_file <- paste0("fits/A3_fit_", resp, ".RData")

  if (file.exists(checkpoint_file)) {
    cat("\nLoading existing fit for:", resp, "\n")
    load(checkpoint_file)
    fits_A3[[resp]] <- fit
    next
  }

  d <- spdat %>% filter(!is.na(.data[[resp]]))
  f <- make_formula(resp)

  cat("\nFitting:", resp, "| N =", nrow(d), "\n")

  if (!is.null(A)) {
    fit <- brm(
      formula  = f,
      data     = d,
      data2    = list(A = A),
      family   = lognormal(),
      prior    = priors_A3,
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
      family   = lognormal(),
      prior    = priors_A3,
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

  fits_A3[[resp]] <- fit
  save(fit, file = checkpoint_file)
  cat(">>> Saved:", checkpoint_file, "<<<\n")
  print(summary(fit))
}

# ============================================================
# SAVE ALL FITS
# ============================================================

save(fits_A3, file = "fits/A3_all_fits.RData")
cat("\nAll A3 models saved to fits/A3_all_fits.RData\n")
