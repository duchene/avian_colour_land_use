# ============================================================
# A1_02_fit_models.R
# Analysis 1: Fit phylogenetic regression models
# Colour ~ land-use proportion scores
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

# Load setup
load("fits/A1_model_setup.RData")
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

fits_A1 <- list()

for (resp in responses) {
  checkpoint_file <- paste0("fits/A1_fit_", resp, ".RData")

  if (file.exists(checkpoint_file)) {
    cat("\nLoading existing fit for:", resp, "\n")
    load(checkpoint_file)
    fits_A1[[resp]] <- fit
    next
  }

  d <- spdat %>% filter(!is.na(.data[[resp]]))
  f <- make_formula(resp)

  # Select family and priors per response
  fam_name <- response_families[[resp]]
  fam <- if (fam_name == "gaussian") gaussian() else lognormal()
  pr  <- if (fam_name == "gaussian") priors_A1_gaussian else priors_A1

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

  fits_A1[[resp]] <- fit
  save(fit, file = checkpoint_file)
  cat(">>> Saved:", checkpoint_file, "<<<\n")
  print(summary(fit))
}

# ============================================================
# SAVE ALL FITS
# ============================================================

save(fits_A1, file = "fits/A1_all_fits.RData")
cat("\nAll A1 models saved to fits/A1_all_fits.RData\n")
