# ============================================================
# A1_02_fit_models.R
# Analysis 1: Fit all brms models with checkpoint saves
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

# Load setup
load("fits/A1_model_setup.RData")
cat("Setup loaded. nthreads =", nthreads, "\n")

# ============================================================
# FITTING FUNCTION
# ============================================================

fit_brms_model <- function(data, response, formula, family_name,
                           priors, nthreads) {

  d <- data %>%
    filter(!is.na(.data[[response]])) %>%
    mutate(y = .data[[response]])

  fam <- switch(family_name,
    "lognormal" = lognormal(),
    "student"   = student(),
    "gaussian"  = gaussian(),
    gaussian()
  )

  # Student-t models need more careful sampling
  if (family_name == "student") {
    iter <- 4000
    warmup <- 2000
    adapt_delta <- 0.95
    max_treedepth <- 12
  } else {
    iter <- 3000
    warmup <- 1500
    adapt_delta <- 0.90
    max_treedepth <- 10
  }

  cat("\nFitting:", response, "| Family:", family_name, "| N =", nrow(d), "\n")
  cat("Settings: iter =", iter, ", adapt_delta =", adapt_delta, "\n")

  brm(
    formula = formula,
    data    = d,
    family  = fam,
    prior   = priors,
    chains  = 4,
    cores   = 4,
    iter    = iter,
    warmup  = warmup,
    threads = threading(nthreads),
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    seed = 42,
    refresh = 200
  )
}

# ============================================================
# FIT REDUCED MODELS (with checkpoints)
# ============================================================

fits_reduced <- list()

for (i in seq_len(nrow(response_specs))) {
  resp <- response_specs$response[i]
  fam  <- response_specs$family[i]

  checkpoint_file <- paste0("fits/A1_fit_", resp, "_reduced.RData")
  if (file.exists(checkpoint_file)) {
    cat("\nLoading existing fit for:", resp, "\n")
    load(checkpoint_file)
    fits_reduced[[resp]] <- fit
    next
  }

  priors <- if (fam == "student") student_priors else base_priors

  fit <- fit_brms_model(
    data = dat, response = resp,
    formula = formula_reduced, family_name = fam,
    priors = priors, nthreads = nthreads
  )

  fits_reduced[[resp]] <- fit

  save(fit, file = checkpoint_file)
  cat(">>> Saved checkpoint:", checkpoint_file, "<<<\n")

  print(summary(fit))
}

# ============================================================
# FIT MINIMAL MODELS (for comparison)
# ============================================================

fits_minimal <- list()

for (i in seq_len(nrow(response_specs))) {
  resp <- response_specs$response[i]
  fam  <- response_specs$family[i]

  checkpoint_file <- paste0("fits/A1_fit_", resp, "_minimal.RData")
  if (file.exists(checkpoint_file)) {
    cat("\nLoading existing fit for:", resp, "minimal\n")
    load(checkpoint_file)
    fits_minimal[[resp]] <- fit
    next
  }

  priors <- if (fam == "student") student_priors else base_priors

  fit <- fit_brms_model(
    data = dat, response = resp,
    formula = formula_minimal, family_name = fam,
    priors = priors, nthreads = nthreads
  )

  fits_minimal[[resp]] <- fit

  save(fit, file = checkpoint_file)
  cat(">>> Saved checkpoint:", checkpoint_file, "<<<\n")

  print(summary(fit))
}

# ============================================================
# SAVE ALL FITS
# ============================================================

save(fits_reduced, fits_minimal, file = "fits/A1_all_fits.RData")
cat("\n\nAll models saved to fits/A1_all_fits.RData\n")
