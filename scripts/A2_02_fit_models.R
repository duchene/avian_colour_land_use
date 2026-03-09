# ============================================================
# A2_02_fit_models.R
# Analysis 2: Fit abundance ~ colour x land-use models
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

# Load setup
load("fits/A2_model_setup.RData")
cat("Setup loaded. nthreads =", nthreads, "\n")

dir.create("fits", showWarnings = FALSE)

# ============================================================
# FIT MODELS
# ============================================================

fits_A2 <- list()

for (cv in colour_vars) {
  checkpoint_file <- paste0("fits/A2_fit_relabund_", cv, ".RData")

  if (file.exists(checkpoint_file)) {
    cat("\nLoading existing fit for:", cv, "\n")
    load(checkpoint_file)
    fits_A2[[cv]] <- fit
    next
  }

  # Prepare data (drop NAs for this colour variable)
  raw_colour <- sub("^z_", "", cv)
  d <- dat %>%
    filter(!is.na(abundance) & !is.na(.data[[raw_colour]])) %>%
    droplevels()

  cat("\nFitting:", cv, "| Family: lognormal | N =", nrow(d), "\n")

  fit <- brm(
    formula  = formulas[[cv]],
    data     = d,
    family   = lognormal(),
    prior    = priors_A2,
    chains   = 4,
    cores    = 4,
    iter     = 3000,
    warmup   = 1500,
    threads  = threading(nthreads),
    control  = list(adapt_delta = 0.90, max_treedepth = 10),
    seed     = 42,
    refresh  = 200
  )

  fits_A2[[cv]] <- fit
  save(fit, file = checkpoint_file)
  cat(">>> Saved:", checkpoint_file, "<<<\n")

  print(summary(fit))
}

# ============================================================
# SAVE ALL FITS
# ============================================================

save(fits_A2, file = "fits/A2_all_fits_relabund.RData")
cat("\nAll A2 models saved to fits/A2_all_fits_relabund.RData\n")
