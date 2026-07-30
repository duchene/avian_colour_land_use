# ============================================================
# A2b_02_fit_models.R
# Analysis 2b: Fit paired-difference models
# Interaction:  diff_abund ~ colour * land-use + (1|SSBS) + (1|Reference)
# Main-effects: diff_abund ~ colour + land-use + (1|SSBS) + (1|Reference)
# Gaussian family
# ============================================================

library(dplyr)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
ncores <- parallel::detectCores()
options(mc.cores = ncores)

# Load setup
load("fits/A2b_model_setup.RData")

# Threading: 4 chains, use remaining cores for within-chain parallelism
nthreads_per_chain <- max(2, floor(ncores / 4))
cat("Detected cores:", ncores, "\n")
cat("Threads per chain:", nthreads_per_chain, "\n")
cat("Total rows:", nrow(dat_final), "\n")

dir.create("fits", showWarnings = FALSE)

# ============================================================
# HELPER: fit or load from checkpoint
# ============================================================
fit_or_load <- function(checkpoint_file, formula, data, prior) {
  if (file.exists(checkpoint_file)) {
    cat("  Loading existing fit:", checkpoint_file, "\n")
    load(checkpoint_file)
    return(fit)
  }

  cat("  Fitting | N =", nrow(data), "\n")
  t0 <- proc.time()

  fit <- brm(
    formula  = formula,
    data     = data,
    family   = gaussian(),
    prior    = prior,
    chains   = 4,
    cores    = 4,
    iter     = 2000,
    warmup   = 1000,
    threads  = threading(nthreads_per_chain, grainsize = 1),
    control  = list(adapt_delta = 0.90, max_treedepth = 10),
    seed     = 42,
    refresh  = 500,
    silent   = 1
  )

  elapsed <- (proc.time() - t0)["elapsed"]
  cat("  Elapsed:", round(elapsed / 60, 1), "min\n")

  save(fit, file = checkpoint_file)
  cat("  Saved:", checkpoint_file, "\n")

  print(summary(fit))
  return(fit)
}

# ============================================================
# FIT INTERACTION MODELS (colour * land-use)
# ============================================================

cat("\n=== INTERACTION MODELS ===\n")
fits_A2b <- list()

for (cv in colour_vars) {
  cat("\n---", cv, "(interaction) ---\n")
  raw_colour <- sub("^z_", "", cv)
  d <- dat_final %>%
    filter(!is.na(diff_abund) & !is.na(.data[[raw_colour]])) %>%
    droplevels()

  fits_A2b[[cv]] <- fit_or_load(
    checkpoint_file = paste0("fits/A2b_fit_", cv, ".RData"),
    formula = formulas[[cv]],
    data = d,
    prior = priors_A2b
  )
}

save(fits_A2b, file = "fits/A2b_all_fits.RData")
cat("\nInteraction models saved to fits/A2b_all_fits.RData\n")

# ============================================================
# FIT MAIN-EFFECTS MODELS (colour + land-use, no interaction)
# ============================================================

cat("\n=== MAIN-EFFECTS MODELS ===\n")

formulas_main <- setNames(
  lapply(colour_vars, function(cv) {
    bf(as.formula(paste0(
      "diff_abund ~ ", cv, " + Predominant_simple + (1 | SSBS) + (1 | Reference)"
    )))
  }),
  colour_vars
)

fits_A2b_main <- list()

for (cv in colour_vars) {
  cat("\n---", cv, "(main effects) ---\n")
  raw_colour <- sub("^z_", "", cv)
  d <- dat_final %>%
    filter(!is.na(diff_abund) & !is.na(.data[[raw_colour]])) %>%
    droplevels()

  fits_A2b_main[[cv]] <- fit_or_load(
    checkpoint_file = paste0("fits/A2b_fit_main_", cv, ".RData"),
    formula = formulas_main[[cv]],
    data = d,
    prior = priors_A2b
  )
}

save(fits_A2b_main, file = "fits/A2b_all_fits_main.RData")
cat("\nMain-effects models saved to fits/A2b_all_fits_main.RData\n")

cat("\nAll A2b fitting complete.\n")
