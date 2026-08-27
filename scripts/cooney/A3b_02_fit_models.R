# ============================================================
# A3b_02_fit_models.R
# Fit A3b models: 4 responses x {base, covariate} = 8 models.
# No phylogeny -> fast relative to A2c/A2d.
# ============================================================

library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/A3b_model_setup.RData")

fit_one <- function(resp, variant) {
  fam  <- if (response_families[[resp]] == "lognormal") lognormal() else gaussian()
  pri  <- if (response_families[[resp]] == "lognormal") priors_lognormal else priors_gaussian
  form <- formulas[[resp]][[variant]]
  tag  <- paste0(resp, "_", variant)
  f    <- paste0("fits/A3b_fit_", tag, ".RData")

  if (file.exists(f)) { cat("skip (exists):", tag, "\n"); load(f); return(fit) }

  cat("\n=== fitting", tag, "===\n")
  fit <- brm(form, data = dat, family = fam, prior = pri,
             chains = 4, iter = 3000, warmup = 1500,
             cores = 4, threads = threading(nthreads),
             control = list(adapt_delta = 0.9, max_treedepth = 10),
             seed = 1, refresh = 200)
  save(fit, file = f)
  cat("saved", f, "\n")
  fit
}

fits <- list()
for (resp in responses) {
  for (variant in c("base", "covariate")) {
    fits[[paste0(resp, "_", variant)]] <- fit_one(resp, variant)
  }
}
save(fits, file = "fits/A3b_all_fits.RData")
cat("\nAll A3b fits complete. Saved to fits/A3b_all_fits.RData\n")
