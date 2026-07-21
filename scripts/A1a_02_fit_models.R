# ============================================================
# A1a_02_fit_models.R
# Fit A1a models: 4 responses x 1 variant (full-interaction) = 4 models.
# No phylogeny -> fast relative to A2c/A2d. Sampler settings match A1b
# so the A1a-vs-A1b LOO comparison differs only in the fixed effects.
# ============================================================

library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/A1a_model_setup.RData")

fit_one <- function(resp) {
  fam  <- if (response_families[[resp]] == "lognormal") lognormal() else gaussian()
  pri  <- if (response_families[[resp]] == "lognormal") priors_lognormal else priors_gaussian
  form <- formulas[[resp]]
  f    <- paste0("fits/A1a_fit_", resp, ".RData")

  if (file.exists(f)) { cat("skip (exists):", resp, "\n"); load(f); return(fit) }

  cat("\n=== fitting", resp, "===\n")
  # adapt_delta 0.9 matches A1b and the original A1 (both converged with
  # this fixed structure / this random structure respectively). If any
  # divergences appear with the combined model, raise to 0.95.
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
for (resp in responses) fits[[resp]] <- fit_one(resp)

save(fits, file = "fits/A1a_all_fits.RData")
cat("\nAll A1a fits complete. Saved to fits/A1a_all_fits.RData\n")
