# ============================================================
# B3a_02_fit_models.R
# Fit B3a models: 4 Dale responses x 1 variant (full-interaction) = 4
# models. No phylogeny. Sampler settings match A3a / B3b so the
# B3a-vs-B3b LOO comparison differs only in the fixed effects.
# ============================================================

library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/B3a_model_setup.RData")

fit_one <- function(resp) {
  fam  <- if (response_families[[resp]] == "lognormal") lognormal() else gaussian()
  pri  <- if (response_families[[resp]] == "lognormal") priors_lognormal else priors_gaussian
  form <- formulas[[resp]]
  f    <- paste0("fits/B3a_fit_", resp, ".RData")

  if (file.exists(f)) { cat("skip (exists):", resp, "\n"); load(f); return(fit) }

  cat("\n=== fitting", resp, "===\n")
  # adapt_delta 0.9 matches A3a. Raise to 0.95 if any divergences appear.
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

save(fits, file = "fits/B3a_all_fits.RData")
cat("\nAll B3a fits complete. Saved to fits/B3a_all_fits.RData\n")
