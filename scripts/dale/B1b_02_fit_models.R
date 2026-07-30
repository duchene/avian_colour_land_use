# ============================================================
# B1b_02_fit_models.R
# Fit B1b models: 4 Dale responses x {base, covariate} = 8 models.
# No phylogeny. Mirrors A1b_02.
# ============================================================

library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/B1b_model_setup.RData")

fit_one <- function(resp, variant) {
  fam  <- if (response_families[[resp]] == "lognormal") lognormal() else gaussian()
  pri  <- if (response_families[[resp]] == "lognormal") priors_lognormal else priors_gaussian
  form <- formulas[[resp]][[variant]]
  tag  <- paste0(resp, "_", variant)
  f    <- paste0("fits/B1b_fit_", tag, ".RData")

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
save(fits, file = "fits/B1b_all_fits.RData")
cat("\nAll B1b fits complete. Saved to fits/B1b_all_fits.RData\n")
