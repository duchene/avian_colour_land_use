# ============================================================
# A3a_02_fit_models.R
# Fit A3a models: 4 responses x 1 variant (full-interaction) = 4 models.
# No phylogeny -> fast relative to A2c/A2d. Sampler settings match A3b
# so the A3a-vs-A3b LOO comparison differs only in the fixed effects.
# ============================================================

# Attach rstan before any saved brmsfit is touched: loaded implicitly (by
# methods dispatch on the stanfit slot) it fails in the GenomeDK container.
suppressMessages(library(rstan))
source("scripts/00_config.R")   # family_prior(), RESP_MAIN
library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/A3a_model_setup.RData")
# Threads per chain from THIS job's allocation (the setup .RData stores the
# value of whatever job ran the setup).
if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK")))
  nthreads <- max(1, floor(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) / 4))

fit_one <- function(resp) {
  fp   <- family_prior(resp, "record")   # 2026-11: one source for family + priors
  fam  <- fp$family
  pri  <- fp$prior
  form <- formulas[[resp]]
  f    <- paste0("fits/A3a_fit_", resp, ".RData")

  if (file.exists(f)) { cat("skip (exists):", resp, "\n"); load(f); return(fit) }

  cat("\n=== fitting", resp, "===\n")
  # adapt_delta 0.9 matches A3b and the original A3 (both converged with
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

save(fits, file = "fits/A3a_all_fits.RData")
cat("\nAll A3a fits complete. Saved to fits/A3a_all_fits.RData\n")
