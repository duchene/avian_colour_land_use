# ============================================================
# A3b_02_fit_models.R
# Fit A3b models: 4 responses x {base, covariate} = 8 models.
# No phylogeny -> fast relative to A2c/A2d.
# ============================================================

# Attach rstan before any saved brmsfit is touched: loaded implicitly (by
# methods dispatch on the stanfit slot) it fails in the GenomeDK container.
suppressMessages(library(rstan))
source("scripts/00_config.R")   # family_prior(), RESP_MAIN
library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")

load("fits/A3b_model_setup.RData")
# Threads per chain from THIS job's allocation (the setup .RData stores the
# value of whatever job ran the setup).
if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK")))
  nthreads <- max(1, floor(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) / 4))

fit_one <- function(resp, variant) {
  fp   <- family_prior(resp, "record")   # 2026-11: one source for family + priors
  fam  <- fp$family
  pri  <- fp$prior
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
