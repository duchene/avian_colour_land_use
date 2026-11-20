# ============================================================
# A2d_colours_02_fit.R
# A2d (abundance change against primary vegetation) refitted once per
# Delhey et al. 2023 plumage colour, each colour as its own predictor.
#
# Identical to the reported A2d (A2d_01_setup.R / A2d_02_fit_models.R):
# same records, formula, priors, phylogeny and sampler settings. Only the
# colour predictor changes: the proportion of male plumage in one colour,
# standardised (z), in place of male colourfulness or dichromatism. No
# ratios or composites are formed (co-author decision, 2026-11-13).
# Records of species without Delhey data are dropped from these fits.
#
# ONLY=<colour> fits one colour (one Slurm array task each, see
# hpc/14_a2d_colours.sh); without ONLY every missing fit is run.
# Fits: fits/A2d_fit_z_col_<colour>.RData
# ============================================================

# Attach rstan before any saved brmsfit is touched (GenomeDK container).
suppressMessages(library(rstan))
suppressPackageStartupMessages({ library(tidyverse); library(brms); library(cmdstanr) })
options(brms.backend = "cmdstanr")
source("scripts/cooney/delhey_colours.R")

load("fits/A2d_model_setup.RData")   # dat_final, A, priors_A2d, nthreads, ...
if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK")))
  nthreads <- max(1, floor(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) / 4))

dat_final <- add_delhey_colours(dat_final, dat_final$phylo)
COL_VARS  <- paste0("z_col_", DELHEY_COLOURS)
for (cc in DELHEY_COLOURS) dat_final[[paste0("z_col_", cc)]] <- (dat_final[[cc]] - mean(dat_final[[cc]], na.rm = TRUE)) / sd(dat_final[[cc]], na.rm = TRUE)

# Sampler length as for the reported A2d; A2D_ITER / A2D_WARMUP lengthen
# it (black was refitted with iter 7500, warmup 2500 after a max Rhat of
# 1.02 at the default length). FORCE=1 refits over an existing fit.
ITER   <- as.integer(Sys.getenv("A2D_ITER", "3500"))
WARMUP <- as.integer(Sys.getenv("A2D_WARMUP", "1500"))

fit_one <- function(cv) {
  f <- paste0("fits/A2d_fit_", cv, ".RData")
  if (file.exists(f) && !nzchar(Sys.getenv("FORCE"))) { cat("skip (exists):", cv, "\n"); return(invisible(NULL)) }
  d_cv <- dat_final %>% filter(!is.na(.data[[cv]]))
  d_cv$phylo <- droplevels(d_cv$phylo)
  tips <- levels(d_cv$phylo)
  form <- bf(as.formula(paste0(
    "diff_abund ~ ", cv, " * Biome4 + ", cv, " * Predominant_simple",
    " + Biome4 * Predominant_simple",
    " + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))")))
  cat("\n=== fitting A2d", cv, "| n =", nrow(d_cv), "| tips =", length(tips), "===\n")
  fit <- brm(form, data = d_cv, data2 = list(A = A[tips, tips]),
             family = gaussian(), prior = priors_A2d,
             chains = 4, iter = ITER, warmup = WARMUP,
             cores = 4, threads = threading(nthreads),
             control = list(adapt_delta = 0.95, max_treedepth = 12),
             seed = 1, refresh = 200)
  save(fit, file = f)
  cat("saved", f, "\n")
}

only <- Sys.getenv("ONLY")
todo <- if (nzchar(only)) intersect(COL_VARS, paste0("z_col_", strsplit(only, ",")[[1]])) else COL_VARS
for (cv in todo) fit_one(cv)
