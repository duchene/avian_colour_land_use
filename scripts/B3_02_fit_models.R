# ============================================================
# B3_02_fit_models.R
# Analysis B3: fit phylogenetic regression models.
# Dale colour ~ land-use proportion scores. Mirrors A3_02.
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

load("fits/B3_model_setup.RData")
cat("Setup loaded. nthreads =", nthreads, "\n")
cat("Species:", nrow(spdat), " | predictors:", paste(lu_predictors, collapse = ", "), "\n")

dir.create("fits", showWarnings = FALSE)

predictor_string <- paste(lu_predictors, collapse = " + ")
make_formula <- function(response)
  bf(as.formula(paste0(response, " ~ 0 + ", predictor_string,
                       " + (1 | gr(phylo, cov = A))")))

fits_B3 <- list()
for (resp in responses) {
  checkpoint_file <- paste0("fits/B3_fit_", resp, ".RData")
  if (file.exists(checkpoint_file)) {
    cat("\nLoading existing fit for:", resp, "\n")
    load(checkpoint_file); fits_B3[[resp]] <- fit; next
  }

  d <- spdat %>% filter(!is.na(.data[[resp]]))
  f <- make_formula(resp)
  fam_name <- response_families[[resp]]
  fam <- if (fam_name == "gaussian") gaussian() else lognormal()
  pr  <- if (fam_name == "gaussian") priors_B3_gaussian else priors_B3

  cat("\nFitting:", resp, "| Family:", fam_name, "| N =", nrow(d), "\n")
  fit <- brm(
    formula = f, data = d, data2 = list(A = A),
    family = fam, prior = pr,
    chains = 4, cores = 4, iter = 4000, warmup = 2000,
    threads = threading(nthreads),
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    seed = 42, refresh = 200
  )
  fits_B3[[resp]] <- fit
  save(fit, file = checkpoint_file)
  cat(">>> Saved:", checkpoint_file, "<<<\n")
  print(summary(fit))
}

save(fits_B3, file = "fits/B3_all_fits.RData")
cat("\nAll B3 models saved to fits/B3_all_fits.RData\n")
