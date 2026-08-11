# ============================================================
# B1a_01_setup.R
# Dale-colour counterpart of A1a.
#
# B1b covariate plus land-use x trophic niche and land-use x body
# mass interactions, so it is a nested superset of B1b and the
# B1a_03 LOO mirrors the A1a-vs-A1b comparison. Identical to A1a in
# everything but the colour source.
# ============================================================

source("scripts/00_config.R")
library(tidyverse)
library(brms)
library(cmdstanr)

dat <- load_community_data()

responses         <- RESP_DALE
response_families <- RESP_FAMILY[responses]

make_formula <- function(resp) {
  bf(as.formula(paste0(
    resp, " ~ Biome4 * Predominant_simple + Trophic.Niche + z_logMass",
    " + Predominant_simple:z_logMass",
    " + Predominant_simple:Trophic.Niche",
    " + (1|SS) + (1|SSB) + (1|SSBS)")))
}
formulas <- setNames(lapply(responses, make_formula), responses)

priors_lognormal <- model_priors("lognormal")
priors_gaussian  <- model_priors("gaussian")

dir.create("fits", showWarnings = FALSE)
save(dat, responses, response_families, formulas,
     priors_lognormal, priors_gaussian, nthreads,
     file = "fits/B1a_model_setup.RData")
cat("\nB1a setup complete. Saved to fits/B1a_model_setup.RData\n")
