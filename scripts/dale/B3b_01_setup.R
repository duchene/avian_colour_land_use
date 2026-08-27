# ============================================================
# B3b_01_setup.R
# Dale-colour counterpart of A3b.
#
# Identical to A3b in data, priors, sampler and random structure.
# The only difference is the colour source: the four Dale et al.
# (2015) plumage metrics in place of the Cooney UVS metrics.
# B3b is the primary community model of the Dale (B) family and the
# reduced model in the B3a LOO.
# ============================================================

source("scripts/00_config.R")
library(tidyverse)
library(brms)
library(cmdstanr)

dat <- load_community_data()

responses         <- RESP_DALE
response_families <- RESP_FAMILY[responses]

make_formulas <- function(resp) {
  list(
    base      = bf(as.formula(paste0(
      resp, " ~ Biome4 * Predominant_simple + (1|SS) + (1|SSB) + (1|SSBS)"))),
    covariate = bf(as.formula(paste0(
      resp, " ~ Biome4 * Predominant_simple + Trophic.Niche + z_logMass",
      " + (1|SS) + (1|SSB) + (1|SSBS)")))
  )
}
formulas <- setNames(lapply(responses, make_formulas), responses)

priors_lognormal <- model_priors("lognormal")
priors_gaussian  <- model_priors("gaussian")

dir.create("fits", showWarnings = FALSE)
save(dat, responses, response_families, formulas,
     priors_lognormal, priors_gaussian, nthreads,
     file = "fits/B3b_model_setup.RData")
cat("\nB3b setup complete. Saved to fits/B3b_model_setup.RData\n")
