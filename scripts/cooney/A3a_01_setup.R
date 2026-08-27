# ============================================================
# A3a_01_setup.R
# Analysis 3a: A3b plus land-use x trophic niche and land-use x
#              body mass interactions (NO phylogeny, same rationale).
#
# A3a is a clean nested SUPERSET of the A3b covariate model:
#
#   A3b covariate : Biome4 * Predominant_simple + Trophic.Niche + z_logMass
#   A3a (this)    : A3b covariate + Predominant_simple:z_logMass
#                                 + Predominant_simple:Trophic.Niche
#
# Same rows, priors, sampler and random structure as A3b (both call
# load_community_data()), so the A3a_03 LOO directly tests whether
# guild- and size-specific land-use responses add predictive value.
# Biome x land use is reported from A3b; A3a contributes the
# trophic/mass story.
# ============================================================

source("scripts/00_config.R")
library(tidyverse)
library(brms)
library(cmdstanr)

dat <- load_community_data()

responses         <- RESP_COONEY
response_families <- RESP_FAMILY[responses]

# Term order matches A3b for the shared biome x land-use interaction
# (Biome4<b>:Predominant_simple<lu>). The two added interactions are
# written Predominant_simple-first so their fixef names read
# Predominant_simple<lu>:Trophic.Niche<t> and
# Predominant_simple<lu>:z_logMass.
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
     file = "fits/A3a_model_setup.RData")
cat("\nA3a setup complete. Saved to fits/A3a_model_setup.RData\n")
