# ============================================================
# A1b_01_setup.R
# Analysis 1b: Community colour/dichromatism ~ biome x land use
#              (NO phylogeny — a per-species phylogenetic term on a
#               within-species-constant response is non-identifiable;
#               see ANALYSIS_NOTES.md.)
# Full PREDICTS random hierarchy (1|SS)+(1|SSB)+(1|SSBS).
# Two variants per response: base, and +Trophic.Niche +z_logMass.
#
# A1b is the primary community-colour model. Biome x land use is
# reported from here.
# ============================================================

source("scripts/00_config.R")
library(tidyverse)
library(brms)
library(cmdstanr)

dat <- load_community_data()
check_cells(dat)

responses         <- RESP_COONEY
response_families <- RESP_FAMILY[responses]

# Biome4 * Predominant_simple expands to both main effects plus their
# 2-way interaction. Trophic.Niche and z_logMass enter as main effects
# only; the land-use interactions with them are A1a's job.
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
     file = "fits/A1b_model_setup.RData")
cat("\nA1b setup complete. Saved to fits/A1b_model_setup.RData\n")
