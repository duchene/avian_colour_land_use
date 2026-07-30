# ============================================================
# A1a_01_setup.R
# Analysis 1a: Community colour/dichromatism ~ biome x land-use
#              PLUS land-use x trophic niche and land-use x body mass
#              interactions (NO phylogeny — same rationale as A1b).
#
# A1a is the REFIT of the original A1 with A1b's random structure,
# data, collapse, priors and sampler. It is therefore a clean nested
# SUPERSET of the A1b covariate model:
#
#   A1b covariate : Biome4 * Predominant_simple + Trophic.Niche + z_logMass
#   A1a (this)    : A1b covariate + Predominant_simple:z_logMass
#                                 + Predominant_simple:Trophic.Niche
#
# The only difference from A1b is those two extra interaction families,
# so LOO (A1a_03) directly tests whether guild- and size-specific
# land-use responses add predictive value beyond A1b. Biome x land-use
# is still reported from A1b; A1a contributes the trophic/mass story.
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
nthreads <- max(2, floor(parallel::detectCores() / 4))
message("Detected cores: ", parallel::detectCores(), " | threads/chain: ", nthreads)

# ------------------------------------------------------------
# DATA  (identical prep to A1b_01_setup.R so the two models share
#        the exact same rows — required for a fair LOO comparison)
# ------------------------------------------------------------
dat <- read.csv("data/present.csv", stringsAsFactors = FALSE)

# Sparse-cell collapse: Temperate Open x Plantation forest (42 records)
# -> Secondary (see ANALYSIS_NOTES.md common design decisions).
n_coll <- sum(dat$Biome4 == "Temperate Open" &
              dat$Predominant_simple == "Plantation forest", na.rm = TRUE)
dat$Predominant_simple[dat$Biome4 == "Temperate Open" &
                       dat$Predominant_simple == "Plantation forest"] <- "Secondary"
cat("Collapsed", n_coll, "Temperate Open x Plantation records -> Secondary\n")

dat <- dat %>%
  mutate(
    Biome4             = relevel(factor(Biome4), ref = "Tropical Forest"),
    Predominant_simple = relevel(factor(Predominant_simple), ref = "Primary vegetation"),
    Trophic.Niche      = factor(Trophic.Niche),
    SS   = factor(SS),
    SSB  = factor(SSB),
    SSBS = factor(SSBS),
    z_logMass = as.numeric(scale(log(Mass)))
  )

cat("Rows:", nrow(dat), " | SS:", nlevels(dat$SS),
    " SSB:", nlevels(dat$SSB), " SSBS:", nlevels(dat$SSBS), "\n")

# ------------------------------------------------------------
# RESPONSES / FAMILIES  (as A1b)
# ------------------------------------------------------------
responses <- c("meancolcooney", "dichrocooney", "malecolcooney", "dichrodiff")
response_families <- c(meancolcooney = "lognormal", dichrocooney = "lognormal",
                       malecolcooney = "lognormal", dichrodiff = "gaussian")

# ------------------------------------------------------------
# FORMULA  (one variant per response: the full-interaction model)
# Term order matches A1b for the shared biome x land-use interaction
# (Biome4<b>:Predominant_simple<lu>); the two added interactions are
# written Predominant_simple-first so their fixef names read
# Predominant_simple<lu>:Trophic.Niche<t> and
# Predominant_simple<lu>:z_logMass (as in the original A1).
# ------------------------------------------------------------
make_formula <- function(resp) {
  bf(as.formula(paste0(
    resp, " ~ Biome4 * Predominant_simple + Trophic.Niche + z_logMass",
    " + Predominant_simple:z_logMass",
    " + Predominant_simple:Trophic.Niche",
    " + (1|SS) + (1|SSB) + (1|SSBS)")))
}
formulas <- setNames(lapply(responses, make_formula), responses)

# ------------------------------------------------------------
# PRIORS (identical to A1b)
# ------------------------------------------------------------
priors_lognormal <- c(
  prior(normal(0, 2),   class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma"))
priors_gaussian <- c(
  prior(normal(0, 50),  class = "Intercept"),
  prior(normal(0, 20),  class = "b"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma"))

dir.create("fits", showWarnings = FALSE)
save(dat, responses, response_families, formulas,
     priors_lognormal, priors_gaussian, nthreads,
     file = "fits/A1a_model_setup.RData")
cat("\nA1a setup complete. Saved to fits/A1a_model_setup.RData\n")
