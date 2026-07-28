# ============================================================
# B1b_01_setup.R
# Analysis B1b: Dale-colour counterpart of A1b.
# Community colour/dichromatism ~ biome x land-use (NO phylogeny),
# full PREDICTS random hierarchy (1|SS)+(1|SSB)+(1|SSBS).
# Two variants per response: base, and +Trophic.Niche +z_logMass.
#
# Identical to A1b except the colour source: the four Dale metrics
# in place of the Cooney UVS metrics. B1b is the primary community
# model of the Dale (B) family and the reduced model in the B1a LOO.
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
nthreads <- max(2, floor(parallel::detectCores() / 4))
message("Detected cores: ", parallel::detectCores(), " | threads/chain: ", nthreads)

# ------------------------------------------------------------
# DATA  (identical prep to A1b_01_setup.R)
# ------------------------------------------------------------
dat <- read.csv("data/present.csv", stringsAsFactors = FALSE)

n_coll <- sum(dat$Biome4 == "Temperate Open" &
              dat$Predominant_simple == "Plantation forest", na.rm = TRUE)
dat$Predominant_simple[dat$Biome4 == "Temperate Open" &
                       dat$Predominant_simple == "Plantation forest"] <- "Secondary"
cat("Collapsed", n_coll, "Temperate Open x Plantation records -> Secondary\n")

dat <- dat %>%
  mutate(
    Biome4             = relevel(factor(Biome4), ref = "Tropical Forest"),
    Predominant_simple = relevel(factor(Predominant_simple), ref = "Primary vegetation"),
    SS   = factor(SS),
    SSB  = factor(SSB),
    SSBS = factor(SSBS),
    z_logMass = as.numeric(scale(log(Mass)))
  )

cat("Rows:", nrow(dat), " | SS:", nlevels(dat$SS),
    " SSB:", nlevels(dat$SSB), " SSBS:", nlevels(dat$SSBS), "\n")

# ------------------------------------------------------------
# RESPONSES / FAMILIES  (Dale metrics)
# ------------------------------------------------------------
responses <- c("meancoldale", "dichrodale", "malecoldale", "dichrodiffdale")
response_families <- c(meancoldale = "lognormal", dichrodale = "lognormal",
                       malecoldale = "lognormal", dichrodiffdale = "gaussian")

# ------------------------------------------------------------
# FORMULAS  (base + covariate variant, per response)
# ------------------------------------------------------------
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
     file = "fits/B1b_model_setup.RData")
cat("\nB1b setup complete. Saved to fits/B1b_model_setup.RData\n")
