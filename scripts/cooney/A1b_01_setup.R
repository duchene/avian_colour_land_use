# ============================================================
# A1b_01_setup.R
# Analysis 1b: Community colour/dichromatism ~ biome x land-use
#              (NO phylogeny — see ANALYSIS_NOTES.md "A1b note on
#               complexity": a per-species phylo term on a within-
#               species-constant response is non-identifiable.)
# Extends A1: single Biome4 x Predominant_simple interaction +
#             full PREDICTS random hierarchy (1|SS)+(1|SSB)+(1|SSBS).
# Two variants per response: base, and +Trophic.Niche +z_logMass.
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
nthreads <- max(2, floor(parallel::detectCores() / 4))
message("Detected cores: ", parallel::detectCores(), " | threads/chain: ", nthreads)

# ------------------------------------------------------------
# DATA
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
    SS   = factor(SS),
    SSB  = factor(SSB),
    SSBS = factor(SSBS),
    z_logMass = as.numeric(scale(log(Mass)))
  )

cat("Rows:", nrow(dat), " | SS:", nlevels(dat$SS),
    " SSB:", nlevels(dat$SSB), " SSBS:", nlevels(dat$SSBS), "\n")

# ------------------------------------------------------------
# RESPONSES / FAMILIES
# ------------------------------------------------------------
responses <- c("meancolcooney", "dichrocooney", "malecolcooney", "dichrodiff")
response_families <- c(meancolcooney = "lognormal", dichrocooney = "lognormal",
                       malecolcooney = "lognormal", dichrodiff = "gaussian")

# ------------------------------------------------------------
# FORMULAS  (base + covariate variant, per response)
# Note: Biome4 * Predominant_simple = both main effects + their
# 2-way interaction. The empty Temperate Open : Plantation cell
# yields a prior-only coefficient (excluded at the summary stage).
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
# PRIORS (match A1)
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

save(dat, responses, response_families, formulas,
     priors_lognormal, priors_gaussian, nthreads,
     file = "fits/A1b_model_setup.RData")
cat("\nA1b setup complete. Saved to fits/A1b_model_setup.RData\n")
