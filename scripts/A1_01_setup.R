# ============================================================
# A1_01_setup.R
# Analysis 1: Colour ~ Land-use
# Load packages, data, and define model formulas/priors
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)

# ============================================================
# EFFICIENCY SETTINGS
# ============================================================
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
nthreads <- max(2, floor(parallel::detectCores() / 4))

message("Detected cores: ", parallel::detectCores())
message("Threads per chain: ", nthreads)

# ============================================================
# LOAD AND PREPARE DATA
# ============================================================
dat <- read.csv("data/present.csv", stringsAsFactors = FALSE)

dat <- dat %>%
  mutate(
    Predominant_simple = factor(Predominant_simple),
    Habitat            = factor(Habitat),
    Biome              = factor(Biome),
    Trophic.Niche      = factor(Trophic.Niche),
    SSBS               = factor(SSBS)
  )

# Mass transformation (log is essential - ranges 4.9 to 644.3)
dat <- dat %>%
  mutate(
    logMass = log(Mass),
    z_logMass = as.numeric(scale(logMass))
  )

cat("Data loaded. Rows:", nrow(dat), "\n")

# ============================================================
# RESPONSE SPECIFICATIONS
# ============================================================
# Transform/family decisions:
# - meancolcooney: Lognormal (heavily right-skewed, all positive)
# - dichrocooney: Lognormal (positive ratio, mildly right-skewed)
# - malecolcooney: Lognormal (right-skewed, all positive)

responses <- c("meancolcooney", "dichrocooney", "malecolcooney")

response_specs <- tibble(
  response = responses,
  family = c("lognormal", "lognormal", "lognormal")
)

cat("\nResponse specifications:\n")
print(response_specs)

# ============================================================
# MODEL FORMULAS
# ============================================================

# REDUCED model (recommended - drops Biome interaction)
formula_reduced <- bf(
  y ~ Predominant_simple + z_logMass + Habitat + Biome + Trophic.Niche +
      Predominant_simple:z_logMass +
      Predominant_simple:Habitat +
      Predominant_simple:Trophic.Niche +
      (1 | SSBS)
)

# MINIMAL model (only Mass and Trophic.Niche interactions)
formula_minimal <- bf(
  y ~ Predominant_simple + z_logMass + Habitat + Biome + Trophic.Niche +
      Predominant_simple:z_logMass +
      Predominant_simple:Trophic.Niche +
      (1 | SSBS)
)

# MAIN effects only (baseline for comparison)
formula_main <- bf(
  y ~ Predominant_simple + z_logMass + Habitat + Biome + Trophic.Niche +
      (1 | SSBS)
)

# ============================================================
# PRIORS
# ============================================================

base_priors <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)

# Stronger priors for student-t models (helps convergence)
student_priors <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 0.3), class = "b"),
  prior(exponential(3), class = "sd"),
  prior(exponential(3), class = "sigma"),
  prior(gamma(2, 0.1), class = "nu")
)

# ============================================================
# SAVE SETUP
# ============================================================
save(dat, responses, response_specs,
     formula_reduced, formula_minimal, formula_main,
     base_priors, student_priors, nthreads,
     file = "fits/A1_model_setup.RData")

cat("\nSetup complete. Saved to fits/A1_model_setup.RData\n")
