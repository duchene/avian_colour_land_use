# ============================================================
# A2_01_setup.R
# Analysis 2: Abundance ~ Colour x Land-use
# Setup: packages, data, formulas, priors
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
    Predominant_simple = relevel(factor(Predominant_simple), ref = "Primary vegetation"),
    SSBS = factor(SSBS),
    abundance = Effort_corrected_measurement
  )

# Standardize colour predictors for comparable effect sizes
dat <- dat %>%
  mutate(
    z_meancolcooney = as.numeric(scale(meancolcooney)),
    z_dichrocooney  = as.numeric(scale(dichrocooney)),
    z_malecolcooney = as.numeric(scale(malecolcooney))
  )

cat("Data loaded. Rows:", nrow(dat), "\n")
cat("Abundance summary:\n")
print(summary(dat$abundance))

# ============================================================
# RESPONSE AND COLOUR PREDICTORS
# ============================================================

colour_vars <- c("z_meancolcooney", "z_dichrocooney", "z_malecolcooney")

cat("\nColour predictor sample sizes (non-NA abundance + colour):\n")
for (cv in colour_vars) {
  n <- sum(!is.na(dat$abundance) & !is.na(dat[[cv]]))
  cat("  ", cv, ":", n, "\n")
}

# ============================================================
# MODEL FORMULAS
# ============================================================

# One formula per colour variable: abundance ~ colour * land-use + (1|SSBS)
formulas <- setNames(
  lapply(colour_vars, function(cv) {
    bf(as.formula(paste0("abundance ~ ", cv, " * Predominant_simple + (1 | SSBS)")))
  }),
  colour_vars
)

# ============================================================
# PRIORS
# ============================================================

# Lognormal family: priors on log-scale
priors_A2 <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)

# ============================================================
# SAVE SETUP
# ============================================================
save(dat, colour_vars, formulas, priors_A2, nthreads,
     file = "fits/A2_model_setup.RData")

cat("\nA2 setup complete. Saved to fits/A2_model_setup.RData\n")
