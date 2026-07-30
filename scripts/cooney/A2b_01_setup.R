# ============================================================
# A2b_01_setup.R
# Analysis 2b: Paired-difference approach
# Response: species relative abundance at modified SSBS minus
#           mean relative abundance at primary-veg SSBS
#           (within same study), excluding zero differences.
# Predictors: colour (malecolcooney, dichrodiff) x land use.
# Gaussian family.
# ============================================================

library(dplyr)
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
nthreads <- max(2, floor(parallel::detectCores() / 4))

message("Detected cores: ", parallel::detectCores())
message("Threads per chain: ", nthreads)

# ============================================================
# LOAD FULL RAW DATA
# ============================================================
fulldat <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"))
cat("Full dataset rows:", nrow(fulldat), "\n")

# ============================================================
# COMPUTE COLOUR VARIABLES
# ============================================================
fulldat$malecolcooney <- fulldat$LociUVS_male_cooney
fulldat$dichrodiff    <- fulldat$LociUVS_male_cooney - fulldat$LociUVS_female_cooney

# ============================================================
# FILTER: studies with Primary vegetation + ≥1 other land use
# ============================================================
prim_refs <- fulldat %>%
  group_by(Reference) %>%
  summarise(
    has_primary = any(Predominant_simple == "Primary vegetation"),
    n_lu = n_distinct(Predominant_simple),
    .groups = "drop"
  ) %>%
  filter(has_primary, n_lu > 1) %>%
  pull(Reference)

dat <- fulldat %>%
  filter(Reference %in% prim_refs) %>%
  mutate(SSBS = factor(SSBS))

cat("After primary-veg study filter:", nrow(dat), "rows\n")

# ============================================================
# RELATIVE ABUNDANCE PER SSBS
# ============================================================
ssbs_totals <- dat %>%
  group_by(SSBS) %>%
  summarise(total_N_SSBS = sum(Effort_corrected_measurement, na.rm = TRUE),
            .groups = "drop")

dat <- dat %>%
  left_join(ssbs_totals, by = "SSBS") %>%
  mutate(relabund = Effort_corrected_measurement / total_N_SSBS)

# ============================================================
# COMPUTE PRIMARY-VEGETATION BASELINE PER SPECIES × STUDY
# ============================================================
# Mean relative abundance across all primary-veg SSBS within each study
primary_baseline <- dat %>%
  filter(Predominant_simple == "Primary vegetation") %>%
  group_by(Reference, Best_guess_binomial) %>%
  summarise(primary_mean = mean(relabund, na.rm = TRUE),
            n_primary_sites = n(),
            .groups = "drop")

cat("\nPrimary baselines computed for",
    nrow(primary_baseline), "species × study combinations\n")

# ============================================================
# PAIRED DIFFERENCES: modified SSBS minus primary baseline
# ============================================================
# Keep only non-primary records, join with primary baseline
modified <- dat %>%
  filter(Predominant_simple != "Primary vegetation") %>%
  inner_join(primary_baseline, by = c("Reference", "Best_guess_binomial")) %>%
  mutate(
    diff_abund = relabund - primary_mean,
    Predominant_simple = factor(Predominant_simple)
  )

cat("\nPaired-difference dataset before zero exclusion:", nrow(modified), "rows\n")
cat("  Exactly zero:", sum(modified$diff_abund == 0), "\n")
cat("  Non-zero:    ", sum(modified$diff_abund != 0), "\n")

# Exclude zero differences (species absent from both primary and modified)
dat_final <- modified %>%
  filter(diff_abund != 0)

cat("\nFinal dataset after excluding zeros:", nrow(dat_final), "rows\n")
cat("  Negative (lost abundance):", sum(dat_final$diff_abund < 0), "\n")
cat("  Positive (gained abundance):", sum(dat_final$diff_abund > 0), "\n")

cat("\nResponse summary (diff_abund):\n")
print(summary(dat_final$diff_abund))

cat("\nLand-use breakdown:\n")
print(table(dat_final$Predominant_simple))

cat("\nUnique studies:", n_distinct(dat_final$Reference), "\n")
cat("Unique SSBS:", n_distinct(dat_final$SSBS), "\n")
cat("Unique species:", n_distinct(dat_final$Best_guess_binomial), "\n")

# ============================================================
# STANDARDIZE COLOUR PREDICTORS
# ============================================================
dat_final <- dat_final %>%
  mutate(
    z_malecolcooney = as.numeric(scale(malecolcooney)),
    z_dichrodiff    = as.numeric(scale(dichrodiff))
  )

colour_vars <- c("z_malecolcooney", "z_dichrodiff")

cat("\nSample sizes (non-NA diff_abund + colour):\n")
for (cv in colour_vars) {
  n <- sum(!is.na(dat_final$diff_abund) & !is.na(dat_final[[cv]]))
  cat("  ", cv, ":", n, "\n")
}

# ============================================================
# MODEL FORMULAS
# ============================================================
# diff_abund ~ colour * land-use + (1 | SSBS) + (1 | Reference)
# Reference random intercept accounts for study-level variation
# in how differences are distributed.
formulas <- setNames(
  lapply(colour_vars, function(cv) {
    bf(as.formula(paste0(
      "diff_abund ~ ", cv, " * Predominant_simple + (1 | SSBS) + (1 | Reference)"
    )))
  }),
  colour_vars
)

# ============================================================
# PRIORS (Gaussian family — response on relative-abundance scale)
# ============================================================
priors_A2b <- c(
  prior(normal(0, 0.5),  class = "Intercept"),
  prior(normal(0, 0.2),  class = "b"),
  prior(exponential(5),  class = "sd"),
  prior(exponential(5),  class = "sigma")
)

# ============================================================
# SAVE
# ============================================================
save(dat_final, colour_vars, formulas, priors_A2b, nthreads,
     file = "fits/A2b_model_setup.RData")

cat("\nA2b setup complete. Saved to fits/A2b_model_setup.RData\n")
