# ============================================================
# A3_01_setup.R
# Analysis 3: Phylogenetic regression
# Colour ~ land-use proportion scores (species-level)
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)
library(ape)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
nthreads <- max(2, floor(parallel::detectCores() / 4))

# ============================================================
# LOAD RAW DATA
# ============================================================

fulldat <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"),
                    stringsAsFactors = FALSE)

n_total <- nrow(fulldat)
cat("Total records:", n_total, "\n")
cat("Unique species:", length(unique(fulldat$Best_guess_binomial)), "\n")

# ============================================================
# COMPUTE LAND-USE PROPORTION SCORES (species-level)
# ============================================================
# For each species x land-use level: count of records / total dataset size

lu_counts <- fulldat %>%
  count(Best_guess_binomial, Predominant_simple) %>%
  mutate(proportion = n / n_total)

# Pivot to one column per land-use level
lu_wide <- lu_counts %>%
  select(-n) %>%
  pivot_wider(
    names_from = Predominant_simple,
    values_from = proportion,
    values_fill = 0,
    names_prefix = "prop_"
  )

# Clean column names (remove spaces)
names(lu_wide) <- gsub(" ", "_", names(lu_wide))

cat("\nLand-use proportion columns:\n")
cat(names(lu_wide)[-1], sep = "\n")

# ============================================================
# COMPUTE SPECIES-LEVEL COLOUR METRICS
# ============================================================
# Colour is a species-level trait (same value per species),
# so take first non-NA value per species.

species_colour <- fulldat %>%
  group_by(Best_guess_binomial) %>%
  summarise(
    jetz_sp        = first(jetz_sp),
    meancolcooney  = first(na.omit(
      (LociUVS_male_cooney + LociUVS_female_cooney) / 2
    )),
    dichrocooney   = first(na.omit(
      LociUVS_male_cooney / LociUVS_female_cooney
    )),
    malecolcooney  = first(na.omit(LociUVS_male_cooney)),
    Mass           = first(na.omit(Mass)),
    Trophic.Niche  = first(na.omit(Trophic.Niche)),
    .groups = "drop"
  )

# ============================================================
# MERGE INTO SPECIES-LEVEL DATASET
# ============================================================

spdat <- left_join(lu_wide, species_colour, by = "Best_guess_binomial")

cat("\nSpecies-level dataset:", nrow(spdat), "species\n")

# Check completeness
cat("\nNon-NA counts:\n")
cat("  meancolcooney:", sum(!is.na(spdat$meancolcooney)), "\n")
cat("  dichrocooney:", sum(!is.na(spdat$dichrocooney)), "\n")
cat("  malecolcooney:", sum(!is.na(spdat$malecolcooney)), "\n")
cat("  Mass:", sum(!is.na(spdat$Mass)), "\n")

# ============================================================
# LOAD AND TRIM PHYLOGENETIC TREE
# ============================================================
# BigBirdTree (BBtree2) from:
# https://github.com/evolucionario/BigBirdTree

tree <- read.tree("data/BBtree2.tre")
cat("\nTree loaded:", length(tree$tip.label), "tips\n")

# ============================================================
# TAXONOMIC MATCHING
# ============================================================
# Tree uses Jetz-era taxonomy. Match species using:
#   1. jetz_sp column (primary — matches Jetz taxonomy)
#   2. Best_guess_binomial with underscores (fallback)

spdat$binom_u <- gsub(" ", "_", spdat$Best_guess_binomial)

spdat$phylo <- ifelse(
  spdat$jetz_sp %in% tree$tip.label, spdat$jetz_sp,
  ifelse(spdat$binom_u %in% tree$tip.label, spdat$binom_u, NA)
)

n_matched <- sum(!is.na(spdat$phylo))
n_jetz    <- sum(spdat$jetz_sp %in% tree$tip.label)
n_binom   <- n_matched - n_jetz

cat("\nTaxonomic matching:\n")
cat("  Matched via jetz_sp:", n_jetz, "\n")
cat("  Additional via Best_guess_binomial:", n_binom, "\n")
cat("  Total matched:", n_matched, "of", nrow(spdat), "\n")
cat("  Unmatched (dropped):", sum(is.na(spdat$phylo)), "\n")

# Drop unmatched species
spdat <- spdat %>% filter(!is.na(phylo))

# Trim tree to matched species
tree <- drop.tip(tree, setdiff(tree$tip.label, spdat$phylo))
cat("Trimmed tree:", length(tree$tip.label), "tips\n")

# Drop helper column
spdat <- spdat %>% select(-binom_u)

cat("Final dataset:", nrow(spdat), "species\n")

# Phylogenetic covariance matrix (correlation form for brms)
A <- vcv.phylo(tree, corr = TRUE)

# ============================================================
# RESPONSE VARIABLES
# ============================================================

responses <- c("meancolcooney", "dichrocooney", "malecolcooney")

# Land-use proportion predictor names
lu_predictors <- names(lu_wide)[-1]  # all prop_ columns
cat("\nPredictors:", paste(lu_predictors, collapse = " + "), "\n")

# ============================================================
# PRIORS
# ============================================================

priors_A3 <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)

# ============================================================
# SAVE SETUP
# ============================================================

save(spdat, responses, lu_predictors, priors_A3, nthreads, A, tree,
     file = "fits/A3_model_setup.RData")

cat("\nA3 setup complete. Saved to fits/A3_model_setup.RData\n")
