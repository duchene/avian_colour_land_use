# ============================================================
# B2c_01_setup.R
# Dale-colour counterpart of A2c.
#
# Relative abundance ~ colour x biome x land use (2-way only;
# phylogenetic). Identical to A2c in data, matching, priors and
# sampler; only the two standardised colour predictors differ.
# ============================================================

source("scripts/00_config.R")
library(tidyverse)
library(brms)
library(cmdstanr)
library(ape)

dat <- load_community_data()
dat <- add_relabund(dat)

# ------------------------------------------------------------
# PHYLOGENETIC MATCHING
# ------------------------------------------------------------
tree <- read.tree("data/BBtree2.tre")
sp   <- match_to_tree(dat, tree)
dat  <- merge(dat, sp, by = "Best_guess_binomial")
dat$phylo <- factor(dat$phylo)

tree <- drop.tip(tree, setdiff(tree$tip.label, unique(as.character(dat$phylo))))
A <- vcv.phylo(tree, corr = TRUE)

# ------------------------------------------------------------
# STANDARDISED COLOUR PREDICTORS (Dale)
# ------------------------------------------------------------
dat$z_malecoldale    <- zscore(dat$malecoldale)
dat$z_dichrodiffdale <- zscore(dat$dichrodiffdale)
colour_vars <- COLOUR_DALE

cat("Records:", nrow(dat), " | species/tips:", length(tree$tip.label),
    " | SS:", nlevels(droplevels(dat$SS)), "\n")
check_cells(dat)

# ------------------------------------------------------------
# FORMULAS (all 2-way interactions; NO 3-way) + priors
# ------------------------------------------------------------
formulas <- setNames(lapply(colour_vars, function(cv) bf(as.formula(paste0(
  "relabund ~ ", cv, " * Biome4 + ", cv, " * Predominant_simple",
  " + Biome4 * Predominant_simple",
  " + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))")))), colour_vars)

priors_B2c <- model_priors("lognormal")

dir.create("fits", showWarnings = FALSE)
save(dat, A, colour_vars, formulas, priors_B2c, nthreads,
     file = "fits/B2c_model_setup.RData")
cat("\nB2c setup complete. Saved to fits/B2c_model_setup.RData\n")
