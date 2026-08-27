# ============================================================
# A2c_01_setup.R
# Analysis 2c: Relative abundance ~ colour x biome x land use
#              (all 2-way interactions, NO 3-way; phylogenetic).
#
# The phylogenetic term is identifiable here because the response
# (relative abundance) varies within species, unlike colour in A3b.
#
# Phylogeny-colour confound: colour is strongly conserved, so the
# phylogenetic random effect competes with colour for the same
# species-level variance. Colour fixed effects should be read as
# effects BEYOND phylogeny. Biome and land use vary within species
# and are not subject to this.
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
# STANDARDISED COLOUR PREDICTORS
# ------------------------------------------------------------
dat$z_malecolcooney <- zscore(dat$malecolcooney)
dat$z_dichrodiff    <- zscore(dat$dichrodiff)
colour_vars <- COLOUR_COONEY

cat("Records:", nrow(dat), " | species/tips:", length(tree$tip.label),
    " | SS:", nlevels(droplevels(dat$SS)), "\n")
check_cells(dat)

# ------------------------------------------------------------
# FORMULAS (all 2-way interactions; NO 3-way — see A2e) + priors
# ------------------------------------------------------------
formulas <- setNames(lapply(colour_vars, function(cv) bf(as.formula(paste0(
  "relabund ~ ", cv, " * Biome4 + ", cv, " * Predominant_simple",
  " + Biome4 * Predominant_simple",
  " + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))")))), colour_vars)

priors_A2c <- model_priors("lognormal")

dir.create("fits", showWarnings = FALSE)
save(dat, A, colour_vars, formulas, priors_A2c, nthreads,
     file = "fits/A2c_model_setup.RData")
cat("\nA2c setup complete. Saved to fits/A2c_model_setup.RData\n")
