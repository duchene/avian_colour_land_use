# ============================================================
# B2d_01_setup.R
# Dale-colour counterpart of A2d.
#
# Paired-difference abundance change ~ colour x biome x land use
# (2-way; phylogenetic). Identical to A2d in data construction,
# matching, priors and sampler; only the two standardised colour
# predictors differ.
# ============================================================

source("scripts/00_config.R")
library(tidyverse)
library(brms)
library(cmdstanr)
library(ape)

# ------------------------------------------------------------
# RAW DATA -> paired differences
# ------------------------------------------------------------
fulldat <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"),
                    stringsAsFactors = FALSE)
fulldat$malecoldale    <- fulldat$Male_plumage_score_dale
fulldat$dichrodiffdale <- fulldat$Male_plumage_score_dale -
                          fulldat$Female_plumage_score_dale

dat_final <- build_paired_differences(fulldat)
check_cells(dat_final)

# ------------------------------------------------------------
# PHYLOGENETIC MATCHING
# ------------------------------------------------------------
tree <- read.tree("data/BBtree2.tre")
sp   <- match_to_tree(dat_final, tree)
dat_final <- merge(dat_final, sp, by = "Best_guess_binomial")
dat_final$phylo <- factor(dat_final$phylo)

tree <- drop.tip(tree, setdiff(tree$tip.label, unique(as.character(dat_final$phylo))))
A <- vcv.phylo(tree, corr = TRUE)

# ------------------------------------------------------------
# STANDARDISED COLOUR PREDICTORS (Dale)
# ------------------------------------------------------------
dat_final$z_malecoldale    <- zscore(dat_final$malecoldale)
dat_final$z_dichrodiffdale <- zscore(dat_final$dichrodiffdale)
colour_vars <- COLOUR_DALE

cat("\nFinal B2d set:", nrow(dat_final), "records |",
    length(tree$tip.label), "species/tips | SS:",
    nlevels(droplevels(dat_final$SS)), "\n")

# ------------------------------------------------------------
# FORMULAS (all 2-way; NO 3-way) + priors
# ------------------------------------------------------------
formulas <- setNames(lapply(colour_vars, function(cv) bf(as.formula(paste0(
  "diff_abund ~ ", cv, " * Biome4 + ", cv, " * Predominant_simple",
  " + Biome4 * Predominant_simple",
  " + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))")))), colour_vars)

priors_B2d <- model_priors("paired_diff")

dir.create("fits", showWarnings = FALSE)
save(dat_final, A, colour_vars, formulas, priors_B2d, nthreads,
     file = "fits/B2d_model_setup.RData")
cat("\nB2d setup complete. Saved to fits/B2d_model_setup.RData\n")
