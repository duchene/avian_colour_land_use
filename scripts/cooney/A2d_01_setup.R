# ============================================================
# A2d_01_setup.R
# Analysis 2d: Paired-difference abundance change
#              ~ colour x biome x land use (2-way; phylogenetic).
#
# Each non-primary record is expressed as the species' relative
# abundance minus its own mean relative abundance in primary
# vegetation within the same study. That anchors every species to
# its own baseline and admits absences, which the presence-only
# A2c cannot do. Primary vegetation is absorbed into the response,
# so land use has four levels with Secondary as reference.
#
# Study unit is the canonical SS throughout, for both the
# qualifying filter and the baseline.
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
fulldat$malecolcooney <- fulldat$LociUVS_male_cooney
fulldat$dichrodiff    <- fulldat$LociUVS_male_cooney - fulldat$LociUVS_female_cooney

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
# STANDARDISED COLOUR PREDICTORS
# ------------------------------------------------------------
dat_final$z_malecolcooney <- zscore(dat_final$malecolcooney)
dat_final$z_dichrodiff    <- zscore(dat_final$dichrodiff)
colour_vars <- COLOUR_COONEY

cat("\nFinal A2d set:", nrow(dat_final), "records |",
    length(tree$tip.label), "species/tips | SS:",
    nlevels(droplevels(dat_final$SS)), "\n")

# ------------------------------------------------------------
# FORMULAS (all 2-way; NO 3-way — see A2f) + priors
# ------------------------------------------------------------
formulas <- setNames(lapply(colour_vars, function(cv) bf(as.formula(paste0(
  "diff_abund ~ ", cv, " * Biome4 + ", cv, " * Predominant_simple",
  " + Biome4 * Predominant_simple",
  " + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))")))), colour_vars)

priors_A2d <- model_priors("paired_diff")

dir.create("fits", showWarnings = FALSE)
save(dat_final, A, colour_vars, formulas, priors_A2d, nthreads,
     file = "fits/A2d_model_setup.RData")
cat("\nA2d setup complete. Saved to fits/A2d_model_setup.RData\n")
