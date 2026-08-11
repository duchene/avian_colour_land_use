# ============================================================
# A3_01_setup.R
# Analysis 3: Phylogenetic regression of species-level colour on
#             land-use association.
#
# One row per species. Colour is a species-level trait, so unlike
# A1b a phylogenetic random effect IS identifiable here: the
# covariance matrix constrains it and sigma captures the iid
# residual.
#
# Predictors are within-species land-use proportions (records of a
# species in land use j / total records of that species), which sum
# to 1 per species. The model therefore drops the intercept, and
# each coefficient reads as the expected colour of a species found
# exclusively in that land use.
#
# OPEN DECISION (flagged 2026-08-11, unchanged pending a call):
# proportions are computed from the FULL raw file, so they include
# studies with no primary vegetation and records in the dropped
# Temperate Open biome. Every other analysis works from the
# filtered present.csv. The rationale is that a species' habitat
# association is a property of the species, best estimated from all
# available records rather than the modelling subset. If that is not
# wanted, swap the raw read for present.csv. See ANALYSIS_NOTES.md.
# ============================================================

source("scripts/00_config.R")
library(tidyverse)
library(brms)
library(cmdstanr)
library(ape)

fulldat <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"),
                    stringsAsFactors = FALSE)
cat("Total records:", nrow(fulldat),
    "| species:", length(unique(fulldat$Best_guess_binomial)), "\n")

# ------------------------------------------------------------
# WITHIN-SPECIES LAND-USE PROPORTIONS
# ------------------------------------------------------------
lu_wide <- fulldat %>%
  count(Best_guess_binomial, Predominant_simple) %>%
  group_by(Best_guess_binomial) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = Predominant_simple, values_from = proportion,
              values_fill = 0, names_prefix = "prop_")
names(lu_wide) <- gsub(" ", "_", names(lu_wide))

# ------------------------------------------------------------
# SPECIES-LEVEL COLOUR (constant within species: take first non-NA)
# ------------------------------------------------------------
species_colour <- fulldat %>%
  group_by(Best_guess_binomial) %>%
  summarise(jetz_sp       = first(jetz_sp),
            meancolcooney = first(na.omit((LociUVS_male_cooney +
                                           LociUVS_female_cooney) / 2)),
            dichrocooney  = first(na.omit(LociUVS_male_cooney /
                                          LociUVS_female_cooney)),
            malecolcooney = first(na.omit(LociUVS_male_cooney)),
            dichrodiff    = first(na.omit(LociUVS_male_cooney -
                                          LociUVS_female_cooney)),
            Mass          = first(na.omit(Mass)),
            Trophic.Niche = first(na.omit(Trophic.Niche)),
            .groups = "drop")

spdat <- left_join(lu_wide, species_colour, by = "Best_guess_binomial")
cat("Species-level dataset:", nrow(spdat), "species\n")
for (v in RESP_COONEY)
  cat("  non-NA", v, ":", sum(!is.na(spdat[[v]])), "\n")

# ------------------------------------------------------------
# TREE MATCHING
# ------------------------------------------------------------
# Species absent from BBtree2 are dropped, never grafted onto a
# congener: grafting invents branch lengths and distorts A.
tree <- read.tree("data/BBtree2.tre")
cat("Tree loaded:", length(tree$tip.label), "tips\n")

sp    <- match_to_tree(spdat, tree)
spdat <- merge(spdat, sp, by = "Best_guess_binomial")

tree <- drop.tip(tree, setdiff(tree$tip.label, spdat$phylo))
A <- vcv.phylo(tree, corr = TRUE)
cat("Final dataset:", nrow(spdat), "species |", length(tree$tip.label), "tips\n")

# ------------------------------------------------------------
# RESPONSES, PREDICTORS, PRIORS
# ------------------------------------------------------------
responses         <- RESP_COONEY
response_families <- RESP_FAMILY[responses]
lu_predictors     <- grep("^prop_", names(spdat), value = TRUE)
cat("\nPredictors:", paste(lu_predictors, collapse = " + "), "\n")

priors_A3          <- model_priors("phylo_lognormal")
priors_A3_gaussian <- model_priors("phylo_gaussian")

dir.create("fits", showWarnings = FALSE)
save(spdat, responses, lu_predictors, priors_A3, priors_A3_gaussian,
     response_families, nthreads, A, tree,
     file = "fits/A3_model_setup.RData")
cat("\nA3 setup complete. Saved to fits/A3_model_setup.RData\n")
