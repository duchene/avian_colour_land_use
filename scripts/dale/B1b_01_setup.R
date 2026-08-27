# ============================================================
# B1b_01_setup.R
# Dale-colour counterpart of A1.
#
# Species-level phylogenetic regression of colour on within-species
# land-use proportions. Identical to A1 in construction, matching
# and priors; only the four colour metrics differ. Dale coverage is
# broader than Cooney, so more species survive the colour-NA filter.
#
# Carries the same OPEN DECISION as A1: proportions come from the
# full raw file rather than present.csv. See A1_01_setup.R.
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
# THE ONE DIFFERENCE FROM B1, matching A1b. B1 counts every row of the
# raw extract, and 84% of those rows are structural zeros from the
# near-complete site by species matrix, so the score reflects how a study
# spread its sites rather than where the bird was found. Built from
# detections only here.
detdat <- fulldat %>% filter(Effort_corrected_measurement != 0)
cat("Detections:", nrow(detdat), "of", nrow(fulldat), "rows |",
    length(unique(detdat$Best_guess_binomial)), "species\n")

lu_wide <- detdat %>%
  count(Best_guess_binomial, Predominant_simple) %>%
  group_by(Best_guess_binomial) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = Predominant_simple, values_from = proportion,
              values_fill = 0, names_prefix = "prop_")
names(lu_wide) <- gsub(" ", "_", names(lu_wide))

# ------------------------------------------------------------
# SPECIES-LEVEL COLOUR (Dale)
# ------------------------------------------------------------
species_colour <- detdat %>%
  group_by(Best_guess_binomial) %>%
  summarise(jetz_sp        = first(jetz_sp),
            meancoldale    = first(na.omit((Male_plumage_score_dale +
                                            Female_plumage_score_dale) / 2)),
            dichrodale     = first(na.omit(Male_plumage_score_dale /
                                           Female_plumage_score_dale)),
            malecoldale    = first(na.omit(Male_plumage_score_dale)),
            dichrodiffdale = first(na.omit(Male_plumage_score_dale -
                                           Female_plumage_score_dale)),
            Mass           = first(na.omit(Mass)),
            Trophic.Niche  = first(na.omit(Trophic.Niche)),
            .groups = "drop")

spdat <- left_join(lu_wide, species_colour, by = "Best_guess_binomial")
cat("Species-level dataset:", nrow(spdat), "species\n")
for (v in RESP_DALE)
  cat("  non-NA", v, ":", sum(!is.na(spdat[[v]])), "\n")

# ------------------------------------------------------------
# TREE MATCHING
# ------------------------------------------------------------
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
responses         <- RESP_DALE
response_families <- RESP_FAMILY[responses]
lu_predictors     <- grep("^prop_", names(spdat), value = TRUE)
cat("\nPredictors:", paste(lu_predictors, collapse = " + "), "\n")

priors_B1b          <- model_priors("phylo_lognormal")
priors_B1bb_gaussian <- model_priors("phylo_gaussian")

dir.create("fits", showWarnings = FALSE)
save(spdat, responses, lu_predictors, priors_B1b, priors_B1bb_gaussian,
     response_families, nthreads, A, tree,
     file = "fits/B1b_model_setup.RData")
cat("\nB1 setup complete. Saved to fits/B1b_model_setup.RData\n")
