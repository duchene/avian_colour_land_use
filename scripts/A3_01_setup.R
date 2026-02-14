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
# For each species x land-use level: count of records in that
# land use / total records of that species. This gives within-
# species relative proportions (sum to 1 per species), removing
# the confound with overall species prevalence.
#
# Because they sum to 1, we drop Cropland as the reference level
# (matching A1) to avoid perfect multicollinearity.

lu_counts <- fulldat %>%
  count(Best_guess_binomial, Predominant_simple) %>%
  group_by(Best_guess_binomial) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

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

# Drop Cropland (reference level) to avoid compositional constraint
lu_wide <- lu_wide %>% select(-prop_Cropland)

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
# BBtree2 uses Jetz-era taxonomy but is a 9,072-tip subset.
# Many species in the data have valid Jetz names that are
# simply absent from this tree. We match in three stages:
#   1. jetz_sp column (primary)
#   2. Best_guess_binomial with underscores (fallback)
#   3. Known genus-level synonyms (same species, different
#      genus name due to taxonomic reclassification)
#
# We do NOT graft missing species onto congeners because
# that introduces arbitrary branch lengths that distort
# the phylogenetic covariance matrix.

# --- Stage 3 synonym table ---
# Each entry: jetz_sp name (not in tree) -> tree tip name (in tree).
# All are the same biological species under a reclassified genus.
synonyms <- c(
  Taeniopygia_bichenovii    = "Taenopygia_bichenovii",
  Macronous_ptilosus        = "Macronus_ptilosus",
  Macronous_gularis         = "Mixornis_gularis",
  Pygochelidon_cyanoleuca   = "Notiochelidon_cyanoleuca",
  Hodgsonius_phaenicuroides = "Luscinia_phaenicuroides",
  Conostoma_oemodium        = "Paradoxornis_aemodium",
  Hyloctistes_subulatus     = "Automolus_subulatus",
  Dioptrornis_fischeri      = "Melaenornis_fischeri",
  Trichastoma_bicolor       = "Pellorneum_bicolor",
  Trichastoma_celebense     = "Pellorneum_celebense",
  Trichastoma_rostratum     = "Pellorneum_rostratum",
  Speirops_lugubris         = "Zosterops_lugubris",
  Babax_lanceolatus         = "Pterorhinus_lanceolatus",
  Rhinomyias_umbratilis     = "Cyornis_umbratilis",
  Rhopocichla_atriceps      = "Dumetia_atriceps"
)

# Verify all synonym targets are in the tree
stopifnot(all(synonyms %in% tree$tip.label))

# Verify no synonym target is already claimed by another species
# (would create duplicate tips in the dataset)

# --- Stage 1: jetz_sp ---
spdat$binom_u <- gsub(" ", "_", spdat$Best_guess_binomial)

spdat$phylo <- ifelse(
  spdat$jetz_sp %in% tree$tip.label, spdat$jetz_sp, NA
)
n_jetz <- sum(!is.na(spdat$phylo))

# --- Stage 2: Best_guess_binomial ---
spdat$phylo <- ifelse(
  !is.na(spdat$phylo), spdat$phylo,
  ifelse(spdat$binom_u %in% tree$tip.label, spdat$binom_u, NA)
)
n_binom <- sum(!is.na(spdat$phylo)) - n_jetz

# --- Stage 3: synonym lookup ---
unmatched_idx <- which(is.na(spdat$phylo))
syn_match <- synonyms[spdat$jetz_sp[unmatched_idx]]
syn_match <- syn_match[!is.na(syn_match)]

# Check no synonym target duplicates an already-matched tip
already_matched <- spdat$phylo[!is.na(spdat$phylo)]
dupes <- syn_match[syn_match %in% already_matched]
if (length(dupes) > 0) {
  cat("WARNING: dropping", length(dupes), "synonym matches that duplicate existing tips\n")
  syn_match <- syn_match[!(syn_match %in% already_matched)]
}

for (jetz_name in names(syn_match)) {
  idx <- which(spdat$jetz_sp == jetz_name & is.na(spdat$phylo))
  if (length(idx) == 1) spdat$phylo[idx] <- syn_match[[jetz_name]]
}
n_synonym <- sum(!is.na(spdat$phylo)) - n_jetz - n_binom

cat("\nTaxonomic matching:\n")
cat("  Stage 1 — jetz_sp:", n_jetz, "\n")
cat("  Stage 2 — Best_guess_binomial:", n_binom, "\n")
cat("  Stage 3 — synonym lookup:", n_synonym, "\n")
cat("  Total matched:", sum(!is.na(spdat$phylo)), "of", nrow(spdat), "\n")
cat("  Unmatched (dropped):", sum(is.na(spdat$phylo)), "\n")

# --- Write synonym table ---
synonym_table <- data.frame(
  jetz_sp    = names(synonyms),
  tree_tip   = unname(synonyms),
  used       = names(synonyms) %in% names(syn_match),
  stringsAsFactors = FALSE
)
dir.create("results", showWarnings = FALSE)
write.csv(synonym_table, "results/A3_synonym_table.csv", row.names = FALSE)
cat("Synonym table written to results/A3_synonym_table.csv\n")

# Drop unmatched species and helper column
spdat <- spdat %>% filter(!is.na(phylo)) %>% select(-binom_u)

# --- Drop duplicate phylo labels ---
# Some BirdLife species map to the same Jetz species (taxonomic lumping).
# Rather than merging ambiguous cases, we exclude both members of each
# duplicate pair to keep the dataset clean.
if (any(duplicated(spdat$phylo))) {
  dup_labels <- unique(spdat$phylo[duplicated(spdat$phylo)])
  cat("\nDropping", sum(spdat$phylo %in% dup_labels),
      "rows with ambiguous phylo labels:\n")
  for (d in dup_labels) {
    rows <- spdat %>% filter(phylo == d)
    cat("  ", d, ":", paste(rows$Best_guess_binomial, collapse = " / "), "\n")
  }
  spdat <- spdat %>% filter(!(phylo %in% dup_labels))
}

# Trim tree to matched species
tree <- drop.tip(tree, setdiff(tree$tip.label, spdat$phylo))
cat("Trimmed tree:", length(tree$tip.label), "tips\n")

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
