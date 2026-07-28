# ============================================================
# B3_01_setup.R
# Analysis B3: Dale-colour counterpart of A3.
# Phylogenetic regression: colour ~ land-use proportion scores
# (species-level).
#
# Identical to A3 in land-use proportions, taxonomic matching, tree
# handling, priors, families, and sampler. The only difference is the
# colour source: the four Dale metrics (species-level) in place of the
# Cooney UVS metrics.
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
cat("Total records:", nrow(fulldat), "\n")
cat("Unique species:", length(unique(fulldat$Best_guess_binomial)), "\n")

# ============================================================
# LAND-USE PROPORTION SCORES (species-level; identical to A3)
# ============================================================
lu_counts <- fulldat %>%
  count(Best_guess_binomial, Predominant_simple) %>%
  group_by(Best_guess_binomial) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

lu_wide <- lu_counts %>%
  select(-n) %>%
  pivot_wider(names_from = Predominant_simple, values_from = proportion,
              values_fill = 0, names_prefix = "prop_")
names(lu_wide) <- gsub(" ", "_", names(lu_wide))

# ============================================================
# SPECIES-LEVEL DALE COLOUR METRICS
# (first non-NA per species, as A3 does for Cooney)
# ============================================================
species_colour <- fulldat %>%
  group_by(Best_guess_binomial) %>%
  summarise(
    jetz_sp        = first(jetz_sp),
    meancoldale    = first(na.omit(
      (Male_plumage_score_dale + Female_plumage_score_dale) / 2)),
    dichrodale     = first(na.omit(
      Male_plumage_score_dale / Female_plumage_score_dale)),
    malecoldale    = first(na.omit(Male_plumage_score_dale)),
    dichrodiffdale = first(na.omit(
      Male_plumage_score_dale - Female_plumage_score_dale)),
    Mass           = first(na.omit(Mass)),
    Trophic.Niche  = first(na.omit(Trophic.Niche)),
    .groups = "drop"
  )

spdat <- left_join(lu_wide, species_colour, by = "Best_guess_binomial")
cat("\nSpecies-level dataset:", nrow(spdat), "species\n")
cat("Non-NA Dale colour: mean", sum(!is.na(spdat$meancoldale)),
    " male", sum(!is.na(spdat$malecoldale)),
    " ratio", sum(!is.na(spdat$dichrodale)),
    " diff", sum(!is.na(spdat$dichrodiffdale)), "\n")

# ============================================================
# LOAD AND MATCH PHYLOGENY (3-stage; identical to A3)
# ============================================================
tree <- read.tree("data/BBtree2.tre")
cat("\nTree loaded:", length(tree$tip.label), "tips\n")

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
stopifnot(all(synonyms %in% tree$tip.label))

spdat$binom_u <- gsub(" ", "_", spdat$Best_guess_binomial)
spdat$phylo <- ifelse(spdat$jetz_sp %in% tree$tip.label, spdat$jetz_sp, NA)
n_jetz <- sum(!is.na(spdat$phylo))
spdat$phylo <- ifelse(!is.na(spdat$phylo), spdat$phylo,
                      ifelse(spdat$binom_u %in% tree$tip.label, spdat$binom_u, NA))
n_binom <- sum(!is.na(spdat$phylo)) - n_jetz

unmatched_idx <- which(is.na(spdat$phylo))
syn_match <- synonyms[spdat$jetz_sp[unmatched_idx]]
syn_match <- syn_match[!is.na(syn_match)]
already_matched <- spdat$phylo[!is.na(spdat$phylo)]
syn_match <- syn_match[!(syn_match %in% already_matched)]
for (jetz_name in names(syn_match)) {
  idx <- which(spdat$jetz_sp == jetz_name & is.na(spdat$phylo))
  if (length(idx) == 1) spdat$phylo[idx] <- syn_match[[jetz_name]]
}
n_synonym <- sum(!is.na(spdat$phylo)) - n_jetz - n_binom

cat("\nTaxonomic matching: jetz", n_jetz, "| binomial", n_binom,
    "| synonym", n_synonym, "| total", sum(!is.na(spdat$phylo)),
    "of", nrow(spdat), "\n")

synonym_table <- data.frame(
  jetz_sp  = names(synonyms), tree_tip = unname(synonyms),
  used     = names(synonyms) %in% names(syn_match), stringsAsFactors = FALSE)
dir.create("results", showWarnings = FALSE)
write.csv(synonym_table, "results/B3_synonym_table.csv", row.names = FALSE)

spdat <- spdat %>% filter(!is.na(phylo)) %>% select(-binom_u)

if (any(duplicated(spdat$phylo))) {
  dup_labels <- unique(spdat$phylo[duplicated(spdat$phylo)])
  cat("\nDropping", sum(spdat$phylo %in% dup_labels),
      "rows with ambiguous phylo labels\n")
  spdat <- spdat %>% filter(!(phylo %in% dup_labels))
}

tree <- drop.tip(tree, setdiff(tree$tip.label, spdat$phylo))
cat("Trimmed tree:", length(tree$tip.label), "tips | species:", nrow(spdat), "\n")

A <- vcv.phylo(tree, corr = TRUE)

# ============================================================
# RESPONSES / PREDICTORS / PRIORS / FAMILIES
# ============================================================
responses <- c("meancoldale", "dichrodale", "malecoldale", "dichrodiffdale")
lu_predictors <- names(lu_wide)[-1]
cat("\nPredictors:", paste(lu_predictors, collapse = " + "), "\n")

priors_B3 <- c(
  prior(normal(0, 5), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)
priors_B3_gaussian <- c(
  prior(normal(0, 50), class = "b"),
  prior(exponential(0.05), class = "sd"),
  prior(exponential(0.05), class = "sigma")
)
response_families <- c(
  meancoldale    = "lognormal",
  dichrodale     = "lognormal",
  malecoldale    = "lognormal",
  dichrodiffdale = "gaussian"
)

save(spdat, responses, lu_predictors, priors_B3, priors_B3_gaussian,
     response_families, nthreads, A, tree,
     file = "fits/B3_model_setup.RData")
cat("\nB3 setup complete. Saved to fits/B3_model_setup.RData\n")
