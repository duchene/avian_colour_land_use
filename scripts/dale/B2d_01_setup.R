# ============================================================
# B2d_01_setup.R
# Analysis B2d: Dale-colour counterpart of A2d.
# Paired-difference abundance change ~ colour x biome x land-use
# (2-way; phylogenetic).
#
# Identical to A2d in data construction, phylogenetic matching, priors,
# sampler and random structure. The only difference is the colour
# predictors: z_malecoldale and z_dichrodiffdale (Dale), computed inline
# from the raw Dale plumage scores in place of the Cooney Loci scores.
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)
library(ape)

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
nthreads <- max(2, floor(parallel::detectCores() / 4))
message("Detected cores: ", parallel::detectCores(), " | threads/chain: ", nthreads)

zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# ------------------------------------------------------------
# LOAD FULL RAW DATA + Dale colour + Biome4
# ------------------------------------------------------------
fulldat <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"),
                    stringsAsFactors = FALSE)
fulldat$malecoldale    <- fulldat$Male_plumage_score_dale
fulldat$dichrodiffdale <- fulldat$Male_plumage_score_dale -
                          fulldat$Female_plumage_score_dale

# Biome4 — identical mapping to 00_data_preparation.R.
fulldat$Biome4 <- ifelse(
  fulldat$Biome %in% c("Tropical & Subtropical Moist Broadleaf Forests",
                       "Tropical & Subtropical Dry Broadleaf Forests",
                       "Tropical & Subtropical Coniferous Forests"), "Tropical Forest",
  ifelse(fulldat$Biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands",
                              "Mangroves"), "Tropical Open",
  ifelse(fulldat$Biome %in% c("Temperate Broadleaf & Mixed Forests",
                              "Temperate Conifer Forests",
                              "Mediterranean Forests, Woodlands & Scrub"), "Temperate Forest",
         "Temperate Open")))

# ------------------------------------------------------------
# QUALIFYING STUDIES (SS with Primary + >=1 other land use)
# ------------------------------------------------------------
qual_ss <- fulldat %>% group_by(SS) %>%
  summarise(has_primary = any(Predominant_simple == "Primary vegetation"),
            n_lu = n_distinct(Predominant_simple), .groups = "drop") %>%
  filter(has_primary, n_lu > 1) %>% pull(SS)
dat <- fulldat %>% filter(SS %in% qual_ss)
cat("After SS qualifying filter:", nrow(dat), "rows from", length(qual_ss), "studies\n")

# ------------------------------------------------------------
# RELATIVE ABUNDANCE per SSBS, then SS x species primary baseline
# ------------------------------------------------------------
ssbs_tot <- dat %>% group_by(SSBS) %>%
  summarise(tot = sum(Effort_corrected_measurement, na.rm = TRUE), .groups = "drop")
dat <- dat %>% left_join(ssbs_tot, by = "SSBS") %>%
  mutate(relabund = Effort_corrected_measurement / tot)

primary_baseline <- dat %>%
  filter(Predominant_simple == "Primary vegetation") %>%
  group_by(SS, Best_guess_binomial) %>%
  summarise(primary_mean = mean(relabund, na.rm = TRUE), .groups = "drop")
cat("Primary baselines:", nrow(primary_baseline), "SS x species combinations\n")

# ------------------------------------------------------------
# PAIRED DIFFERENCES (non-primary records minus baseline), drop zeros
# ------------------------------------------------------------
modified <- dat %>%
  filter(Predominant_simple != "Primary vegetation") %>%
  inner_join(primary_baseline, by = c("SS", "Best_guess_binomial")) %>%
  mutate(diff_abund = relabund - primary_mean)

dat_final <- modified %>% filter(diff_abund != 0)
cat("Non-zero paired differences:", nrow(dat_final),
    " (neg:", sum(dat_final$diff_abund < 0),
    " pos:", sum(dat_final$diff_abund > 0), ")\n")

# ------------------------------------------------------------
# SPARSE-CELL CHECK + collapse (data-driven, re-verified here)
# ------------------------------------------------------------
cat("\nBiome4 x Predominant_simple (records):\n")
ct <- table(dat_final$Biome4, dat_final$Predominant_simple)
print(ct)
sparse <- which(ct > 0 & ct < 50, arr.ind = TRUE)
if (nrow(sparse)) {
  cat("\nNON-EMPTY cells with <50 records (review):\n")
  for (i in seq_len(nrow(sparse)))
    cat("  ", rownames(ct)[sparse[i,1]], "x", colnames(ct)[sparse[i,2]],
        ":", ct[sparse[i,1], sparse[i,2]], "\n")
}
to_pl <- sum(dat_final$Biome4 == "Temperate Open" &
             dat_final$Predominant_simple == "Plantation forest")
if (to_pl > 0 && to_pl < 50) {
  dat_final$Predominant_simple[dat_final$Biome4 == "Temperate Open" &
    dat_final$Predominant_simple == "Plantation forest"] <- "Secondary"
  cat("\nCollapsed Temperate Open x Plantation (", to_pl, "records) -> Secondary\n")
} else {
  cat("\nTemperate Open x Plantation =", to_pl,
      "records; no collapse applied (not <50).\n")
}

# ------------------------------------------------------------
# FACTORS  (Primary absent -> 4 land-use levels; ref = Secondary)
# ------------------------------------------------------------
dat_final <- dat_final %>%
  mutate(Biome4 = relevel(factor(Biome4), ref = "Tropical Forest"),
         Predominant_simple = relevel(factor(Predominant_simple), ref = "Secondary"),
         SS = factor(SS), SSB = factor(SSB), SSBS = factor(SSBS))

# ------------------------------------------------------------
# PHYLOGENETIC MATCHING (reuse A3 3-stage logic)
# ------------------------------------------------------------
tree <- read.tree("data/BBtree2.tre")
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
  Rhopocichla_atriceps      = "Dumetia_atriceps")
stopifnot(all(synonyms %in% tree$tip.label))

sp <- dat_final %>% distinct(Best_guess_binomial, jetz_sp)
sp$binom_u <- gsub(" ", "_", sp$Best_guess_binomial)
sp$phylo <- ifelse(sp$jetz_sp %in% tree$tip.label, sp$jetz_sp, NA_character_)
sp$phylo <- ifelse(!is.na(sp$phylo), sp$phylo,
                   ifelse(sp$binom_u %in% tree$tip.label, sp$binom_u, NA_character_))
for (jn in names(synonyms))
  sp$phylo[is.na(sp$phylo) & sp$jetz_sp == jn] <- synonyms[[jn]]
sp <- sp %>% filter(!is.na(phylo))
dupes <- unique(sp$phylo[duplicated(sp$phylo)])
if (length(dupes)) { cat("dropping", length(dupes), "ambiguous tips\n")
                     sp <- sp %>% filter(!(phylo %in% dupes)) }

dat_final <- dat_final %>% inner_join(sp %>% select(Best_guess_binomial, phylo),
                                      by = "Best_guess_binomial")
dat_final$phylo <- factor(dat_final$phylo)

# standardized colour predictors (Dale)
dat_final <- dat_final %>% mutate(z_malecoldale   = zscore(malecoldale),
                                  z_dichrodiffdale = zscore(dichrodiffdale))
colour_vars <- c("z_malecoldale", "z_dichrodiffdale")

tree <- drop.tip(tree, setdiff(tree$tip.label, unique(as.character(dat_final$phylo))))
A <- vcv.phylo(tree, corr = TRUE)

cat("\nFinal B2d set:", nrow(dat_final), "records |",
    length(tree$tip.label), "species/tips | SS:",
    nlevels(droplevels(dat_final$SS)), "\n")

# ------------------------------------------------------------
# FORMULAS (all 2-way; NO 3-way) + Gaussian priors (as A2d)
# ------------------------------------------------------------
formulas <- setNames(lapply(colour_vars, function(cv) bf(as.formula(paste0(
  "diff_abund ~ ", cv, " * Biome4 + ", cv, " * Predominant_simple + Biome4 * Predominant_simple",
  " + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))")))), colour_vars)

priors_B2d <- c(prior(normal(0, 0.5), class = "Intercept"),
                prior(normal(0, 0.2), class = "b"),
                prior(exponential(5), class = "sd"),
                prior(exponential(5), class = "sigma"))

save(dat_final, A, colour_vars, formulas, priors_B2d, nthreads,
     file = "fits/B2d_model_setup.RData")
cat("\nB2d setup complete. Saved to fits/B2d_model_setup.RData\n")
