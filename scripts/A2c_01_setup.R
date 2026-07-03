# ============================================================
# A2c_01_setup.R
# Analysis 2c: Relative abundance ~ colour x biome x land-use
#              (2-way interactions only; phylogenetic).
# Extends A2 by adding Biome4 (as main + 2-way interactions),
# the full PREDICTS random hierarchy, and a phylogenetic random
# effect (1|gr(phylo, cov=A)). The phylo term is identifiable
# here because the response (relative abundance) varies within
# species (see ANALYSIS_NOTES.md A1b complexity note for contrast).
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
# DATA + relative abundance per SSBS
# ------------------------------------------------------------
dat <- read.csv("data/present.csv", stringsAsFactors = FALSE)
ssbs_tot <- dat %>% group_by(SSBS) %>%
  summarise(tot = sum(Effort_corrected_measurement, na.rm = TRUE), .groups = "drop")
dat <- dat %>% left_join(ssbs_tot, by = "SSBS") %>%
  mutate(relabund = Effort_corrected_measurement / tot)

# Sparse-cell collapse: Temperate Open x Plantation forest -> Secondary
n_coll <- sum(dat$Biome4 == "Temperate Open" &
              dat$Predominant_simple == "Plantation forest", na.rm = TRUE)
dat$Predominant_simple[dat$Biome4 == "Temperate Open" &
                       dat$Predominant_simple == "Plantation forest"] <- "Secondary"
cat("Collapsed", n_coll, "Temperate Open x Plantation records -> Secondary\n")

dat <- dat %>%
  mutate(Biome4 = relevel(factor(Biome4), ref = "Tropical Forest"),
         Predominant_simple = relevel(factor(Predominant_simple), ref = "Primary vegetation"),
         SS = factor(SS), SSB = factor(SSB), SSBS = factor(SSBS)) %>%
  filter(!is.na(relabund) & relabund > 0)

# ------------------------------------------------------------
# PHYLOGENETIC MATCHING (reuse A3 3-stage logic, joined to records)
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

sp <- dat %>% distinct(Best_guess_binomial, jetz_sp)
sp$binom_u <- gsub(" ", "_", sp$Best_guess_binomial)
sp$phylo <- ifelse(sp$jetz_sp %in% tree$tip.label, sp$jetz_sp, NA_character_)
sp$phylo <- ifelse(!is.na(sp$phylo), sp$phylo,
                   ifelse(sp$binom_u %in% tree$tip.label, sp$binom_u, NA_character_))
for (jn in names(synonyms))
  sp$phylo[is.na(sp$phylo) & sp$jetz_sp == jn] <- synonyms[[jn]]
sp <- sp %>% filter(!is.na(phylo))
dupes <- unique(sp$phylo[duplicated(sp$phylo)])      # tip claimed by >1 species
if (length(dupes)) { cat("dropping", length(dupes), "ambiguous tips\n")
                     sp <- sp %>% filter(!(phylo %in% dupes)) }

dat <- dat %>% inner_join(sp %>% select(Best_guess_binomial, phylo),
                          by = "Best_guess_binomial")
dat$phylo <- factor(dat$phylo)

# ------------------------------------------------------------
# STANDARDIZED COLOUR PREDICTORS
# Primary metrics; meancolcooney/dichrocooney are extensions.
# ------------------------------------------------------------
dat <- dat %>% mutate(z_malecolcooney = zscore(malecolcooney),
                      z_dichrodiff    = zscore(dichrodiff))
colour_vars <- c("z_malecolcooney", "z_dichrodiff")

# Phylogenetic covariance (correlation form); subset per-model at fit time.
tree <- drop.tip(tree, setdiff(tree$tip.label, unique(as.character(dat$phylo))))
A <- vcv.phylo(tree, corr = TRUE)

cat("Records:", nrow(dat), " | species/tips:", length(tree$tip.label),
    " | SS:", nlevels(droplevels(dat$SS)), "\n")
print(table(dat$Biome4, dat$Predominant_simple))

# ------------------------------------------------------------
# FORMULAS (all 2-way interactions; NO 3-way) + priors
# ------------------------------------------------------------
formulas <- setNames(lapply(colour_vars, function(cv) bf(as.formula(paste0(
  "relabund ~ ", cv, " * Biome4 + ", cv, " * Predominant_simple + Biome4 * Predominant_simple",
  " + (1|SS) + (1|SSB) + (1|SSBS) + (1 | gr(phylo, cov = A))")))), colour_vars)

priors_A2c <- c(prior(normal(0, 2),   class = "Intercept"),
                prior(normal(0, 0.5), class = "b"),
                prior(exponential(2), class = "sd"),
                prior(exponential(2), class = "sigma"))

save(dat, A, colour_vars, formulas, priors_A2c, nthreads,
     file = "fits/A2c_model_setup.RData")
cat("\nA2c setup complete. Saved to fits/A2c_model_setup.RData\n")
