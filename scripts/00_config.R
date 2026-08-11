# ============================================================
# 00_config.R  —  single source of truth for shared decisions
# ============================================================
# Sourced first by every other script. Nothing here reads data or
# fits a model. It holds the choices that MUST be identical across
# A1a, A1b, A2c, A2d, A2e, A2f, A3 and their five B counterparts.
#
# RULE: change a decision here and it changes everywhere. Never
# re-implement any of these blocks inside an analysis script. If
# you find yourself copying a biome map, a synonym list, a prior
# set or a factor relevel into a setup script, it belongs here.
# ============================================================

# ------------------------------------------------------------
# Sampler backend and threading
# ------------------------------------------------------------
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
nthreads <- max(2, floor(parallel::detectCores() / 4))

# ------------------------------------------------------------
# Biome policy
# ------------------------------------------------------------
# The 11 WWF biomes present in the raw data collapse to 4 groups.
# Temperate Open is then DROPPED (decision 2026-08-11): on the
# corrected dataset it holds 115 primary records against 36 / 42 /
# 0 / 35 in the four modified land uses, Temperate Open x Cropland
# is empty, and no collapse leaves a defensible contrast. Dropping
# it costs ~0.9% of records and makes biome x land use fully
# estimable with no special-casing anywhere. See ANALYSIS_NOTES.md.
#
# Consequence: Biome4 carries THREE levels, not four.

BIOME_REF  <- "Tropical Forest"
BIOME_DROP <- "Temperate Open"

BIOME_MAP <- c(
  "Tropical & Subtropical Moist Broadleaf Forests"           = "Tropical Forest",
  "Tropical & Subtropical Dry Broadleaf Forests"             = "Tropical Forest",
  "Tropical & Subtropical Coniferous Forests"                = "Tropical Forest",
  "Tropical & Subtropical Grasslands, Savannas & Shrublands" = "Tropical Open",
  "Mangroves"                                                = "Tropical Open",
  "Temperate Broadleaf & Mixed Forests"                      = "Temperate Forest",
  "Temperate Conifer Forests"                                = "Temperate Forest",
  # Mediterranean sclerophyll is genuinely borderline (its WWF name spans
  # "Forests" and "Woodlands & Scrub"). Assigned to Temperate Forest.
  "Mediterranean Forests, Woodlands & Scrub"                 = "Temperate Forest",
  "Temperate Grasslands, Savannas & Shrublands"              = "Temperate Open",
  "Montane Grasslands & Shrublands"                          = "Temperate Open",
  "Tundra"                                                   = "Temperate Open")

# Map raw Biome to Biome4. Errors on any unlisted biome rather than
# silently bucketing it, which is how the old nested-ifelse version
# could quietly reassign records.
assign_biome4 <- function(biome) {
  unknown <- setdiff(unique(biome[!is.na(biome)]), names(BIOME_MAP))
  if (length(unknown))
    stop("Biome values absent from BIOME_MAP in 00_config.R: ",
         paste(unknown, collapse = " | "))
  unname(BIOME_MAP[biome])
}

# Assign Biome4, drop BIOME_DROP, and set the reference level.
# Every setup script calls this exactly once and does nothing else
# to Biome4.
apply_biome_policy <- function(df, verbose = TRUE) {
  if (!"Biome4" %in% names(df)) df$Biome4 <- assign_biome4(df$Biome)
  n_before <- nrow(df)
  df <- df[!is.na(df$Biome4) & df$Biome4 != BIOME_DROP, , drop = FALSE]
  if (verbose)
    cat("Biome policy: dropped", n_before - nrow(df), "records in",
        BIOME_DROP, "->", nrow(df), "retained\n")
  df$Biome4 <- relevel(factor(df$Biome4), ref = BIOME_REF)
  stopifnot(nlevels(df$Biome4) == 3)
  df
}

# ------------------------------------------------------------
# Land-use policy
# ------------------------------------------------------------
LU_REF <- "Primary vegetation"
LU_LEVELS <- c("Primary vegetation", "Secondary", "Plantation forest",
               "Cropland", "Pasture")

set_landuse <- function(df) {
  df$Predominant_simple <- relevel(factor(df$Predominant_simple), ref = LU_REF)
  df
}

# Guard against a biome x land-use cell too thin to estimate. With
# Temperate Open dropped no cell should trip this; it fires loudly
# if the data change again.
check_cells <- function(df, min_records = 50, min_species = 15) {
  ct <- table(df$Biome4, df$Predominant_simple)
  sp <- tapply(df$Best_guess_binomial,
               list(df$Biome4, df$Predominant_simple),
               function(x) length(unique(x)))
  cat("\nBiome4 x land use, records:\n"); print(ct)
  cat("\nBiome4 x land use, species:\n"); print(sp)
  sparse <- which(ct < min_records | is.na(sp) | sp < min_species, arr.ind = TRUE)
  if (nrow(sparse)) {
    warning("Sparse biome x land-use cells (<", min_records, " records or <",
            min_species, " species):\n",
            paste(rownames(ct)[sparse[, 1]], "x", colnames(ct)[sparse[, 2]],
                  collapse = "\n"), call. = FALSE)
  } else {
    cat("\nAll biome x land-use cells estimable.\n")
  }
  invisible(ct)
}

# ------------------------------------------------------------
# Shared data loaders
# ------------------------------------------------------------
# Record-level analytical data with every shared policy applied.
# Used by A1a / A1b / A2c and their B counterparts, so those models
# are guaranteed to see identical rows (required for a fair LOO).
# apply_biome_policy is idempotent: present.csv already excludes
# BIOME_DROP, so this only relevels.
load_community_data <- function(path = "data/present.csv", verbose = TRUE) {
  dat <- read.csv(path, stringsAsFactors = FALSE)
  dat <- apply_biome_policy(dat, verbose = verbose)
  dat <- set_landuse(dat)
  dat$Trophic.Niche <- factor(dat$Trophic.Niche)
  dat$SS   <- factor(dat$SS)
  dat$SSB  <- factor(dat$SSB)
  dat$SSBS <- factor(dat$SSBS)
  dat$z_logMass <- as.numeric(scale(log(dat$Mass)))
  if (verbose)
    cat("Rows:", nrow(dat), "| SS:", nlevels(dat$SS),
        " SSB:", nlevels(dat$SSB), " SSBS:", nlevels(dat$SSBS), "\n")
  dat
}

# Relative abundance per site: a species' effort-corrected count
# over the site total. Used by A2c / B2c.
add_relabund <- function(dat) {
  tot <- tapply(dat$Effort_corrected_measurement, dat$SSBS, sum, na.rm = TRUE)
  dat$relabund <- dat$Effort_corrected_measurement / tot[as.character(dat$SSBS)]
  dat[!is.na(dat$relabund) & dat$relabund > 0, , drop = FALSE]
}

# Paired-difference response for A2d / B2d: each non-primary record's
# relative abundance minus that species' mean relative abundance in
# primary vegetation within the same study.
#
# Built from the RAW file, absences included, because the design needs
# species that disappeared (negative difference) or appeared (positive).
# Differences of exactly zero are uninformative (species absent from
# both) and are dropped. The qualifying filter therefore runs on the
# raw rows here, unlike 00_data_preparation.R which filters presences.
# Primary vegetation is absorbed into the baseline, so the land-use
# factor has four levels with Secondary as reference.
build_paired_differences <- function(fulldat, verbose = TRUE) {
  fulldat <- apply_biome_policy(fulldat, verbose = verbose)

  qual <- tapply(fulldat$Predominant_simple, fulldat$SS,
                 function(x) any(x == LU_REF) && length(unique(x)) > 1)
  keep_ss <- names(qual)[qual]
  dat <- fulldat[fulldat$SS %in% keep_ss, , drop = FALSE]
  if (verbose)
    cat("Qualifying SS filter:", nrow(dat), "rows from", length(keep_ss), "studies\n")

  tot <- tapply(dat$Effort_corrected_measurement, dat$SSBS, sum, na.rm = TRUE)
  dat$relabund <- dat$Effort_corrected_measurement / tot[as.character(dat$SSBS)]

  prim <- dat[dat$Predominant_simple == LU_REF, ]
  key  <- paste(prim$SS, prim$Best_guess_binomial, sep = "\r")
  base <- tapply(prim$relabund, key, mean, na.rm = TRUE)
  if (verbose) cat("Primary baselines:", length(base), "SS x species combinations\n")

  mod <- dat[dat$Predominant_simple != LU_REF, , drop = FALSE]
  mod$primary_mean <- base[paste(mod$SS, mod$Best_guess_binomial, sep = "\r")]
  mod <- mod[!is.na(mod$primary_mean), , drop = FALSE]
  mod$diff_abund <- mod$relabund - mod$primary_mean

  # NB the !is.na guard is load-bearing: base-R `[` propagates NA rows
  # where dplyr::filter() silently drops them. relabund is NA wherever
  # Effort_corrected_measurement or the site total is missing.
  out <- mod[!is.na(mod$diff_abund) & mod$diff_abund != 0, , drop = FALSE]
  if (verbose)
    cat("Non-zero paired differences:", nrow(out),
        " (neg:", sum(out$diff_abund < 0),
        " pos:", sum(out$diff_abund > 0), ")\n")

  out$Predominant_simple <- relevel(factor(out$Predominant_simple), ref = "Secondary")
  out$SS   <- factor(out$SS)
  out$SSB  <- factor(out$SSB)
  out$SSBS <- factor(out$SSBS)
  out
}

# Parameter names of biome x land-use cells with no data, whose
# coefficients are prior-only and must be flagged rather than
# interpreted. Derived from the fitted data, not hardcoded: the
# diagnostics scripts used to name a specific Temperate Open cell,
# which silently became dead code the moment that biome was dropped.
# Returns character(0) when every cell is populated, as it is now.
#
# Cells involving a reference level are skipped: they produce no
# interaction parameter, so there is nothing to flag. ref_lu differs
# between the record-level models (Primary vegetation) and the paired
# difference (Secondary, since primary is absorbed into the response).
empty_cell_params <- function(dat, ref_biome = BIOME_REF, ref_lu = LU_REF) {
  ct  <- table(dat$Biome4, dat$Predominant_simple)
  idx <- which(ct == 0, arr.ind = TRUE)
  if (!nrow(idx)) return(character(0))
  strip <- function(x) gsub("[^A-Za-z0-9]", "", x)
  b <- rownames(ct)[idx[, 1]]
  l <- colnames(ct)[idx[, 2]]
  keep <- b != ref_biome & l != ref_lu
  if (!any(keep)) return(character(0))
  paste0("Biome4", strip(b[keep]), ":Predominant_simple", strip(l[keep]))
}

# ------------------------------------------------------------
# Phylogenetic tree matching
# ------------------------------------------------------------
# Genus reclassifications where the same biological species sits
# under a different name in BBtree2 (Jetz-era taxonomy).
TREE_SYNONYMS <- c(
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

# Three-stage match: jetz_sp, then underscored Best_guess_binomial,
# then the synonym table. Returns a species -> tip lookup with
# ambiguous tips (claimed by >1 species) removed.
match_to_tree <- function(df, tree, verbose = TRUE) {
  stopifnot(all(TREE_SYNONYMS %in% tree$tip.label))
  sp <- unique(df[, c("Best_guess_binomial", "jetz_sp")])
  sp$binom_u <- gsub(" ", "_", sp$Best_guess_binomial)
  sp$phylo <- ifelse(sp$jetz_sp %in% tree$tip.label, sp$jetz_sp, NA_character_)
  sp$phylo <- ifelse(!is.na(sp$phylo), sp$phylo,
                     ifelse(sp$binom_u %in% tree$tip.label, sp$binom_u, NA_character_))
  for (jn in names(TREE_SYNONYMS))
    sp$phylo[is.na(sp$phylo) & sp$jetz_sp == jn] <- TREE_SYNONYMS[[jn]]
  sp <- sp[!is.na(sp$phylo), ]
  dupes <- unique(sp$phylo[duplicated(sp$phylo)])
  if (length(dupes)) {
    if (verbose) cat("Dropping", length(dupes), "ambiguous tips\n")
    sp <- sp[!(sp$phylo %in% dupes), ]
  }
  if (verbose) cat("Tree matching:", nrow(sp), "species matched to tips\n")
  sp[, c("Best_guess_binomial", "phylo")]
}

# ------------------------------------------------------------
# Priors
# ------------------------------------------------------------
# Lazy (a function, not an object) so 00_data_preparation.R can
# source this file without brms installed.
#
#   "lognormal"   meancol / malecol / dichro ratio, and abundance
#   "gaussian"    dichrodiff at the record level
#   "paired_diff" A2d/B2d, tight priors for the +/-0.06 difference scale
#   "phylo_*"     species-level A3/B3, wider fixed-effect priors
model_priors <- function(type = c("lognormal", "gaussian", "paired_diff",
                                  "phylo_lognormal", "phylo_gaussian")) {
  type <- match.arg(type)
  sp <- brms::set_prior
  switch(type,
    lognormal = c(sp("normal(0, 2)",   class = "Intercept"),
                  sp("normal(0, 0.5)", class = "b"),
                  sp("exponential(2)", class = "sd"),
                  sp("exponential(2)", class = "sigma")),
    gaussian  = c(sp("normal(0, 50)",  class = "Intercept"),
                  sp("normal(0, 20)",  class = "b"),
                  sp("exponential(2)", class = "sd"),
                  sp("exponential(2)", class = "sigma")),
    paired_diff = c(sp("normal(0, 0.5)", class = "Intercept"),
                    sp("normal(0, 0.2)", class = "b"),
                    sp("exponential(5)", class = "sd"),
                    sp("exponential(5)", class = "sigma")),
    # No intercept in A3/B3 (the compositional predictors span it), so
    # b coefficients are absolute levels and take a wider prior.
    # NB these are the values the A3 code has always used. ANALYSIS_NOTES
    # previously recorded normal(0,100)/exponential(2) for the Gaussian
    # case, which never matched the script; the code is authoritative.
    phylo_lognormal = c(sp("normal(0, 5)",      class = "b"),
                        sp("exponential(2)",    class = "sd"),
                        sp("exponential(2)",    class = "sigma")),
    phylo_gaussian  = c(sp("normal(0, 50)",     class = "b"),
                        sp("exponential(0.05)", class = "sd"),
                        sp("exponential(0.05)", class = "sigma")))
}

# ------------------------------------------------------------
# Response metrics
# ------------------------------------------------------------
# Cooney UVS (A family) and Dale (B family). RESP_FAMILY maps each
# metric to its brms family; the three positive metrics are
# lognormal, the male-minus-female difference is Gaussian.
RESP_COONEY <- c("meancolcooney", "dichrocooney", "malecolcooney", "dichrodiff")
RESP_DALE   <- c("meancoldale",   "dichrodale",   "malecoldale",   "dichrodiffdale")

RESP_FAMILY <- c(meancolcooney = "lognormal", dichrocooney   = "lognormal",
                 malecolcooney = "lognormal", dichrodiff     = "gaussian",
                 meancoldale   = "lognormal", dichrodale     = "lognormal",
                 malecoldale   = "lognormal", dichrodiffdale = "gaussian")

# Metrics carried into the abundance models (A2c/A2d, B2c/B2d).
COLOUR_COONEY <- c("z_malecolcooney", "z_dichrodiff")
COLOUR_DALE   <- c("z_malecoldale",   "z_dichrodiffdale")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

message("00_config.R loaded | Biome4 = 3 levels (", BIOME_DROP,
        " dropped) | threads/chain: ", nthreads)
