# ============================================================
# A3b_colours.R
# Assemblage-level colour (A3b) refitted once per Delhey et al. (2023)
# plumage colour: does the share of the male plumage in a given colour,
# across the species recorded at a site, differ between land-use types
# within each biome?
#
# Same records, fixed effects (biome x land use, trophic niche, log body
# mass), random effects (study, block, site) and sampler as the reported
# A3b colourfulness model (fits/A3b_fit_malecolcooney_covariate.RData);
# records of species without Delhey data are dropped (species matched on
# their BirdTree name). No phylogeny, as for
# A3b (colour is constant within species; see ANALYSIS_NOTES.md).
#
# Response: each colour proportion on the logit scale after the Smithson &
# Verkuilen (2006) squeeze, as in the MR-PMM (A4_colour_mrpmm.R), then
# z-standardised across species, so effects are in species-level SD.
# Gaussian family; priors N(0, 2) intercept, N(0, 1) coefficients,
# exponential(2) standard deviations.
#
# Caveat: each species is repeated at every site where it was recorded,
# so the effective sample size is nearer the number of species than the
# number of records, and intervals are likely too narrow (as for A3b).
#
# Usage (repo root, compute node):
#   A3C_COL=<colour> Rscript scripts/cooney/A3b_colours.R   fit one colour
#   A3C_COL=summary  Rscript scripts/cooney/A3b_colours.R   summarise
#   A3C_FAM=bernoulli with A3C_COL=<colour>: presence check (see below)
#
# Presence check: for colours absent from many species (rufous, red,
# green, blue), the logit-Gaussian model cannot reproduce the excess of
# zeros (posterior predictive skewness). Those colours are also fitted as
# presence (any of the colour in the male plumage) with a Bernoulli model,
# same fixed and random effects; its land-use contrasts (logit scale) are
# compared in direction with the Gaussian ones.
# Fits:    fits/A3b_fit_col_<colour>.RData, fits/A3b_fit_colpres_<colour>.RData
# Outputs: results/cooney/A3b_colours_convergence.csv
#          results/cooney/A3b_colours_contrasts.csv   (land-use pairs within biome)
#          results/cooney/A3b_colours_ppc.csv
#          results/cooney/A3b_colours_presence_contrasts.csv (presence check)
# ============================================================

source("scripts/cooney/posthoc_helpers.R")   # attaches rstan and brms
source("scripts/cooney/delhey_colours.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
options(brms.backend = "cmdstanr")

COLS <- c("black", "white", "grey", "brown", "rufous", "yellow", "red", "green", "blue")
MODE <- Sys.getenv("A3C_COL", "summary")
FAM  <- Sys.getenv("A3C_FAM", "gaussian"); stopifnot(FAM %in% c("gaussian", "bernoulli"))
PRES <- c("rufous", "red", "green", "blue")
out_fit <- function(cc, fam = "gaussian")
  sprintf(if (fam == "bernoulli") "fits/A3b_fit_colpres_%s.RData" else "fits/A3b_fit_col_%s.RData", cc)

build_data <- function() {
  source("scripts/00_config.R")
  load("fits/A3b_model_setup.RData")            # dat
  e <- new.env(); load("fits/A3b_fit_malecolcooney_covariate.RData", envir = e)
  keep <- complete.cases(dat[, c("malecolcooney", "Biome4", "Predominant_simple",
                                 "Trophic.Niche", "z_logMass", "SS", "SSB", "SSBS")])
  d <- dat[keep, ]
  stopifnot(nrow(d) == nobs(e$fit), isTRUE(all.equal(d$malecolcooney, e$fit$data$malecolcooney)))
  # Delhey names follow the BirdTree (Jetz et al. 2012) taxonomy, like
  # jetz_sp; A3b is not restricted to species on the tree, so match on
  # jetz_sp, then the tree label, then the binomial. Delhey names claimed
  # by two PREDICTS species are ambiguous and dropped.
  dl <- delhey_male()$phylo
  sp <- unique(d[, c("Best_guess_binomial", "jetz_sp")])
  tr <- match_to_tree(d, ape::read.tree("data/BBtree2.tre"), verbose = FALSE, phylo = "claramunt")
  sp$tree <- tr$phylo[match(sp$Best_guess_binomial, tr$Best_guess_binomial)]
  sp$bu <- gsub(" ", "_", sp$Best_guess_binomial)
  sp$key <- ifelse(sp$jetz_sp %in% dl, sp$jetz_sp, ifelse(sp$tree %in% dl, sp$tree,
                   ifelse(sp$bu %in% dl, sp$bu, NA_character_)))
  sp$key[sp$key %in% sp$key[duplicated(sp$key) & !is.na(sp$key)]] <- NA
  stopifnot(!anyDuplicated(sp$Best_guess_binomial))
  d$phylo <- sp$key[match(d$Best_guess_binomial, sp$Best_guess_binomial)]
  d <- add_delhey_colours(d, d$phylo)
  d <- d[!is.na(d$black), ]
  # species-level transform: squeeze, logit, z across species
  spv <- d[!duplicated(d$phylo), c("phylo", COLS)]
  nsq <- nrow(spv)
  for (cc in COLS) {
    y <- qlogis((spv[[cc]] * (nsq - 1) + 0.5) / nsq)
    d[[paste0("y_", cc)]] <- (y[match(d$phylo, spv$phylo)] - mean(y)) / sd(y)
  }
  list(d = d, n_species = nsq)
}

if (MODE != "summary") {
  stopifnot(MODE %in% COLS)
  if (file.exists(out_fit(MODE, FAM)) && !nzchar(Sys.getenv("FORCE"))) quit(save = "no")
  b <- build_data(); d <- b$d
  nthreads <- if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK")))
    max(1, floor(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) / 4)) else 1
  yv <- paste0("y_", MODE)
  if (FAM == "bernoulli") { yv <- paste0("p_", MODE); d[[yv]] <- as.integer(d[[MODE]] > 0) }
  form <- bf(as.formula(paste0(yv, " ~ Biome4 * Predominant_simple + Trophic.Niche + z_logMass",
                               " + (1|SS) + (1|SSB) + (1|SSBS)")))
  pri <- c(set_prior("normal(0, 2)", class = "Intercept"), set_prior("normal(0, 1)", class = "b"),
           set_prior("exponential(2)", class = "sd"), set_prior("exponential(2)", class = "sigma"))
  cat("A3b colour", MODE, "| n =", nrow(d), "records |", b$n_species, "species\n")
  if (FAM == "bernoulli") pri <- pri[pri$class != "sigma", ]
  fit <- brm(form, data = d, family = if (FAM == "bernoulli") bernoulli() else gaussian(), prior = pri,
             chains = 4, iter = 3000, warmup = 1500, cores = 4, threads = threading(nthreads),
             control = list(adapt_delta = 0.9, max_treedepth = 10), seed = 1, refresh = 200)
  save(fit, file = out_fit(MODE, FAM))
  cat("saved", out_fit(MODE, FAM), "\n")
  quit(save = "no")
}

# ---------------- summary ----------------
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
g <- make_grid(Predominant_simple = LU5, Biome4 = BIOME3, Trophic.Niche = TROPHIC6, z_logMass = 0)
conv <- list(); rows <- list(); ppc <- list()
for (cc in COLS) {
  if (!file.exists(out_fit(cc))) { message("missing ", out_fit(cc), ", skipped"); next }
  load(out_fit(cc))
  sm <- summarise_draws(as_draws_array(fit, variable = "^(b_|sd_|sigma)", regex = TRUE),
                        "rhat", "ess_bulk", "ess_tail")
  np <- nuts_params(fit)
  conv[[cc]] <- data.frame(colour = COLOUR_LABEL[[cc]], n_records = nobs(fit),
    max_rhat = max(sm$rhat), min_ess_bulk = min(sm$ess_bulk), min_ess_tail = min(sm$ess_tail),
    n_divergent = sum(np$Value[np$Parameter == "divergent__"]))
  # land-use contrasts within each biome, in species-level SD (no rescaling)
  for (bi in BIOME3)
    rows[[paste(cc, bi)]] <- pairwise(emm_draws(fit, g[g$Biome4 == bi, ], "Predominant_simple"),
                                      1, FALSE, colour = COLOUR_LABEL[[cc]], biome = bi)
  # posterior predictive check: skew and tails of the record-level response
  obs <- fit$data[[1]]
  yrep <- posterior_predict(fit, ndraws = 200)
  st <- function(y) c(sd = sd(y), skew = skew(y), q10 = unname(quantile(y, 0.1)),
                      q90 = unname(quantile(y, 0.9)))
  rs <- t(apply(yrep, 1, st)); os <- st(obs)
  ppc[[cc]] <- data.frame(colour = COLOUR_LABEL[[cc]], stat = names(os), observed = os,
                          replicated_median = apply(rs, 2, median),
                          ppp = colMeans(sweep(rs, 2, os, ">=")), row.names = NULL)
  rm(fit, yrep); gc(verbose = FALSE)
}
# presence check: land-use contrasts on the logit scale of presence
prow <- list()
for (cc in PRES) {
  f <- out_fit(cc, "bernoulli")
  if (!file.exists(f)) { message("missing ", f, ", skipped"); next }
  load(f)
  sm <- summarise_draws(as_draws_array(fit, variable = "^(b_|sd_)", regex = TRUE), "rhat", "ess_bulk")
  conv[[paste0(cc, "_presence")]] <- data.frame(colour = paste(COLOUR_LABEL[[cc]], "(presence)"),
    n_records = nobs(fit), max_rhat = max(sm$rhat), min_ess_bulk = min(sm$ess_bulk),
    min_ess_tail = NA_real_, n_divergent = sum(nuts_params(fit)$Value[nuts_params(fit)$Parameter == "divergent__"]))
  for (bi in BIOME3)
    prow[[paste(cc, bi)]] <- pairwise(emm_draws(fit, g[g$Biome4 == bi, ], "Predominant_simple"),
                                      1, FALSE, colour = COLOUR_LABEL[[cc]], biome = bi)
  rm(fit); gc(verbose = FALSE)
}
conv <- do.call(rbind, conv); ppc <- do.call(rbind, ppc)
tab <- bind_rows(rows) %>%
  mutate(biome = relabel_levels(biome), contrast = relabel_levels(contrast)) %>%
  select(colour, biome, contrast, estimate, lo, hi, pd, credible)
write.csv(conv, "results/cooney/A3b_colours_convergence.csv", row.names = FALSE)
write.csv(tab, "results/cooney/A3b_colours_contrasts.csv", row.names = FALSE)
write.csv(ppc, "results/cooney/A3b_colours_ppc.csv", row.names = FALSE)
print(conv, row.names = FALSE, digits = 3)
print(ppc, row.names = FALSE, digits = 3)
cat("\nCredible land use minus primary vegetation contrasts:\n")
print(as.data.frame(tab %>% filter(credible, grepl("- Primary vegetation$", contrast))), digits = 3)
if (length(prow)) {
  pt <- bind_rows(prow) %>%
    mutate(biome = relabel_levels(biome), contrast = relabel_levels(contrast)) %>%
    select(colour, biome, contrast, estimate, lo, hi, pd, credible)
  cmpd <- merge(tab, pt, by = c("colour", "biome", "contrast"), suffixes = c("_gaussian", "_presence"))
  cmpd$same_direction <- sign(cmpd$estimate_gaussian) == sign(cmpd$estimate_presence)
  write.csv(cmpd, "results/cooney/A3b_colours_presence_contrasts.csv", row.names = FALSE)
  v <- cmpd[grepl("- Primary vegetation$", cmpd$contrast), ]
  cat("\nPresence check, land use minus primary vegetation:", nrow(v), "contrasts;",
      sum(v$same_direction), "same direction;", sum(v$credible_gaussian), "credible (Gaussian),",
      sum(v$credible_presence), "credible (presence),",
      sum(v$credible_gaussian & v$credible_presence & v$same_direction), "credible in both, same sign\n")
  print(v[v$credible_gaussian | v$credible_presence, c("colour", "biome", "contrast", "estimate_gaussian",
          "credible_gaussian", "estimate_presence", "credible_presence")], digits = 3, row.names = FALSE)
}
cat("Saved results/cooney/A3b_colours_{convergence,contrasts,ppc}.csv\n")
