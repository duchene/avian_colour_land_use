# ============================================================
# C2_sex_averaged.R
# Sensitivity of the colourfulness results to using the mean of male and
# female colourfulness instead of male colourfulness.
#
# The Methods state that male colourfulness is strongly linked with the
# mean colourfulness of the sexes, so that both give nearly identical
# results. This script shows it (co-author request: supplementary
# information rather than "results not shown").
#
# Sex-averaged colourfulness = (male LociUVS + female LociUVS) / 2
# (Cooney et al. 2022). Each reported colourfulness model is refitted
# with it, keeping the formula, priors, data rows, sampler and seed of the
# reported fit (brms update with new data):
#   A1b  species level, lognormal, fits/A1b_fit_malecolcooney.RData
#   A2d  abundance change, colourfulness as a z-scored predictor,
#        fits/A2d_fit_z_malecolcooney.RData
#   A3b  assemblage level with covariates,
#        fits/A3b_fit_malecolcooney_covariate.RData
#
# Usage (from the repo root, on a compute node):
#   C2_MODEL=A1b|A2d|A3b Rscript scripts/cooney/C2_sex_averaged.R   fit one
#   C2_MODEL=summary     Rscript scripts/cooney/C2_sex_averaged.R   compare
# Outputs: fits/C2_<model>_meancol.RData,
#          results/cooney/C2_sex_averaged_comparison.csv
# ============================================================

source("scripts/cooney/posthoc_helpers.R")   # attaches rstan and brms
suppressPackageStartupMessages({ library(cmdstanr); library(tidyr) })
options(brms.backend = "cmdstanr")

MODEL <- Sys.getenv("C2_MODEL", "summary")
nthreads <- if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK")))
  max(1, floor(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) / 4)) else 1
SRC <- c(A1b = "fits/A1b_fit_malecolcooney.RData",
         A2d = "fits/A2d_fit_z_malecolcooney.RData",
         A3b = "fits/A3b_fit_malecolcooney_covariate.RData")
CTRL <- list(A1b = list(adapt_delta = 0.95, max_treedepth = 12),
             A2d = list(adapt_delta = 0.95, max_treedepth = 12),
             A3b = list(adapt_delta = 0.90, max_treedepth = 10))
out_fit <- function(m) sprintf("fits/C2_%s_meancol.RData", m)

co <- read.csv("data/external/cooney2022_passerines.csv", stringsAsFactors = FALSE)
co$meancol <- (co$LociUVS_male + co$LociUVS_female) / 2

# The reported fit's own data, with male colourfulness replaced by the
# sex average of the same species. Rows are those of the reported fit.
swap_data <- function(m, fit) {
  d <- fit$data
  if (m == "A1b") {
    # A1b_01_setup.R already stores the sex average (meancolcooney) per
    # species; the male value must reproduce the fitted one
    e <- new.env(); load("fits/A1b_model_setup.RData", envir = e)
    i <- match(as.character(d$phylo), as.character(e$spdat$phylo))
    stopifnot(!anyNA(i), isTRUE(all.equal(e$spdat$malecolcooney[i], d$malecolcooney)))
    mc <- e$spdat$meancolcooney[i]
    stopifnot(!anyNA(mc))
    d$malecolcooney <- mc
  } else {
    e <- new.env(); load(sprintf("fits/%s_model_setup.RData", m), envir = e)
    full <- if (m == "A2d") e$dat_final else e$dat
    full <- full[!is.na(full$malecolcooney), ]
    stopifnot(nrow(full) == nrow(d))
    mc <- (full$LociUVS_male_cooney + full$LociUVS_female_cooney) / 2
    stopifnot(!anyNA(mc), isTRUE(all.equal(full$malecolcooney,
              if (m == "A2d") full$LociUVS_male_cooney else d$malecolcooney)))
    if (m == "A2d") {
      # z-score over the same rows the reported predictor was z-scored on
      stopifnot(isTRUE(all.equal(d$z_malecolcooney, full$z_malecolcooney)))
      d$z_malecolcooney <- (mc - mean(mc)) / sd(mc)
    } else d$malecolcooney <- mc
  }
  d
}

if (MODEL %in% names(SRC)) {
  load(SRC[[MODEL]])
  nd <- swap_data(MODEL, fit)
  cat("C2", MODEL, ": n =", nrow(nd), "| cor(male, sex-averaged) over rows =",
      round(cor(fit$data[[if (MODEL == "A2d") "z_malecolcooney" else "malecolcooney"]],
                nd[[if (MODEL == "A2d") "z_malecolcooney" else "malecolcooney"]]), 3), "\n")
  fit <- update(fit, newdata = nd, recompile = TRUE, cores = 4,
                threads = threading(nthreads), control = CTRL[[MODEL]],
                seed = 42, refresh = 200)
  save(fit, file = out_fit(MODEL))
  cat("saved", out_fit(MODEL), "\n")
  quit(save = "no")
}

# ---------- summary: the reported quantities under both metrics ----------
# A1b: land use minus primary vegetation (coefficient differences).
q_A1b <- function(fit) {
  d <- posterior::as_draws_matrix(fit)
  lus <- c(Secondary = "prop_Secondary", Plantation = "prop_Plantation_forest",
           Cropland = "prop_Cropland", Pasture = "prop_Pasture")
  m <- sapply(lus, function(l) d[, paste0("b_", l)] - d[, "b_prop_Primary_vegetation"])
  colnames(m) <- paste(names(lus), "- Primary vegetation")
  columnwise(m, link_sd(fit), TRUE, quantity = "Land-use association contrast",
             held_fixed = NA_character_)
}
# A2d: colour slope within each land use and biome, and their contrasts.
q_A2d <- function(fit) {
  g <- make_grid(Predominant_simple = LU4, Biome4 = BIOME3); g$z_malecolcooney <- 0
  sdl <- link_sd(fit)
  bind_rows(lapply(c("Predominant_simple", "Biome4"), function(t) {
    s <- slope_draws(fit, g, t, "z_malecolcooney")
    bind_rows(columnwise(s, sdl, FALSE, quantity = "Colour slope", held_fixed = NA_character_),
              pairwise(s, sdl, FALSE, quantity = "Colour slope contrast", held_fixed = NA_character_))
  }))
}
# A3b: land-use contrasts within each biome, and the body-mass slope.
q_A3b <- function(fit) {
  g <- make_grid(Predominant_simple = LU5, Biome4 = BIOME3,
                 Trophic.Niche = TROPHIC6, z_logMass = 0)
  sdl <- link_sd(fit)
  bind_rows(lapply(BIOME3, function(b) {
    pairwise(emm_draws(fit, g[g$Biome4 == b, ], "Predominant_simple"), sdl, TRUE,
             quantity = "Land-use contrast", held_fixed = paste("Biome =", b))
  }), columnwise(slope_draws(fit, transform(g, overall = "All"), "overall", "z_logMass"),
                 sdl, TRUE, quantity = "Body-mass slope", held_fixed = NA_character_))
}
QF <- list(A1b = q_A1b, A2d = q_A2d, A3b = q_A3b)

rows <- list()
for (m in names(SRC)) {
  if (!file.exists(out_fit(m))) { message("missing ", out_fit(m), ", skipped"); next }
  for (v in c(male = SRC[[m]], mean = out_fit(m))) {
    load(v)
    r <- QF[[m]](fit)
    r$model <- m; r$metric <- if (v == SRC[[m]]) "male" else "sex_averaged"
    rows[[paste(m, v)]] <- r
    rm(fit); gc(verbose = FALSE)
  }
}
tab <- bind_rows(rows) %>%
  mutate(held_fixed = relabel_levels(held_fixed), contrast = relabel_levels(contrast)) %>%
  select(model, quantity, held_fixed, contrast, metric, estimate, lo, hi, credible) %>%
  pivot_wider(names_from = metric, values_from = c(estimate, lo, hi, credible)) %>%
  mutate(same_conclusion = credible_male == credible_sex_averaged &
           (!credible_male | sign(estimate_male) == sign(estimate_sex_averaged)))
write.csv(tab, "results/cooney/C2_sex_averaged_comparison.csv", row.names = FALSE)

cat("\nCorrelation of male and sex-averaged colourfulness across",
    sum(complete.cases(co[, c("LociUVS_male", "meancol")])), "passerines (log scale):",
    round(cor(log(co$LociUVS_male), log(co$meancol), use = "complete.obs"), 3), "\n")
print(tab %>% group_by(model) %>%
        summarise(n = n(), same = sum(same_conclusion),
                  credible_male = sum(credible_male), credible_avg = sum(credible_sex_averaged),
                  cor_estimates = cor(estimate_male, estimate_sex_averaged)))
print(as.data.frame(tab %>% filter(!same_conclusion)), digits = 3)
cat("Saved results/cooney/C2_sex_averaged_comparison.csv\n")
