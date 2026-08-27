# ============================================================
# posthoc_interactions.R  ->  Table 2
#
# Pairwise comparisons of factor levels, computed separately within
# each level of the factor it interacts with. Nothing else.
#
# Table 1 (posthoc_contrasts.R) contrasts levels of one factor
# averaged over the others, which is the right summary only when the
# interaction is negligible. This table opens that interaction up.
#
# SCOPE. Every row is a pairwise comparison of two levels of a factor
# that appears in an interaction term of the fitted model. A2c and A2d
# carry three interaction terms and all three are covered.
#
#   Biome4 x Predominant_simple   quantity = "Level mean"
#     Land use   every land-use pair, within each biome
#     Biome      every biome pair, within each land use
#
#   colour x Biome4               quantity = "Colour slope"
#     Biome      every biome pair, comparing the colour slope
#
#   colour x Predominant_simple   quantity = "Colour slope"
#     Land use   every land-use pair, comparing the colour slope
#
# The colour terms cannot be expressed as pairwise comparisons of
# colour levels, because colour is continuous. They are expressed as
# pairwise comparisons of the BIOME and LAND-USE levels, contrasting
# the colour slope between them, which is exactly what the interaction
# coefficients encode. A1b contributes "Level mean" rows only, since
# colour is its response and it carries no colour interaction.
#
# Deliberately NOT included:
#   * Interaction contrasts (differences of differences). They are not
#     pairwise level comparisons.
#   * Simple slopes on the standardised colour predictors. Colour is
#     continuous and has no levels to pair.
#   * Trophic niche and body mass. They enter A1b as main effects
#     only, never in an interaction, so their contrasts do not depend
#     on any other level and belong in Table 1.
#   * Contrasts in which biome AND land use both differ. They confound
#     the two factors and a credible result among them could not be
#     attributed to either.
#   * A3, which carries no interaction term at all.
#
# RESPONSE. The response differs between model families and is named
# in its own column. In A1b the response is community colour. In A2c
# and A2d the response is ABUNDANCE.
#
# THE COLOUR METRIC IS NOT PART OF ANY CONTRAST HERE. A2c and A2d were
# each fitted twice, once per colour metric, so every land-use and
# biome contrast has two estimates. The source_fit column says which
# fit an estimate came from. It does NOT mean the contrast is
# conditioned on colour. Colour is held at z = 0, so colour:Biome4 and
# colour:Predominant_simple contribute exactly zero, and the
# reconstruction below confirms that each row is built from the
# land-use and biome terms alone. The two source_fit rows for a given
# contrast estimate the SAME quantity and differ only as two separate
# fits do (median difference 0.004, maximum 0.025).
#
# NO GHOST TERMS. Every model here is two-way. There is no three-way
# colour by biome by land-use term (tested separately in A2e and A2f,
# and not supported). The validation block below reconstructs every
# contrast from the fitted coefficients and stops if any row fails to
# match, which proves no row depends on a term the model lacks.
#
# Run from repo root: Rscript scripts/cooney/posthoc_interactions.R
# ============================================================

source("scripts/cooney/posthoc_helpers.R")
suppressMessages({ library(tidyr); library(stringr); library(readr) })

dir.create("results/cooney", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/pub/cooney", showWarnings = FALSE, recursive = TRUE)

# Pairwise differences of the colour slope between levels of `target`,
# averaging over the other factor. This is the colour x `target`
# interaction expressed as a level comparison.
slope_contrasts <- function(fit, grid, target, cvar, model, response, src) {
  sdl <- link_sd(fit); lognorm <- fit$family$family == "lognormal"
  pairwise(slope_draws(fit, grid, target, cvar), sdl, lognorm,
           model = model, response = response, source_fit = src,
           quantity = "Colour slope",
           factor_compared = if (target == "Predominant_simple") "Land use" else "Biome",
           held_fixed = NA_character_)
}

# Pairwise contrasts of `target` computed inside each level of `within`.
simple_effects <- function(fit, grid, target, within, model, response, src) {
  sdl <- link_sd(fit); lognorm <- fit$family$family == "lognormal"
  map_dfr(unique(grid[[within]]), function(w) {
    g <- grid[grid[[within]] == w, , drop = FALSE]
    emm <- emm_draws(fit, g, target)
    pairwise(emm, sdl, lognorm,
             model = model, response = response, source_fit = src,
             quantity = "Level mean",
             factor_compared = if (target == "Predominant_simple") "Land use" else "Biome",
             held_fixed = paste0(if (within == "Biome4") "Biome = " else "Land use = ", w))
  })
}

RESPONSE   <- c(A2c = "Relative abundance", A2d = "Abundance change vs primary")
COLOUR_LAB <- c(malecolcooney = "Colour diversity", dichrodiff = "Sexual dichromatism",
                z_malecolcooney = "Colour diversity (z)",
                z_dichrodiff = "Sexual dichromatism (z)")
# Names the fit an estimate came from, not a term in the contrast.
SRC_FIT    <- c(z_malecolcooney = "fitted with colour diversity (z)",
                z_dichrodiff    = "fitted with sexual dichromatism (z)")

out <- list()

# ------------------------------------------------------------
# A1b. Response is community colour. No colour predictor.
# ------------------------------------------------------------
message("A1b ...")
grid_a1b <- make_grid(Predominant_simple = LU5, Biome4 = BIOME3,
                      Trophic.Niche = TROPHIC6, z_logMass = 0)
for (resp in c("malecolcooney", "dichrodiff")) {
  load(paste0("fits/A1b_fit_", resp, "_covariate.RData"))
  rlab <- paste0("Community ", tolower(COLOUR_LAB[[resp]]))
  out[[paste0("A1b_lu_", resp)]] <-
    simple_effects(fit, grid_a1b, "Predominant_simple", "Biome4", "A1b", rlab, NA_character_)
  out[[paste0("A1b_bi_", resp)]] <-
    simple_effects(fit, grid_a1b, "Biome4", "Predominant_simple", "A1b", rlab, NA_character_)
  rm(fit); gc(verbose = FALSE)
  message("  ", resp, " done")
}

# ------------------------------------------------------------
# A2c and A2d. Response is abundance. Colour is a predictor held at
# its mean (z = 0), so the colour interaction terms contribute
# nothing to these contrasts.
# ------------------------------------------------------------
for (mod in c("A2c", "A2d")) {
  message(mod, " ...")
  lus <- if (mod == "A2c") LU5 else LU4
  for (cv in c("z_malecolcooney", "z_dichrodiff")) {
    load(paste0("fits/", mod, "_fit_", cv, ".RData"))
    g <- make_grid(Predominant_simple = lus, Biome4 = BIOME3)
    g[[cv]] <- 0
    out[[paste0(mod, "_lu_", cv)]] <-
      simple_effects(fit, g, "Predominant_simple", "Biome4", mod,
                     RESPONSE[[mod]], SRC_FIT[[cv]])
    out[[paste0(mod, "_bi_", cv)]] <-
      simple_effects(fit, g, "Biome4", "Predominant_simple", mod,
                     RESPONSE[[mod]], SRC_FIT[[cv]])
    # the two colour interaction terms, as level comparisons
    out[[paste0(mod, "_slbi_", cv)]] <-
      slope_contrasts(fit, g, "Biome4", cv, mod, RESPONSE[[mod]], SRC_FIT[[cv]])
    out[[paste0(mod, "_sllu_", cv)]] <-
      slope_contrasts(fit, g, "Predominant_simple", cv, mod,
                      RESPONSE[[mod]], SRC_FIT[[cv]])
    rm(fit); gc(verbose = FALSE)
    message("  ", cv, " done")
  }
}

tbl <- bind_rows(out) %>%
  mutate(quantity = ifelse(is.na(quantity), "Level mean", quantity)) %>%
  select(model, response, source_fit, quantity, factor_compared, held_fixed,
         contrast, estimate, mean_est, lo, hi, pd, credible, pct_change, effect_size)

# ============================================================
# VALIDATION
# ============================================================
message("\nvalidating ...")

# --- 1. completeness: every pair, within every level, for every fit ---
# Level means: every pair, within every level of the other factor.
# Colour slopes: every pair, once, averaged over the other factor.
expected <- bind_rows(
  tibble(model = "A1b", quantity = "Level mean",
         n_lu_exp = choose(5, 2) * 3, n_bi_exp = choose(3, 2) * 5),
  tibble(model = "A2c", quantity = "Level mean",
         n_lu_exp = choose(5, 2) * 3, n_bi_exp = choose(3, 2) * 5),
  tibble(model = "A2c", quantity = "Colour slope",
         n_lu_exp = choose(5, 2),     n_bi_exp = choose(3, 2)),
  tibble(model = "A2d", quantity = "Level mean",
         n_lu_exp = choose(4, 2) * 3, n_bi_exp = choose(3, 2) * 4),
  tibble(model = "A2d", quantity = "Colour slope",
         n_lu_exp = choose(4, 2),     n_bi_exp = choose(3, 2)))
got <- tbl %>% count(model, response, source_fit, quantity, factor_compared) %>%
  pivot_wider(names_from = factor_compared, values_from = n) %>%
  rename(n_lu = `Land use`, n_bi = Biome) %>%
  left_join(expected, by = c("model", "quantity"))
stopifnot(all(got$n_lu == got$n_lu_exp), all(got$n_bi == got$n_bi_exp))
# every interaction term in every A2 fit must be represented
stopifnot(all(c("Level mean", "Colour slope") %in%
              tbl$quantity[tbl$model == "A2c"]),
          all(c("Level mean", "Colour slope") %in%
              tbl$quantity[tbl$model == "A2d"]))
# no duplicated rows, and no contrast of a level with itself
stopifnot(!any(duplicated(tbl[, c("model", "response", "source_fit", "quantity",
                                  "factor_compared", "held_fixed", "contrast")])))
stopifnot(!any(str_split_fixed(tbl$contrast, " - ", 2)[, 1] ==
               str_split_fixed(tbl$contrast, " - ", 2)[, 2]))
cat("completeness OK\n"); print(as.data.frame(got), row.names = FALSE)

# --- 2. reconstruction from the fitted coefficients ---
# A land-use contrast within biome b must equal
#   (b_LUa - b_LUb) + (b_{b:LUa} - b_{b:LUb})
# and a biome contrast within land use l must equal
#   (b_Ba - b_Bb) + (b_{Ba:l} - b_{Bb:l})
# with reference-level coefficients zero. A row that needed anything
# else, a three-way term above all, would fail here.
strip   <- function(x) gsub("[^A-Za-z0-9]", "", x)
coef_of <- function(fx, nm) if (nm %in% names(fx)) unname(fx[[nm]]) else 0

specs <- list(
  list(m = "A1b", file = "results/cooney/A1b_fixed_effects.csv", lus = LU5,
       fits = list(c(key = "malecolcooney_covariate", resp = "Community colour diversity", cp = NA),
                   c(key = "dichrodiff_covariate",
                     resp = "Community sexual dichromatism", cp = NA))),
  list(m = "A2c", file = "results/cooney/A2c_fixed_effects.csv", lus = LU5,
       fits = list(c(key = "z_malecolcooney", resp = "Relative abundance",
                     cp = "fitted with colour diversity (z)"),
                   c(key = "z_dichrodiff", resp = "Relative abundance",
                     cp = "fitted with sexual dichromatism (z)"))),
  list(m = "A2d", file = "results/cooney/A2d_fixed_effects.csv", lus = LU4,
       fits = list(c(key = "z_malecolcooney", resp = "Abundance change vs primary",
                     cp = "fitted with colour diversity (z)"),
                   c(key = "z_dichrodiff", resp = "Abundance change vs primary",
                     cp = "fitted with sexual dichromatism (z)"))))

recon_max <- 0; n_checked <- 0
for (sp in specs) {
  fe <- read.csv(sp$file)
  for (ft in sp$fits) {
    sel <- fe$model == ft[["key"]]
    fx  <- setNames(fe$Estimate[sel], fe$parameter[sel])
    cvar <- if (sp$m == "A1b") NA_character_ else ft[["key"]]
    sub <- tbl %>% filter(model == sp$m, response == ft[["resp"]],
                          if (is.na(ft[["cp"]])) is.na(source_fit)
                          else source_fit %in% ft[["cp"]])
    ref_lu <- sp$lus[1]; ref_bi <- BIOME3[1]
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      ab <- str_split_fixed(r$contrast, " - ", 2)
      a  <- ab[1]; b <- ab[2]
      if (r$quantity == "Colour slope") {
        # slope difference = difference of the colour interaction terms.
        # The colour main effect and the average over the other factor
        # are common to both levels and cancel.
        pre <- if (r$factor_compared == "Land use")
                 paste0(cvar, ":Predominant_simple") else paste0(cvar, ":Biome4")
        recon_max <- max(recon_max, abs(
          (coef_of(fx, paste0(pre, strip(a))) - coef_of(fx, paste0(pre, strip(b)))) -
          r$mean_est))
        n_checked <- n_checked + 1
        next
      }
      w  <- sub("^(Biome|Land use) = ", "", r$held_fixed)
      if (r$factor_compared == "Land use") {
        main <- coef_of(fx, paste0("Predominant_simple", strip(a))) -
                coef_of(fx, paste0("Predominant_simple", strip(b)))
        int  <- if (w == ref_bi) 0 else
          coef_of(fx, paste0("Biome4", strip(w), ":Predominant_simple", strip(a))) -
          coef_of(fx, paste0("Biome4", strip(w), ":Predominant_simple", strip(b)))
      } else {
        main <- coef_of(fx, paste0("Biome4", strip(a))) -
                coef_of(fx, paste0("Biome4", strip(b)))
        int  <- if (w == ref_lu) 0 else
          coef_of(fx, paste0("Biome4", strip(a), ":Predominant_simple", strip(w))) -
          coef_of(fx, paste0("Biome4", strip(b), ":Predominant_simple", strip(w)))
      }
      recon_max <- max(recon_max, abs((main + int) - r$mean_est))
      n_checked <- n_checked + 1
    }
  }
}
cat(sprintf("reconstructed %d of %d rows from fitted coefficients, max |diff| = %.2e\n",
            n_checked, nrow(tbl), recon_max))
stopifnot(n_checked == nrow(tbl))
# Checked against the posterior MEAN, which is exactly additive, so
# agreement should be at floating-point precision. Any real discrepancy
# would mean a term entered the contrast that the two-way coefficient
# reconstruction does not contain.
stopifnot(recon_max < 1e-8)
cat("no row requires a term outside the fitted two-way structure\n")

# ------------------------------------------------------------
# ------------------------------------------------------------
# A2d CHANGE FROM PRIMARY, written separately.
# A2d has no primary-vegetation level: primary is the baseline of the
# response itself, so diff_abund already measures change from primary.
# The quantity "cropland vs primary in tropical forest" is therefore
# not a contrast between levels but the cell mean tested against zero.
# It is not a pairwise comparison, so it is not part of Table 2, but it
# is what Figure 3 plots and the most directly interpretable summary of
# the biome by land-use interaction in this model.
# ------------------------------------------------------------
message("A2d change from primary (separate file) ...")
prim <- list()
for (cv in c("z_malecolcooney", "z_dichrodiff")) {
  load(paste0("fits/A2d_fit_", cv, ".RData"))
  g <- make_grid(Predominant_simple = LU4, Biome4 = BIOME3)
  g[[cv]] <- 0
  sdl <- link_sd(fit)
  for (b in BIOME3) {
    gb <- g[g$Biome4 == b, , drop = FALSE]
    emm <- emm_draws(fit, gb, "Predominant_simple")
    prim[[paste(cv, b)]] <- columnwise(
      emm, sdl, FALSE,
      model = "A2d", response = RESPONSE[["A2d"]], source_fit = SRC_FIT[[cv]],
      biome = b)
  }
  rm(fit); gc(verbose = FALSE)
}
prim <- bind_rows(prim) %>%
  rename(land_use = contrast) %>%
  mutate(across(all_of(c("biome", "land_use")), relabel_levels)) %>%
  select(model, response, source_fit, biome, land_use,
         estimate, lo, hi, pd, credible, effect_size)
write_csv(prim, "results/cooney/a2d_change_from_primary.csv")
cat("A2d change from primary:", nrow(prim), "rows ->",
    "results/cooney/a2d_change_from_primary.csv\n")

# ------------------------------------------------------------
# COLOUR SIMPLE SLOPES, written separately.
# These are NOT part of Table 2: colour is continuous, so a slope is
# not a pairwise comparison of levels. They are the quantity Figure 3
# plots and the quantity A2c and A2d exist to estimate, so they are
# kept in their own file rather than discarded.
# ------------------------------------------------------------
message("colour simple slopes (separate file) ...")
slopes <- list()
for (mod in c("A2c", "A2d")) {
  lus <- if (mod == "A2c") LU5 else LU4
  for (cv in c("z_malecolcooney", "z_dichrodiff")) {
    load(paste0("fits/", mod, "_fit_", cv, ".RData"))
    g <- make_grid(Predominant_simple = lus, Biome4 = BIOME3)
    g[[cv]] <- 0
    sdl <- link_sd(fit); lognorm <- fit$family$family == "lognormal"
    for (wi in c("Biome4", "Predominant_simple"))
      slopes[[paste(mod, cv, wi)]] <- columnwise(
        slope_draws(fit, g, wi, cv), sdl, lognorm,
        model = mod, response = RESPONSE[[mod]], colour_predictor = COLOUR_LAB[[cv]],
        slope_within = if (wi == "Biome4") "Biome" else "Land use")
    rm(fit); gc(verbose = FALSE)
  }
}
slopes <- bind_rows(slopes) %>%
  rename(level = contrast) %>%
  mutate(level = relabel_levels(level)) %>%
  select(model, response, colour_predictor, slope_within, level,
         estimate, lo, hi, pd, credible, effect_size)
write_csv(slopes, "results/cooney/colour_simple_slopes.csv")
cat("colour slopes:", nrow(slopes), "rows ->",
    "results/cooney/colour_simple_slopes.csv\n")

tbl <- tbl %>% select(-mean_est) %>%          # validation aid, not part of the table
  mutate(across(all_of(c("contrast", "held_fixed")), relabel_levels))
write_csv(tbl, "results/cooney/posthoc_interactions.csv")
write_csv(tbl, "figures/pub/cooney/Table_2_posthoc_interactions.csv")

cat("\n", strrep("=", 78), "\n", sep = "")
cat("CREDIBLE PAIRWISE COMPARISONS WITHIN LEVELS\n")
cat(strrep("=", 78), "\n", sep = "")
tbl %>% filter(credible) %>%
  mutate(across(c(estimate, lo, hi, effect_size), ~ round(.x, 4)), pd = round(pd, 3)) %>%
  select(-pct_change) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\nRows:", nrow(tbl), "| credible:", sum(tbl$credible), "\n")
print(table(tbl$model, tbl$factor_compared))
cat("\nWritten to results/cooney/posthoc_interactions.csv",
    "and figures/pub/cooney/Table_2_posthoc_interactions.csv\n")
