# ============================================================
# posthoc_helpers.R
# Shared machinery for the two post-hoc contrast tables.
# Sourced by posthoc_contrasts.R (Table 1, main effects) and
# posthoc_interactions.R (Table 2, simple effects and slopes).
#
# Everything works on posterior draws of the linear predictor with
# random effects excluded (re_formula = NA), so contrasts are
# population-level. No multiplicity correction is applied or needed:
# these are posterior probabilities of a difference, not repeated
# frequentist tests.
# ============================================================

suppressMessages({ library(brms); library(dplyr); library(purrr); library(tibble) })

# SD of the response on the scale the linear predictor lives on.
# Contrasts divided by this are standardised effect sizes.
link_sd <- function(fit) {
  y <- fit$data[[1]]
  if (fit$family$family == "lognormal") sd(log(y)) else sd(y)
}

# Balanced reference grid over the supplied factor levels, numeric
# covariates pinned at the value given.
make_grid <- function(...) {
  expand.grid(..., stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
}

# Draws x levels matrix of estimated marginal means for `target`,
# averaging the linear predictor over every other cell of the grid.
emm_draws <- function(fit, grid, target) {
  lp <- posterior_linpred(fit, newdata = grid, re_formula = NA)
  lv <- unique(grid[[target]])
  out <- vapply(lv, function(l) rowMeans(lp[, grid[[target]] == l, drop = FALSE]),
                numeric(nrow(lp)))
  colnames(out) <- lv
  out
}

# Draws x levels matrix of the slope on `cvar` per 1 SD, within each
# level of `target`. The model is linear in cvar, so the slope is the
# difference of the linear predictor at cvar = +0.5 and -0.5.
slope_draws <- function(fit, grid, target, cvar) {
  g_hi <- g_lo <- grid
  g_hi[[cvar]] <- 0.5
  g_lo[[cvar]] <- -0.5
  lp <- posterior_linpred(fit, newdata = g_hi, re_formula = NA) -
        posterior_linpred(fit, newdata = g_lo, re_formula = NA)
  lv <- unique(grid[[target]])
  out <- vapply(lv, function(l) rowMeans(lp[, grid[[target]] == l, drop = FALSE]),
                numeric(nrow(lp)))
  colnames(out) <- lv
  out
}

# Summarise one vector of contrast draws.
summarise_draws_vec <- function(d, sd_link, lognormal) {
  q <- quantile(d, c(0.025, 0.5, 0.975))
  # mean_est is carried for validation only. The posterior mean of a
  # linear combination equals the same combination of the parameter
  # means, so it reconstructs exactly from fixef(); the median does
  # not. Reporting scripts drop this column.
  tibble(estimate = unname(q[2]), mean_est = mean(d),
         lo = unname(q[1]), hi = unname(q[3]),
         pd = max(mean(d > 0), mean(d < 0)),
         credible = unname(sign(q[1]) == sign(q[3])),
         pct_change = if (lognormal) (exp(unname(q[2])) - 1) * 100 else NA_real_,
         effect_size = unname(q[2]) / sd_link)
}

# Every pairwise difference of the columns of an EMM matrix.
pairwise <- function(emm, sd_link, lognormal, ...) {
  lv <- colnames(emm)
  combos <- t(combn(seq_along(lv), 2))
  map_dfr(seq_len(nrow(combos)), function(k) {
    i <- combos[k, 1]; j <- combos[k, 2]
    # j - i, so a negative estimate means level j sits below level i
    bind_cols(tibble(..., contrast = paste(lv[j], "-", lv[i])),
              summarise_draws_vec(emm[, j] - emm[, i], sd_link, lognormal))
  })
}

# Each column of a matrix summarised on its own (for slopes).
columnwise <- function(mat, sd_link, lognormal, ...) {
  map_dfr(colnames(mat), function(l)
    bind_cols(tibble(..., contrast = l),
              summarise_draws_vec(mat[, l], sd_link, lognormal)))
}

# Display names for reporting. Grids must be built from the DATA level
# names (LU5/LU4/BIOME3 below), so relabelling happens at write time.
DISPLAY <- c("Tropical Forest"  = "Tropical closed",
             "Temperate Forest" = "Temperate closed",
             "Tropical Open"    = "Tropical open",
             "Plantation forest" = "Plantation")
relabel_levels <- function(x) stringr::str_replace_all(as.character(x), DISPLAY)

RESP_LABEL <- c(malecolcooney = "Male colourfulness", dichrodiff = "Sexual dichromatism",
                z_malecolcooney = "Male colourfulness", z_dichrodiff = "Sexual dichromatism")

LU5 <- c("Primary vegetation", "Secondary", "Plantation forest", "Cropland", "Pasture")
LU4 <- c("Secondary", "Plantation forest", "Cropland", "Pasture")
BIOME3 <- c("Tropical Forest", "Temperate Forest", "Tropical Open")
TROPHIC6 <- c("Aquatic predator", "Frugivore", "Granivore", "Invertivore",
              "Nectarivore", "Omnivore")
