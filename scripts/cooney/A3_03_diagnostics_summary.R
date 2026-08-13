# ============================================================
# A3_03_diagnostics_summary.R
# Analysis 3: Diagnostics, summaries, and visualizations
# ============================================================

library(tidyverse)
library(brms)
library(bayesplot)

# Load data and fits
load("fits/A3_model_setup.RData")
load("fits/A3_all_fits.RData")

dir.create("results/cooney", showWarnings = FALSE)
dir.create("figures/cooney", showWarnings = FALSE)

# ============================================================
# CONVERGENCE DIAGNOSTICS
# ============================================================

cat("\n", rep("=", 60), "\n")
cat("CONVERGENCE DIAGNOSTICS\n")
cat(rep("=", 60), "\n")

check_convergence <- function(fit) {
  s <- summary(fit)
  fe <- s$fixed
  max_rhat   <- max(fe$Rhat, na.rm = TRUE)
  min_ess_b  <- min(fe$Bulk_ESS, na.rm = TRUE)
  min_ess_t  <- min(fe$Tail_ESS, na.rm = TRUE)
  np <- nuts_params(fit)
  n_div <- sum(np$Value[np$Parameter == "divergent__"])

  list(max_rhat = max_rhat, min_ess_bulk = min_ess_b,
       min_ess_tail = min_ess_t, n_divergent = n_div,
       converged = max_rhat < 1.01 & min_ess_b > 400 & n_div == 0)
}

diag_summary <- map_dfr(names(fits_A3), function(resp) {
  diag <- check_convergence(fits_A3[[resp]])
  cat("\n===", resp, "===\n")
  cat("Max Rhat:", round(diag$max_rhat, 4), "\n")
  cat("Min ESS bulk:", round(diag$min_ess_bulk), "\n")
  cat("Min ESS tail:", round(diag$min_ess_tail), "\n")
  cat("Divergences:", diag$n_divergent, "\n")
  cat("Converged:", ifelse(diag$converged, "YES", "NO"), "\n")

  tibble(
    response = resp,
    max_rhat = diag$max_rhat,
    min_ess_bulk = diag$min_ess_bulk,
    min_ess_tail = diag$min_ess_tail,
    n_divergent = diag$n_divergent,
    converged = diag$converged
  )
})

print(diag_summary)
write_csv(diag_summary, "results/cooney/A3_convergence.csv")

# ============================================================
# MODEL SUMMARIES AND R2
# ============================================================

cat("\n", rep("=", 60), "\n")
cat("MODEL SUMMARIES\n")
cat(rep("=", 60), "\n")

r2_summary <- map_dfr(names(fits_A3), function(resp) {
  cat("\n\n==========", resp, "==========\n")
  print(summary(fits_A3[[resp]]))

  r2 <- bayes_R2(fits_A3[[resp]])
  cat("\nBayes R2:\n")
  print(r2)

  tibble(
    response = resp,
    R2_estimate = r2[1, "Estimate"],
    R2_error = r2[1, "Est.Error"],
    R2_lower = r2[1, "Q2.5"],
    R2_upper = r2[1, "Q97.5"]
  )
})

write_csv(r2_summary, "results/cooney/A3_r2_summary.csv")

# ============================================================
# FIXED EFFECTS TABLE
# ============================================================

fixed_effects <- map_dfr(names(fits_A3), function(resp) {
  fe <- as.data.frame(fixef(fits_A3[[resp]]))
  fe$parameter <- rownames(fe)
  fe$response <- resp
  fe %>% select(response, parameter, everything())
})

write_csv(fixed_effects, "results/cooney/A3_fixed_effects.csv")
cat("\nFixed effects saved to results/A3_fixed_effects.csv\n")

# ============================================================
# FOREST PLOTS
# ============================================================

cat("\nGenerating forest plots...\n")

for (resp in names(fits_A3)) {
  df <- fixed_effects %>%
    filter(response == resp, parameter != "Intercept") %>%
    mutate(
      significant = sign(Q2.5) == sign(Q97.5),
      param_short = gsub("prop_", "", parameter)
    )

  p <- ggplot(df, aes(x = Estimate, y = reorder(param_short, Estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, color = significant), height = 0.2) +
    geom_point(aes(color = significant), size = 2.5) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray50"),
                       name = "95% CI\nexcludes 0") +
    labs(title = paste("Land-use Proportion Effects:", resp),
         subtitle = "Phylogenetic regression",
         x = "Effect Estimate (95% CI)",
         y = "Land-use Proportion Predictor") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(paste0("figures/cooney/A3_forest_", resp, ".png"), p,
         width = 7, height = 4, dpi = 150)
}

# ============================================================
# POSTERIOR PREDICTIVE CHECKS
# ============================================================

cat("\nGenerating posterior predictive checks...\n")

for (resp in names(fits_A3)) {
  p <- pp_check(fits_A3[[resp]], ndraws = 50, type = "dens_overlay") +
    ggtitle(paste("PP Check:", resp)) +
    theme_minimal()

  ggsave(paste0("figures/cooney/A3_pp_check_", resp, ".png"), p,
         width = 7, height = 5, dpi = 150)
}

# ============================================================
# LAND-USE CONTRASTS AGAINST PRIMARY VEGETATION
# ============================================================
# The five compositional coefficients are strongly correlated in the
# posterior, so comparing their marginal credible intervals understates
# the evidence for a difference. Contrast the posterior draws instead.
# This is the only valid way to ask whether species associated with a
# given land use differ in colour from primary-vegetation species.

cat("\nLand-use contrasts against primary vegetation...\n")

lu_contrast <- c("prop_Cropland", "prop_Pasture",
                 "prop_Plantation_forest", "prop_Secondary")

contrasts_tbl <- map_dfr(names(fits_A3), function(resp) {
  d  <- posterior::as_draws_matrix(fits_A3[[resp]])
  bp <- d[, "b_prop_Primary_vegetation"]
  map_dfr(lu_contrast, function(lu) {
    delta <- d[, paste0("b_", lu)] - bp
    q <- quantile(delta, c(0.025, 0.5, 0.975))
    tibble(response = resp,
           contrast = paste0(sub("prop_", "", lu), " - Primary vegetation"),
           Estimate = unname(q[2]), Q2.5 = unname(q[1]), Q97.5 = unname(q[3]),
           prob_positive = mean(delta > 0),
           credible = unname(sign(q[1]) == sign(q[3])))
  })
})
write_csv(contrasts_tbl, "results/cooney/A3_landuse_contrasts.csv")
print(as.data.frame(contrasts_tbl), digits = 3)

cat("\nA3 analysis complete. Results in results/ and figures/\n")
