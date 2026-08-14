# ============================================================
# B3_03_diagnostics_summary.R
# Analysis B3: diagnostics, summaries, forest and pp-check plots.
# Mirrors A3_03. Outputs results/B3_{convergence,r2_summary,
# fixed_effects}.csv and figures/B3_forest_*, figures/B3_pp_check_*.
# ============================================================

library(tidyverse)
library(brms)
library(bayesplot)

load("fits/B3_model_setup.RData")
load("fits/B3_all_fits.RData")
dir.create("results/dale", showWarnings = FALSE)
dir.create("figures/dale", showWarnings = FALSE)

check_convergence <- function(fit) {
  s <- summary(fit); fe <- s$fixed
  np <- nuts_params(fit)
  n_div <- sum(np$Value[np$Parameter == "divergent__"])
  max_rhat <- max(fe$Rhat, na.rm = TRUE)
  min_ess_b <- min(fe$Bulk_ESS, na.rm = TRUE)
  list(max_rhat = max_rhat, min_ess_bulk = min_ess_b,
       min_ess_tail = min(fe$Tail_ESS, na.rm = TRUE), n_divergent = n_div,
       converged = max_rhat < 1.01 & min_ess_b > 400 & n_div == 0)
}

diag_summary <- map_dfr(names(fits_B3), function(resp) {
  d <- check_convergence(fits_B3[[resp]])
  tibble(response = resp, max_rhat = d$max_rhat, min_ess_bulk = d$min_ess_bulk,
         min_ess_tail = d$min_ess_tail, n_divergent = d$n_divergent,
         converged = d$converged)
})
print(diag_summary)
write_csv(diag_summary, "results/dale/B3_convergence.csv")

r2_summary <- map_dfr(names(fits_B3), function(resp) {
  r2 <- bayes_R2(fits_B3[[resp]])
  tibble(response = resp, R2_estimate = r2[1, "Estimate"], R2_error = r2[1, "Est.Error"],
         R2_lower = r2[1, "Q2.5"], R2_upper = r2[1, "Q97.5"])
})
write_csv(r2_summary, "results/dale/B3_r2_summary.csv")
cat("\nBayesian R2:\n"); print(r2_summary)

fixed_effects <- map_dfr(names(fits_B3), function(resp) {
  fe <- as.data.frame(fixef(fits_B3[[resp]]))
  fe$parameter <- rownames(fe); fe$response <- resp
  fe %>% select(response, parameter, everything())
})
write_csv(fixed_effects, "results/dale/B3_fixed_effects.csv")
cat("\nFixed effects saved to results/B3_fixed_effects.csv\n")

for (resp in names(fits_B3)) {
  df <- fixed_effects %>%
    filter(response == resp, parameter != "Intercept") %>%
    mutate(significant = sign(Q2.5) == sign(Q97.5),
           param_short = gsub("prop_", "", parameter))
  p <- ggplot(df, aes(x = Estimate, y = reorder(param_short, Estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, color = significant), height = 0.2) +
    geom_point(aes(color = significant), size = 2.5) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray50"),
                       name = "95% CI\nexcludes 0") +
    labs(title = paste("Land-use proportion effects (Dale):", resp),
         subtitle = "Phylogenetic regression",
         x = "Effect estimate (95% CI)", y = "Land-use proportion predictor") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(paste0("figures/dale/B3_forest_", resp, ".png"), p, width = 7, height = 4, dpi = 150)
}

for (resp in names(fits_B3)) {
  p <- pp_check(fits_B3[[resp]], ndraws = 50, type = "dens_overlay") +
    ggtitle(paste("PP check (Dale):", resp)) + theme_minimal()
  ggsave(paste0("figures/dale/B3_pp_check_", resp, ".png"), p, width = 7, height = 5, dpi = 150)
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

contrasts_tbl <- map_dfr(names(fits_B3), function(resp) {
  d  <- posterior::as_draws_matrix(fits_B3[[resp]])
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
write_csv(contrasts_tbl, "results/dale/B3_landuse_contrasts.csv")
print(as.data.frame(contrasts_tbl), digits = 3)

cat("\nB3 analysis complete. Results in results/ and figures/\n")
