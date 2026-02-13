# ============================================================
# A2_03_diagnostics_summary.R
# Analysis 2: Diagnostics, summaries, and visualizations
# ============================================================

library(tidyverse)
library(brms)
library(bayesplot)
library(loo)

# Load data and fits
load("fits/A2_model_setup.RData")
load("fits/A2_all_fits.RData")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

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

diag_summary <- map_dfr(names(fits_A2), function(cv) {
  diag <- check_convergence(fits_A2[[cv]])
  cat("\n===", cv, "===\n")
  cat("Max Rhat:", round(diag$max_rhat, 4), "\n")
  cat("Min ESS bulk:", round(diag$min_ess_bulk), "\n")
  cat("Min ESS tail:", round(diag$min_ess_tail), "\n")
  cat("Divergences:", diag$n_divergent, "\n")
  cat("Converged:", ifelse(diag$converged, "YES", "NO"), "\n")

  tibble(
    colour_predictor = cv,
    max_rhat = diag$max_rhat,
    min_ess_bulk = diag$min_ess_bulk,
    min_ess_tail = diag$min_ess_tail,
    n_divergent = diag$n_divergent,
    converged = diag$converged
  )
})

print(diag_summary)
write_csv(diag_summary, "results/A2_convergence.csv")

# ============================================================
# MODEL SUMMARIES AND R2
# ============================================================

cat("\n", rep("=", 60), "\n")
cat("MODEL SUMMARIES\n")
cat(rep("=", 60), "\n")

r2_summary <- map_dfr(names(fits_A2), function(cv) {
  cat("\n\n==========", cv, "==========\n")
  print(summary(fits_A2[[cv]]))

  r2 <- bayes_R2(fits_A2[[cv]])
  cat("\nBayes R2:\n")
  print(r2)

  tibble(
    colour_predictor = cv,
    R2_estimate = r2[1, "Estimate"],
    R2_error = r2[1, "Est.Error"],
    R2_lower = r2[1, "Q2.5"],
    R2_upper = r2[1, "Q97.5"]
  )
})

cat("\n", rep("=", 60), "\n")
cat("R2 SUMMARY\n")
cat(rep("=", 60), "\n")
print(r2_summary)
write_csv(r2_summary, "results/A2_r2_summary.csv")

# ============================================================
# FIXED EFFECTS TABLE
# ============================================================

fixed_effects <- map_dfr(names(fits_A2), function(cv) {
  fe <- as.data.frame(fixef(fits_A2[[cv]]))
  fe$parameter <- rownames(fe)
  fe$colour_predictor <- cv
  fe %>% select(colour_predictor, parameter, everything())
})

write_csv(fixed_effects, "results/A2_fixed_effects.csv")
cat("\nFixed effects saved to results/A2_fixed_effects.csv\n")

# ============================================================
# POSTERIOR PREDICTIVE CHECKS
# ============================================================

cat("\nGenerating posterior predictive checks...\n")

for (cv in names(fits_A2)) {
  p <- pp_check(fits_A2[[cv]], ndraws = 50, type = "dens_overlay") +
    ggtitle(paste("PP Check:", cv)) +
    coord_cartesian(xlim = c(-5, 10)) +
    theme_minimal()

  ggsave(paste0("figures/A2_pp_check_", cv, ".png"), p,
         width = 7, height = 5, dpi = 150)
}

# ============================================================
# CONDITIONAL EFFECTS: Colour x Land-use interaction
# ============================================================

cat("\nGenerating conditional effects plots...\n")

for (cv in names(fits_A2)) {
  ce <- conditional_effects(fits_A2[[cv]],
                            effects = paste0(cv, ":Predominant_simple"),
                            resolution = 100)
  p <- plot(ce, plot = FALSE)[[1]] +
    labs(title = paste("Abundance ~", cv, "x Land-use"),
         x = paste(cv, "(standardized)"),
         y = "Predicted Abundance") +
    theme_minimal()

  ggsave(paste0("figures/A2_conditional_", cv, ".png"), p,
         width = 8, height = 5, dpi = 150)
}

# ============================================================
# FOREST PLOT: Interaction coefficients
# ============================================================

cat("\nGenerating forest plots...\n")

interaction_effects <- fixed_effects %>%
  filter(grepl(":", parameter)) %>%
  mutate(
    land_use = str_extract(parameter, "Predominant_simple[A-Za-z]+"),
    land_use = str_remove(land_use, "Predominant_simple"),
    significant = sign(Q2.5) == sign(Q97.5)
  )

for (cv in names(fits_A2)) {
  df <- interaction_effects %>% filter(colour_predictor == cv)

  p <- ggplot(df, aes(x = Estimate, y = reorder(land_use, Estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, color = significant), height = 0.2) +
    geom_point(aes(color = significant), size = 2.5) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray50"),
                       name = "95% CI\nexcludes 0") +
    labs(title = paste("Colour x Land-use Interaction:", cv),
         x = "Interaction Effect (95% CI)", y = "Land-use Type") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(paste0("figures/A2_forest_", cv, ".png"), p,
         width = 7, height = 4, dpi = 150)
}

# ============================================================
# COMPARISON ACROSS COLOUR VARIABLES
# ============================================================

cat("\nGenerating comparison plot...\n")

colour_main <- fixed_effects %>%
  filter(grepl("^z_", parameter) & !grepl(":", parameter)) %>%
  mutate(significant = sign(Q2.5) == sign(Q97.5))

p_compare <- ggplot(colour_main, aes(x = Estimate, y = colour_predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, color = significant), height = 0.2) +
  geom_point(aes(color = significant), size = 3) +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray50"),
                     name = "95% CI\nexcludes 0") +
  labs(title = "Effect of Colour on Abundance (across metrics)",
       subtitle = "Main effect of standardized colour predictor",
       x = "Effect on log(Abundance)", y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/A2_colour_comparison.png", p_compare,
       width = 7, height = 4, dpi = 150)

cat("\nA2 analysis complete. Results in results/ and figures/\n")
