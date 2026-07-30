# ============================================================
# A2b_03_diagnostics_summary.R
# Analysis 2b: Diagnostics, summaries, and visualizations
# Paired-difference response (modified minus primary baseline)
# ============================================================

library(dplyr)
library(ggplot2)
library(brms)
library(bayesplot)

# Load data and fits
load("fits/A2b_model_setup.RData")
load("fits/A2b_all_fits.RData")

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

diag_rows <- list()
for (cv in names(fits_A2b)) {
  diag <- check_convergence(fits_A2b[[cv]])
  cat("\n===", cv, "===\n")
  cat("Max Rhat:", round(diag$max_rhat, 4), "\n")
  cat("Min ESS bulk:", round(diag$min_ess_bulk), "\n")
  cat("Min ESS tail:", round(diag$min_ess_tail), "\n")
  cat("Divergences:", diag$n_divergent, "\n")
  cat("Converged:", ifelse(diag$converged, "YES", "NO"), "\n")

  diag_rows[[cv]] <- data.frame(
    colour_predictor = cv,
    max_rhat = diag$max_rhat,
    min_ess_bulk = diag$min_ess_bulk,
    min_ess_tail = diag$min_ess_tail,
    n_divergent = diag$n_divergent,
    converged = diag$converged,
    stringsAsFactors = FALSE
  )
}

diag_summary <- do.call(rbind, diag_rows)
print(diag_summary)
write.csv(diag_summary, "results/cooney/A2b_convergence.csv", row.names = FALSE)

# ============================================================
# MODEL SUMMARIES AND R2
# ============================================================

cat("\n", rep("=", 60), "\n")
cat("MODEL SUMMARIES\n")
cat(rep("=", 60), "\n")

r2_rows <- list()
for (cv in names(fits_A2b)) {
  cat("\n\n==========", cv, "==========\n")
  print(summary(fits_A2b[[cv]]))

  r2 <- bayes_R2(fits_A2b[[cv]])
  cat("\nBayes R2:\n")
  print(r2)

  r2_rows[[cv]] <- data.frame(
    colour_predictor = cv,
    R2_estimate = r2[1, "Estimate"],
    R2_error = r2[1, "Est.Error"],
    R2_lower = r2[1, "Q2.5"],
    R2_upper = r2[1, "Q97.5"],
    stringsAsFactors = FALSE
  )
}

r2_summary <- do.call(rbind, r2_rows)
cat("\n", rep("=", 60), "\n")
cat("R2 SUMMARY\n")
cat(rep("=", 60), "\n")
print(r2_summary)
write.csv(r2_summary, "results/cooney/A2b_r2_summary.csv", row.names = FALSE)

# ============================================================
# FIXED EFFECTS TABLE
# ============================================================

fe_rows <- list()
for (cv in names(fits_A2b)) {
  fe <- as.data.frame(fixef(fits_A2b[[cv]]))
  fe$parameter <- rownames(fe)
  fe$colour_predictor <- cv
  fe_rows[[cv]] <- fe
}

fixed_effects <- do.call(rbind, fe_rows)
rownames(fixed_effects) <- NULL
fixed_effects <- fixed_effects[, c("colour_predictor", "parameter",
                                    "Estimate", "Est.Error",
                                    "Q2.5", "Q97.5")]

write.csv(fixed_effects, "results/cooney/A2b_fixed_effects.csv", row.names = FALSE)
cat("\nFixed effects saved to results/A2b_fixed_effects.csv\n")

# ============================================================
# POSTERIOR PREDICTIVE CHECKS
# ============================================================

cat("\nGenerating posterior predictive checks...\n")

for (cv in names(fits_A2b)) {
  p <- pp_check(fits_A2b[[cv]], ndraws = 50, type = "dens_overlay") +
    ggtitle(paste("PP Check:", cv, "(paired difference)")) +
    theme_minimal()

  ggsave(paste0("figures/cooney/A2b_pp_check_", cv, ".png"), p,
         width = 7, height = 5, dpi = 150)
}

# ============================================================
# CONDITIONAL EFFECTS: Colour x Land-use interaction
# ============================================================

cat("\nGenerating conditional effects plots...\n")

for (cv in names(fits_A2b)) {
  ce <- conditional_effects(fits_A2b[[cv]],
                            effects = paste0(cv, ":Predominant_simple"),
                            resolution = 100)
  p <- plot(ce, plot = FALSE)[[1]] +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = paste("Abundance change ~", cv, "x Land-use"),
         x = paste(cv, "(standardized)"),
         y = "Predicted change in relative abundance\n(vs primary vegetation)") +
    theme_minimal()

  ggsave(paste0("figures/cooney/A2b_conditional_", cv, ".png"), p,
         width = 8, height = 5, dpi = 150)
}

# ============================================================
# FOREST PLOT: Interaction coefficients
# ============================================================

cat("\nGenerating forest plots...\n")

interaction_effects <- fixed_effects[grepl(":", fixed_effects$parameter), ]
interaction_effects$land_use <- sub(".*Predominant_simple", "",
                                     interaction_effects$parameter)
interaction_effects$significant <- sign(interaction_effects$Q2.5) ==
                                   sign(interaction_effects$Q97.5)

for (cv in names(fits_A2b)) {
  df <- interaction_effects[interaction_effects$colour_predictor == cv, ]

  p <- ggplot(df, aes(x = Estimate,
                       y = reorder(land_use, Estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, color = significant),
                   height = 0.2) +
    geom_point(aes(color = significant), size = 2.5) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray50"),
                       name = "95% CI\nexcludes 0") +
    labs(title = paste("Colour x Land-use Interaction:", cv),
         subtitle = "Effect on abundance change (modified - primary)",
         x = "Interaction Effect (95% CI)", y = "Land-use Type") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(paste0("figures/cooney/A2b_forest_", cv, ".png"), p,
         width = 7, height = 4, dpi = 150)
}

# ============================================================
# MAIN EFFECTS: Colour and land-use
# ============================================================

cat("\nGenerating main effects plots...\n")

# Land-use main effects (intercept = mean diff for reference land use
# at mean colour; other levels are offsets)
lu_effects <- fixed_effects[grepl("Predominant_simple",
                                   fixed_effects$parameter) &
                             !grepl(":", fixed_effects$parameter), ]
lu_effects$land_use <- sub("Predominant_simple", "", lu_effects$parameter)
lu_effects$significant <- sign(lu_effects$Q2.5) == sign(lu_effects$Q97.5)

for (cv in names(fits_A2b)) {
  df <- lu_effects[lu_effects$colour_predictor == cv, ]

  p <- ggplot(df, aes(x = Estimate, y = reorder(land_use, Estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, color = significant),
                   height = 0.2) +
    geom_point(aes(color = significant), size = 2.5) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray50"),
                       name = "95% CI\nexcludes 0") +
    labs(title = paste("Land-use main effects:", cv),
         subtitle = "Mean abundance change at mean colour",
         x = "Effect on abundance change (95% CI)", y = "Land-use Type") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(paste0("figures/cooney/A2b_lu_effects_", cv, ".png"), p,
         width = 7, height = 4, dpi = 150)
}

cat("\nA2b analysis complete. Results in results/ and figures/\n")
