# ============================================================
# A1_05_summarize.R
# Analysis 1: Generate summary tables and figures
# ============================================================

library(tidyverse)
library(brms)
library(bayesplot)

# Load data
load("fits/A1_model_setup.RData")
load("fits/A1_all_fits.RData")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ============================================================
# MODEL SUMMARIES
# ============================================================

cat("\n", rep("=", 60), "\n")
cat("MODEL SUMMARIES (REDUCED MODELS)\n")
cat(rep("=", 60), "\n")

model_summaries <- list()

for (resp in responses) {
  cat("\n\n==========", resp, "==========\n")
  print(summary(fits_reduced[[resp]]))

  r2 <- bayes_R2(fits_reduced[[resp]])
  cat("\nBayes R2:\n")
  print(r2)

  model_summaries[[resp]] <- list(
    summary = summary(fits_reduced[[resp]]),
    r2 = r2
  )
}

# ============================================================
# RANDOM EFFECT VARIANCE
# ============================================================

cat("\n", rep("=", 60), "\n")
cat("RANDOM EFFECT (SSBS) VARIANCE ESTIMATES\n")
cat(rep("=", 60), "\n")

re_variance <- map_dfr(responses, function(resp) {
  re <- summary(fits_reduced[[resp]])$random$SSBS
  tibble(
    response = resp,
    sd_estimate = re$Estimate,
    sd_error = re$Est.Error,
    sd_lower = re$`l-95% CI`,
    sd_upper = re$`u-95% CI`
  )
})

print(re_variance)
write_csv(re_variance, "results/A1_random_effect_variance.csv")

# ============================================================
# FIXED EFFECTS TABLE
# ============================================================

extract_fixed_effects <- function(fit, resp_name) {
  fe <- as.data.frame(summary(fit)$fixed)
  fe$parameter <- rownames(fe)
  fe$response <- resp_name
  fe %>% select(response, parameter, everything())
}

fixed_effects <- map_dfr(responses, function(resp) {
  extract_fixed_effects(fits_reduced[[resp]], resp)
})

write_csv(fixed_effects, "results/A1_fixed_effects_all.csv")
cat("\nFixed effects saved to results/A1_fixed_effects_all.csv\n")

# ============================================================
# R2 SUMMARY
# ============================================================

r2_summary <- map_dfr(responses, function(resp) {
  r2 <- bayes_R2(fits_reduced[[resp]])
  tibble(
    response = resp,
    R2_estimate = r2[1, "Estimate"],
    R2_error = r2[1, "Est.Error"],
    R2_lower = r2[1, "Q2.5"],
    R2_upper = r2[1, "Q97.5"]
  )
})

cat("\n", rep("=", 60), "\n")
cat("BAYES R2 SUMMARY\n")
cat(rep("=", 60), "\n")
print(r2_summary)

write_csv(r2_summary, "results/A1_bayes_r2_summary.csv")

# ============================================================
# COEFFICIENT PLOTS
# ============================================================

cat("\nGenerating coefficient plots...\n")

for (resp in responses) {
  pars <- variables(fits_reduced[[resp]])
  pars_to_plot <- pars[grepl("^b_Predominant_simple[A-Z]|^b_z_logMass$", pars)]

  if (length(pars_to_plot) > 0) {
    p <- mcmc_intervals(fits_reduced[[resp]], pars = pars_to_plot) +
      ggtitle(paste("Main Effects:", resp)) +
      theme_minimal()

    ggsave(paste0("figures/A1_coef_plot_", resp, ".png"), p,
           width = 8, height = 6, dpi = 150)
  }
}

# ============================================================
# INTERACTION HEATMAPS
# ============================================================

cat("\nGenerating interaction heatmaps...\n")

extract_trophic_interactions <- function(fit, resp_name) {
  fe <- as.data.frame(fixef(fit))
  fe$param <- rownames(fe)

  fe %>%
    filter(grepl("Predominant_simple.*:Trophic.Niche", param)) %>%
    mutate(
      land_use = str_extract(param, "Predominant_simple[A-Za-z]+"),
      land_use = str_remove(land_use, "Predominant_simple"),
      trophic = str_extract(param, "Trophic.Niche[A-Za-z]+"),
      trophic = str_remove(trophic, "Trophic.Niche"),
      response = resp_name
    ) %>%
    select(response, land_use, trophic, Estimate, Q2.5, Q97.5)
}

extract_biome_interactions <- function(fit, resp_name) {
  fe <- as.data.frame(fixef(fit))
  fe$param <- rownames(fe)

  fe %>%
    filter(grepl("Predominant_simple.*:Biome4", param)) %>%
    mutate(
      land_use = str_extract(param, "Predominant_simple[A-Za-z]+"),
      land_use = str_remove(land_use, "Predominant_simple"),
      biome = str_extract(param, "Biome4[A-Za-z ]+"),
      biome = str_remove(biome, "Biome4"),
      response = resp_name
    ) %>%
    select(response, land_use, biome, Estimate, Q2.5, Q97.5)
}

# Extract interactions for all responses
trophic_effects <- map_dfr(responses, ~extract_trophic_interactions(fits_reduced[[.x]], .x))
biome_effects <- map_dfr(responses, ~extract_biome_interactions(fits_reduced[[.x]], .x))

# Land-use x Trophic heatmap (combined)
p_trophic <- ggplot(trophic_effects, aes(x = trophic, y = land_use, fill = Estimate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Estimate)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Effect") +
  facet_wrap(~response, ncol = 2) +
  labs(title = "Land-use x Trophic Niche Interactions",
       x = "Trophic Niche", y = "Land-use Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold"))

ggsave("figures/A1_heatmap_trophic_interactions.png", p_trophic,
       width = 12, height = 8, dpi = 150)

# Land-use x Biome4 heatmap (combined)
p_biome <- ggplot(biome_effects, aes(x = biome, y = land_use, fill = Estimate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Estimate)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Effect") +
  facet_wrap(~response, ncol = 2) +
  labs(title = "Land-use x Biome Interactions",
       x = "Biome", y = "Land-use Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold"))

ggsave("figures/A1_heatmap_biome_interactions.png", p_biome,
       width = 12, height = 8, dpi = 150)

# Save interaction data
write_csv(trophic_effects, "results/A1_trophic_interactions.csv")
write_csv(biome_effects, "results/A1_biome_interactions.csv")

# ============================================================
# CONDITIONAL EFFECTS (Land-use x Mass)
# ============================================================

cat("\nGenerating conditional effects plots...\n")

for (resp in responses) {
  ce_mass <- conditional_effects(fits_reduced[[resp]],
                                  effects = "z_logMass:Predominant_simple",
                                  resolution = 100)
  p <- plot(ce_mass, plot = FALSE)[[1]] +
    labs(title = paste(resp, "- Body Mass x Land-use"),
         x = "Standardized log(Mass)", y = "Predicted Value") +
    theme_minimal()

  ggsave(paste0("figures/A1_conditional_mass_", resp, ".png"), p,
         width = 8, height = 5, dpi = 150)
}

# ============================================================
# FOREST PLOT OF MAIN EFFECTS
# ============================================================

cat("\nGenerating forest plots...\n")

main_effects_df <- map_dfr(responses, function(resp) {
  fe <- as.data.frame(fixef(fits_reduced[[resp]]))
  fe$param <- rownames(fe)
  fe %>%
    filter(!grepl(":", param), param != "Intercept") %>%
    mutate(
      response = resp,
      significant = sign(Q2.5) == sign(Q97.5),
      category = case_when(
        grepl("Predominant", param) ~ "Land-use",
        grepl("Biome4", param) ~ "Biome",
        grepl("Trophic", param) ~ "Trophic",
        param == "z_logMass" ~ "Mass",
        TRUE ~ "Other"
      ),
      param_short = str_remove_all(param, "Predominant_simple|Biome4|Trophic.Niche")
    )
})

for (resp in responses) {
  df <- main_effects_df %>%
    filter(response == resp, category %in% c("Land-use", "Trophic", "Mass"))

  p <- ggplot(df, aes(x = Estimate, y = reorder(param_short, Estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, color = significant), height = 0.2) +
    geom_point(aes(color = significant), size = 2) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray50"),
                       name = "95% CI\nexcludes 0") +
    facet_wrap(~category, scales = "free_y", ncol = 1) +
    labs(title = paste("Main Effects:", resp),
         x = "Effect Estimate (95% CI)", y = NULL) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

  ggsave(paste0("figures/A1_forest_", resp, ".png"), p,
         width = 7, height = 8, dpi = 150)
}

cat("\nSummary complete. Results in results/ and figures/\n")
