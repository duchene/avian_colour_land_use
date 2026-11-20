# ============================================================
# supp_tables.R
# Supplementary Tables S1 and S2, formatted from saved result tables.
#
# Table S1  male versus sex-averaged colourfulness: every reported
#           colourfulness quantity of A1b, A2d and A3b under both metrics
#           (results/cooney/C2_sex_averaged_comparison.csv).
# Table S2  abundance change (A2d) fitted once per colour: slope of
#           abundance change per SD of the colour proportion, overall and
#           within each land use and biome, with each fit's convergence
#           (results/cooney/A2d_colours_{slopes,convergence}.csv).
# Table S3  every phylogenetic quantity under the Claramunt et al. (2025)
#           and Jetz et al. (2012) trees (results/cooney/T1_tree_comparison.csv).
#
# Outputs: figures/pub/cooney/Table_S1_sex_averaged_colourfulness.csv
#          figures/pub/cooney/Table_S2_colour_abundance_slopes.csv
#          figures/pub/cooney/Table_S3_tree_comparison.csv
# ============================================================

out <- "figures/pub/cooney"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
r4 <- function(v) signif(v, 3)

c2 <- read.csv("results/cooney/C2_sex_averaged_comparison.csv", stringsAsFactors = FALSE)
model_lab <- c(A1b = "Species level (A1b)", A2d = "Abundance change (A2d)",
               A3b = "Assemblage level (A3b)")
s1 <- data.frame(
  model = model_lab[c2$model], quantity = c2$quantity,
  within = ifelse(is.na(c2$held_fixed), "", c2$held_fixed), contrast = c2$contrast,
  male_estimate = r4(c2$estimate_male), male_lo = r4(c2$lo_male), male_hi = r4(c2$hi_male),
  male_credible = c2$credible_male,
  sex_averaged_estimate = r4(c2$estimate_sex_averaged), sex_averaged_lo = r4(c2$lo_sex_averaged),
  sex_averaged_hi = r4(c2$hi_sex_averaged), sex_averaged_credible = c2$credible_sex_averaged,
  same_conclusion = c2$same_conclusion)
write.csv(s1, file.path(out, "Table_S1_sex_averaged_colourfulness.csv"), row.names = FALSE)

sl <- read.csv("results/cooney/A2d_colours_slopes.csv", stringsAsFactors = FALSE)
cv <- read.csv("results/cooney/A2d_colours_convergence.csv", stringsAsFactors = FALSE)
s2 <- merge(sl, cv[, c("colour", "max_rhat", "min_ess_bulk")], by = "colour", sort = FALSE)
s2 <- s2[order(match(s2$colour, cv$colour), match(s2$slope_within, c("Overall", "Land use", "Biome"))), ]
s2 <- data.frame(colour = s2$colour, slope_within = s2$slope_within, level = s2$level,
                 estimate = r4(s2$estimate), lo = r4(s2$lo), hi = r4(s2$hi), pd = round(s2$pd, 3),
                 credible = s2$credible, effect_size = r4(s2$effect_size),
                 fit_max_rhat = round(s2$max_rhat, 3), fit_min_ess_bulk = round(s2$min_ess_bulk))
write.csv(s2, file.path(out, "Table_S2_colour_abundance_slopes.csv"), row.names = FALSE)
t1 <- read.csv("results/cooney/T1_tree_comparison.csv", stringsAsFactors = FALSE)
analysis_lab <- c(
  A1b_landuse_contrasts = "Species level (A1b): land use minus primary vegetation",
  A1b_r2_summary = "Species level (A1b): Bayesian R2",
  posthoc_contrasts = "Abundance change (A2d): main-effect contrasts",
  posthoc_interactions = "Abundance change (A2d): contrasts within levels",
  colour_simple_slopes = "Abundance change (A2d): colourfulness and dichromatism slopes",
  A2d_colours_slopes = "Abundance change (A2d): slopes per colour",
  A4_colour_landuse_slopes = "MR-PMM (A4): colour on land-use slopes (SD)",
  A4_correlations = "MR-PMM (A4): correlations",
  A4_phylogenetic_signal = "MR-PMM (A4): phylogenetic signal",
  C1_sample_bias = "Sample bias (C1): PGLS difference (credible = p < 0.05)")
stopifnot(all(t1$table %in% names(analysis_lab)))
s3 <- data.frame(analysis = analysis_lab[t1$table], quantity = t1$key,
                 claramunt_estimate = r4(t1$est_claramunt), claramunt_lo = r4(t1$lo_claramunt),
                 claramunt_hi = r4(t1$hi_claramunt), claramunt_credible = t1$cred_claramunt,
                 jetz_estimate = r4(t1$est_jetz), jetz_lo = r4(t1$lo_jetz), jetz_hi = r4(t1$hi_jetz),
                 jetz_credible = t1$cred_jetz, same_conclusion = t1$same_conclusion)
s3 <- s3[order(match(t1$table, names(analysis_lab))), ]
write.csv(s3, file.path(out, "Table_S3_tree_comparison.csv"), row.names = FALSE)
cat("Saved", file.path(out, c("Table_S1_sex_averaged_colourfulness.csv",
                              "Table_S2_colour_abundance_slopes.csv",
                              "Table_S3_tree_comparison.csv")), sep = "\n")
