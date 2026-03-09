# ============================================================
# 00_descriptive_figures.R
# Descriptive figures: colour distributions by land use,
# trophic niche, and their combination
#
# All plots use log-scale y-axes, matching the lognormal
# regression models. For dichromatism (a ratio), a reference
# line at 1 (equal male/female colour) is shown.
# ============================================================

library(tidyverse)

# ============================================================
# LOAD DATA
# ============================================================

dat <- read.csv("data/present.csv", stringsAsFactors = FALSE)
cat("Records:", nrow(dat), "\n")

# Order land use with Primary vegetation first
dat$Predominant_simple <- factor(
  dat$Predominant_simple,
  levels = c("Primary vegetation", "Secondary", "Cropland",
             "Pasture", "Plantation forest")
)

# Colour palette for land-use types
lu_cols <- c(
  "Primary vegetation" = "#228B22",
  "Secondary"          = "#9ACD32",
  "Cropland"           = "#DAA520",
  "Pasture"            = "#CD853F",
  "Plantation forest"  = "#8B4513"
)

dir.create("figures", showWarnings = FALSE)

# ============================================================
# HELPERS
# ============================================================

colour_vars <- c(
  meancolcooney  = "Mean plumage colour (LociUVS)",
  malecolcooney  = "Male plumage colour (LociUVS)",
  dichrocooney   = "Sexual dichromatism (male/female)",
  dichrodiff     = "Sexual dichromatism (male - female)"
)

# Sensible y-axis per variable
add_scale <- function(p, varname) {
  if (varname == "dichrocooney") {
    p <- p +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
      scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4),
                    labels = c("0.5", "0.75", "1", "1.5", "2", "3", "4"))
  } else if (varname == "dichrodiff") {
    # Difference can be negative — use linear scale with reference at 0
    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40")
  } else {
    p <- p +
      scale_y_log10(breaks = c(30, 50, 75, 100, 150, 200, 300, 500),
                    labels = c("30", "50", "75", "100", "150", "200", "300", "500"))
  }
  p
}

# ============================================================
# 1. BOXPLOTS: COLOUR ~ LAND USE
# ============================================================

for (v in names(colour_vars)) {
  d <- dat %>% filter(!is.na(.data[[v]]), if (v != "dichrodiff") .data[[v]] > 0 else TRUE)

  p <- ggplot(d, aes(x = Predominant_simple, y = .data[[v]],
                      fill = Predominant_simple)) +
    geom_violin(alpha = 0.4, linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.7) +
    scale_fill_manual(values = lu_cols) +
    labs(x = NULL, y = colour_vars[v],
         title = paste(colour_vars[v], "by land use")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))

  p <- add_scale(p, v)

  ggsave(paste0("figures/desc_landuse_", v, ".png"), p,
         width = 7, height = 5, dpi = 200)
  cat("Saved:", v, "by land use\n")
}

# ============================================================
# 2. BOXPLOTS: COLOUR ~ TROPHIC NICHE
# ============================================================

# Order trophic niches by median mean colour
niche_order <- dat %>%
  filter(!is.na(meancolcooney), !is.na(Trophic.Niche)) %>%
  group_by(Trophic.Niche) %>%
  summarise(med = median(meancolcooney)) %>%
  arrange(desc(med)) %>%
  pull(Trophic.Niche)

dat$Trophic.Niche <- factor(dat$Trophic.Niche, levels = niche_order)

for (v in names(colour_vars)) {
  d <- dat %>% filter(!is.na(.data[[v]]), !is.na(Trophic.Niche),
                      if (v != "dichrodiff") .data[[v]] > 0 else TRUE)

  p <- ggplot(d, aes(x = Trophic.Niche, y = .data[[v]],
                      fill = Trophic.Niche)) +
    geom_violin(alpha = 0.4, linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.7) +
    labs(x = NULL, y = colour_vars[v],
         title = paste(colour_vars[v], "by trophic niche")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1))

  p <- add_scale(p, v)

  ggsave(paste0("figures/desc_trophic_", v, ".png"), p,
         width = 8, height = 5, dpi = 200)
  cat("Saved:", v, "by trophic niche\n")
}

# ============================================================
# 3. FACETED: COLOUR ~ LAND USE x TROPHIC NICHE
# ============================================================
# Focus on the 4 most common trophic niches

top_niches <- dat %>%
  filter(!is.na(Trophic.Niche)) %>%
  count(Trophic.Niche, sort = TRUE) %>%
  head(4) %>%
  pull(Trophic.Niche)

d_facet <- dat %>%
  filter(!is.na(meancolcooney), Trophic.Niche %in% top_niches,
         meancolcooney > 0) %>%
  mutate(Trophic.Niche = factor(Trophic.Niche, levels = top_niches))

p <- ggplot(d_facet, aes(x = Predominant_simple, y = meancolcooney,
                          fill = Predominant_simple)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  facet_wrap(~ Trophic.Niche) +
  scale_fill_manual(values = lu_cols) +
  scale_y_log10() +
  labs(x = NULL, y = "Mean plumage colour (LociUVS, log scale)",
       title = "Mean colour by land use and trophic niche") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("figures/desc_landuse_x_trophic_meancolcooney.png", p,
       width = 10, height = 7, dpi = 200)
cat("Saved: land use x trophic niche panel\n")

# ============================================================
# 4. SPECIES-LEVEL: MEAN COLOUR BY LAND USE
# ============================================================
# Each species appears once, assigned to its most frequent land use

sp_summary <- dat %>%
  filter(!is.na(meancolcooney)) %>%
  count(Best_guess_binomial, Predominant_simple, meancolcooney,
        malecolcooney, dichrocooney, dichrodiff) %>%
  group_by(Best_guess_binomial) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("\nSpecies-level: assigned", nrow(sp_summary),
    "species to dominant land use\n")
print(table(sp_summary$Predominant_simple))

for (v in names(colour_vars)) {
  d <- sp_summary %>% filter(!is.na(.data[[v]]), if (v != "dichrodiff") .data[[v]] > 0 else TRUE)

  p <- ggplot(d, aes(x = Predominant_simple, y = .data[[v]],
                      fill = Predominant_simple)) +
    geom_violin(alpha = 0.4, linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.7) +
    scale_fill_manual(values = lu_cols) +
    labs(x = NULL, y = colour_vars[v],
         title = paste(colour_vars[v], "by dominant land use (species-level)")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))

  p <- add_scale(p, v)

  ggsave(paste0("figures/desc_species_landuse_", v, ".png"), p,
         width = 7, height = 5, dpi = 200)
  cat("Saved: species-level", v, "\n")
}

# ============================================================
# 5. BODY MASS vs COLOUR (scatter, log-log)
# ============================================================

d <- dat %>%
  filter(!is.na(meancolcooney), !is.na(Mass), Mass > 0,
         meancolcooney > 0) %>%
  distinct(Best_guess_binomial, .keep_all = TRUE)

p <- ggplot(d, aes(x = log(Mass), y = meancolcooney,
                    colour = Predominant_simple)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = lu_cols, name = "Land use") +
  scale_y_log10(breaks = c(30, 50, 75, 100, 150, 200, 300, 500),
                labels = c("30", "50", "75", "100", "150", "200", "300", "500")) +
  labs(x = "log(Body mass)", y = "Mean plumage colour (LociUVS, log scale)",
       title = "Body mass vs. mean colour by land use (species-level)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("figures/desc_mass_vs_colour.png", p,
       width = 8, height = 6, dpi = 200)
cat("Saved: mass vs colour scatter\n")

cat("\nAll descriptive figures complete.\n")
