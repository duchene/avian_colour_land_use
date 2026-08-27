# ============================================================
# pub_figures.R
# Publication figures for the avian colour x land-use manuscript.
#
# Figures 1 and 2 are introductory, made separately and held outside
# this repository. This script produces 3, 5 and 6.
# Figure 3  A1  species level
# Figure 5  A2d abundance change
# Figure 6  A3b community colour
#
# Each figure is two rows by three columns. A row is one response:
# male colourfulness on top, sexual dichromatism below.
# Column 1 carries every coefficient of that model in a single
# forest, main effects and interactions together, separated by
# rules. Columns 2 and 3 carry the raw data behind the main effects
# and behind the interactions respectively.
#
# Mean colour and the dichromatism ratio go to the supplement.
#
# Pairwise contrasts are in Table 1 (posthoc_contrasts.R) and
# Table 2 (posthoc_interactions.R).
#
# Run from repo root: Rscript scripts/cooney/pub_figures.R
# ============================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr)
  library(patchwork); library(scales); library(ape); library(phytools)
})

theme_set(theme_bw(base_size = 9))
mm2in <- function(x) x / 25.4
OUT <- "figures/pub/cooney"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 400)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
  message("Saved ", name)
}
try_fig <- function(expr)
  tryCatch(expr, error = function(e) message("  !! FAILED: ", conditionMessage(e)))

# ---- shared palettes and labels ---------------------------------------------
lu_levels <- c("Primary vegetation", "Secondary", "Plantation", "Cropland", "Pasture")
lu_short  <- c("Primary vegetation" = "Primary", "Secondary" = "Secondary",
               "Plantation" = "Plantation", "Cropland" = "Cropland",
               "Pasture" = "Pasture")
lu_display <- setNames(lu_levels, lu_levels)
lu_pal <- c("Primary vegetation" = "#1B7837", "Secondary" = "#7FBC41",
            "Plantation" = "#FEC44F", "Cropland" = "#D95F0E",
            "Pasture" = "#8C3B08")
# Display names are now canonical: the post-hoc tables are written with
# them. Only the raw data files still carry the original level names, so
# they are converted on read.
biome_levels <- c("Tropical closed", "Temperate closed", "Tropical open")
relabel <- function(x) str_replace_all(as.character(x),
  c("Tropical Forest" = "Tropical closed", "Temperate Forest" = "Temperate closed",
    "Tropical Open" = "Tropical open", "Plantation forest" = "Plantation"))
biome_display <- setNames(biome_levels, biome_levels)
biome_pal <- c("Tropical closed" = "#1F5C99", "Temperate closed" = "#7FB1DE",
               "Tropical open" = "#C9A227")

resp_levels <- c("Male colourfulness", "Sexual dichromatism")
resp_pal <- c("Male colourfulness" = "#1F5C99", "Sexual dichromatism" = "#B2456E")
resp_ramp <- function(rsp)
  setNames(colorRampPalette(c("#FFFFFF", resp_pal[[rsp]]))(5)[c(2, 3, 5)],
           c("low", "mid", "high"))

# Range the drawn boxplot whiskers actually cover. geom_boxplot with
# outlier.shape = NA hides the outliers but ggplot still scales the axis
# to them, so a handful of extreme species stretch the panel and squash
# every box into a band. This returns the 1.5 x IQR reach per group.
whisker_lim <- function(v, g) {
  r <- tapply(v, g, function(z) {
    z <- z[!is.na(z)]
    if (!length(z)) return(c(NA, NA))
    q <- quantile(z, c(0.25, 0.75)); w <- 1.5 * diff(q)
    c(min(z[z >= q[1] - w]), max(z[z <= q[2] + w]))
  })
  range(unlist(r), na.rm = TRUE)
}

pretty_lu <- function(x) str_replace(x, "Plantationforest", "Plantation forest")
pretty_biome <- function(x) x |>
  str_replace("TemperateForest", "Temperate Forest") |>
  str_replace("TropicalOpen", "Tropical Open") |>
  str_replace("TropicalForest", "Tropical Forest")

# Grouped forest of posterior CONTRASTS, not model coefficients. A
# coefficient is conditional on the reference levels of everything
# else, so a land-use coefficient describes tropical forest alone.
# The contrasts plotted here are the quantities the text reports and
# the raw panels beside them show. `d` needs columns estimate, lo,
# hi, credible, label and block, and blocks become left-hand strips
# sized to their row count.
forest_panel <- function(d, fill_col, title, xlab) {
  # A label can appear in more than one block ("Cropland" is both a
  # colour-slope row and a land-use contrast), so key the axis by block
  # and label together and strip the prefix back off for display.
  d <- d[order(d$block, d$label), ]
  k <- paste0(as.integer(d$block), "\u001f", as.character(d$label))
  d$key <- factor(k, levels = unique(k))
  ggplot(d, aes(estimate, key)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, alpha = credible),
                   height = 0, linewidth = 0.55, colour = fill_col) +
    geom_point(aes(alpha = credible, shape = credible),
               size = 1.9, colour = fill_col, fill = fill_col) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
    scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 1), guide = "none") +
    facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y",
               labeller = labeller(block = function(x)
                 str_wrap(str_remove(x, " \\(primary baseline\\)"), 9))) +
    scale_y_discrete(labels = function(x) sub("^[0-9]+\u001f", "", x)) +
    labs(x = xlab, y = NULL, title = title) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.y = element_text(size = 5.2),
          axis.text.x = element_text(size = 5.2, angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(size = 6),
          strip.placement = "outside",
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text.y.left = element_text(angle = 90, size = 4.6, face = "bold",
                                           margin = margin(1, 1, 1, 1)),
          panel.spacing.y = unit(1.2, "pt"),
          plot.title = element_text(size = 6.4, face = "bold"))
}

# Assemble a contrast frame from the two post-hoc tables.
read_posthoc <- function() {
  list(t1 = read.csv("results/cooney/posthoc_contrasts.csv"),
       t2 = read.csv("results/cooney/posthoc_interactions.csv"),
       sl = read.csv("results/cooney/colour_simple_slopes.csv"),
       cp = read.csv("results/cooney/a2d_change_from_primary.csv"))
}
mk <- function(d, block, label) {
  d$block <- block; d$label <- label
  d[, c("estimate", "lo", "hi", "credible", "block", "label")]
}

raw_theme <- theme(panel.grid.minor = element_blank(),
                   axis.text.x = element_text(size = 6.4, angle = 30,
                                              hjust = 1, vjust = 1),
                   plot.title = element_text(size = 8.5, face = "bold"),
                   legend.key.size = unit(3.2, "mm"),
                   legend.title = element_text(size = 6.5),
                   legend.text = element_text(size = 6))

# ============================================================
# FIGURE 2 - A3b community colour
# a/d  every coefficient, main effects then biome x land use
# b/e  observed response by land use           (main effects)
# c/f  observed response by land use and biome (interaction)
# ============================================================
message("Figure 6 ...")
try_fig({
  ph <- read_posthoc()

  # Panel a/c: the Biome4 x Predominant_simple interaction, shown as the
  # land-use contrast against primary vegetation computed WITHIN each
  # biome. The interaction is the difference between the three blocks.
  # Body mass is the one main effect retained.
  lu_in_biome <- ph$t2 |>
    filter(model == "A3b", quantity == "Level mean", factor_compared == "Land use",
           str_detect(held_fixed, "^Biome = "),
           str_detect(contrast, " - Primary vegetation$")) |>
    mutate(response = str_remove(response, "^Community ") |>
             str_replace("^male colourfulness$", "Male colourfulness") |>
             str_replace("^sexual dichromatism$", "Sexual dichromatism"),
           label = str_remove(contrast, " - Primary vegetation$"),
           block = factor(str_remove(held_fixed, "^Biome = "), levels = biome_levels))

  mass_c <- read.csv("results/cooney/A3b_fixed_effects.csv") |>
    filter(model %in% c("malecolcooney_covariate", "dichrodiff_covariate"),
           parameter == "z_logMass") |>
    transmute(response = c(malecolcooney_covariate = "Male colourfulness",
                           dichrodiff_covariate = "Sexual dichromatism")[model],
              estimate = Estimate, lo = Q2.5, hi = Q97.5,
              credible = sign(Q2.5) == sign(Q97.5),
              label = "Body mass (+1 SD)", block = "Mass")
  stopifnot(nrow(mass_c) == 2)

  a1b <- bind_rows(
      lu_in_biome |> select(response, estimate, lo, hi, credible, label, block),
      mass_c      |> select(response, estimate, lo, hi, credible, label, block)) |>
    mutate(response = factor(response, levels = resp_levels),
           credible = as.character(credible),
           block = factor(block, levels = c(biome_levels, "Mass")),
           label = factor(label, levels = rev(c(lu_levels[-1], "Body mass (+1 SD)"))))
  levels(a1b$block) <- c(paste0(biome_levels, "\n(primary baseline)"), "Slope")

  # Panel b/d: the same 15 cells in the observed data.
  raw <- read.csv("data/present.csv") |>
    select(Predominant_simple, Biome4, malecolcooney, dichrodiff) |>
    pivot_longer(c(malecolcooney, dichrodiff), names_to = "resp_raw", values_to = "value") |>
    filter(!is.na(value)) |>
    mutate(Predominant_simple = factor(relabel(Predominant_simple), levels = lu_levels),
           Biome4 = factor(relabel(Biome4), levels = biome_levels),
           response = factor(ifelse(resp_raw == "malecolcooney",
                                    "Male colourfulness", "Sexual dichromatism"),
                             levels = resp_levels))
  grp <- interaction(raw$Predominant_simple, raw$Biome4, drop = TRUE)
  ylim_dich <- whisker_lim(raw$value[raw$response == "Sexual dichromatism"],
                           grp[raw$response == "Sexual dichromatism"])
  ylim_mal  <- whisker_lim(raw$value[raw$response == "Male colourfulness"],
                           grp[raw$response == "Male colourfulness"])
  pad <- function(r, f = 0.04) r + c(-1, 1) * diff(r) * f
  ylim_dich <- pad(ylim_dich); ylim_mal <- pad(ylim_mal)

  panel_row <- function(rsp, tags, show_x) {
    col  <- resp_pal[[rsp]]
    logy <- rsp == "Male colourfulness"
    ylab <- paste0(rsp, if (logy) " (log scale)" else "")

    # Raw panel first, land use and biome as the two factors, boxplots
    # so the 15 cells stay legible and comparable by quartile.
    pa <- ggplot(filter(raw, response == rsp),
                 aes(Predominant_simple, value, fill = Biome4)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.28, colour = "grey25",
                   position = position_dodge(width = 0.8), width = 0.72) +
      scale_fill_manual(values = biome_pal, labels = biome_display, name = "Biome") +
      scale_x_discrete(labels = lu_short) +
      labs(x = if (show_x) "Land use" else NULL, y = ylab,
           title = paste0("(", tags[1], ") Observed, by land use and biome")) +
      raw_theme + theme(legend.position = "bottom",
                        legend.direction = "horizontal",
                        legend.justification = "center",
                        legend.title = element_blank(),
                        legend.text = element_text(size = 5.6),
                        legend.key.size = unit(2.6, "mm"),
                        legend.margin = margin(0, 0, 0, 0))
    pa <- if (logy) pa + scale_y_continuous(trans = "log10") +
                      coord_cartesian(ylim = ylim_mal)
          else pa + coord_cartesian(ylim = ylim_dich)

    pb <- forest_panel(filter(a1b, response == rsp), col,
                       paste0("(", tags[2], ") Contrasts"),
                       if (show_x) "Contrast (95% CrI)" else NULL)
    list(pa, pb)
  }

  r1 <- panel_row("Male colourfulness", c("a", "b"), FALSE)
  r2 <- panel_row("Sexual dichromatism", c("c", "d"), TRUE)

  # Each row is 80 mm tall, so a forest about 20 mm wide runs at roughly
  # 4:1. A tall narrow forest keeps every interval within one eye
  # movement of its label, and the estimates here span barely 0.1 units,
  # so horizontal room buys nothing.
  # wrap_plots per row, not `|` with a trailing plot_layout: the `/`
  # operator absorbs the outer plot_layout and the widths never reach
  # the nested rows, so the forest stays as wide as the boxplot.
  row1 <- wrap_plots(r1[[1]], r1[[2]], nrow = 1, widths = c(1, 0.20))
  row2 <- wrap_plots(r2[[1]], r2[[2]], nrow = 1, widths = c(1, 0.20))
  fig6 <- wrap_plots(row1, row2, ncol = 1)
  save_fig(fig6, "Figure_6_A3b_community_colour", mm2in(150), mm2in(160))
})

# ============================================================
# FIGURE 3 - A2d abundance change
# The response is always abundance change against a species' own
# primary-vegetation baseline. Colour is a PREDICTOR here, so rows
# are named for the model rather than for the response.
# a/d  every coefficient, main effects then interactions
# b/e  observed change by land use                 (main effects)
# c/f  observed change by biome and colour tercile (interaction)
# ============================================================
message("Figure 5 ...")
try_fig({
  ph <- read_posthoc()
  lu_d <- c("Secondary", "Plantation", "Cropland", "Pasture")
  cp_lab <- c("Male colourfulness" = "fitted with male colourfulness (z)",
              "Sexual dichromatism" = "fitted with sexual dichromatism (z)")

  t2 <- ph$t2 |> filter(model == "A2d") |>
    mutate(response = names(cp_lab)[match(source_fit, cp_lab)])

  # All three interaction terms of A2d, as pairwise level comparisons.
  sl_bi <- t2 |> filter(quantity == "Colour slope", factor_compared == "Biome") |>
    mutate(label = relabel(sub("^(.*) - (.*)$", "\\1 (\\2 baseline)", contrast)),
           block = "Colour slope\nbetween biomes")
  sl_lu <- t2 |> filter(quantity == "Colour slope", factor_compared == "Land use") |>
    mutate(label = relabel(sub("^(.*) - (.*)$", "\\1 (\\2 baseline)", contrast)),
           block = "Colour slope\nbetween land uses")
  # A2d carries no primary-vegetation level: primary is the baseline of
  # the response, so the raw panel's y axis already reads "vs primary".
  # The cell mean per land use would be the change from primary, but its
  # interval is ten times wider than any contrast (median width 0.098
  # against 0.0086) because it inherits the population intercept. Those
  # are kept in results/cooney/a2d_change_from_primary.csv instead.
  lu_bi <- t2 |> filter(quantity == "Level mean", factor_compared == "Land use",
                        str_detect(held_fixed, "^Biome = "),
                        str_detect(contrast, " - Secondary$")) |>
    mutate(label = relabel(str_remove(contrast, " - Secondary$")),
           block = paste0(relabel(str_remove(held_fixed, "^Biome = ")), "\n(secondary baseline)"))

  a2d <- bind_rows(sl_bi, sl_lu, lu_bi) |>
    select(response, estimate, lo, hi, credible, label, block) |>
    mutate(response = factor(response, levels = resp_levels),
           credible = as.character(credible),
           block = factor(block, levels = c("Colour slope\nbetween biomes",
                                            "Colour slope\nbetween land uses",
                                            paste0(biome_levels, "\n(secondary baseline)"))))
  levels(a2d$block) <- str_wrap(str_replace_all(levels(a2d$block), "\n", " "), 12)
  ord <- c(unique(sl_bi$label), unique(sl_lu$label), relabel(lu_d[-1]))
  a2d$label <- factor(a2d$label, levels = rev(unique(ord)))
  # The forests are narrowed to make room for the raw panels. Wrapping
  # these labels collides at 18 rows, so the baseline named inside the
  # brackets is abbreviated instead and each row stays on one line.
  abbrev <- c("Tropical closed baseline"  = "Tr. closed base",
              "Temperate closed baseline" = "Te. closed base",
              "Tropical open baseline"    = "Tr. open base",
              "Secondary baseline"        = "Secondary base",
              "Plantation baseline"       = "Plantation base",
              "Cropland baseline"         = "Cropland base")
  levels(a2d$label) <- str_replace_all(levels(a2d$label), abbrev)

  # ---- raw data for all three interaction terms -----------------------
  #
  # colour x biome      partial residuals, faceted by biome
  # colour x land use   partial residuals, faceted by land use
  # biome x land use    abundance change by both factors, colour ignored
  #
  # A plain plot of abundance change against colour shows nothing, and if
  # anything trends the wrong way (Spearman rho of colour decile against
  # mean change is -0.22, -0.12 and -0.01 by biome, against fitted slopes
  # of -0.0006, +0.0094 and +0.0079). That is not a binning problem. The
  # colour effect is a within-study, beyond-lineage quantity and does not
  # exist in the raw margin, which is dominated by which studies sampled
  # which species.
  #
  # The colour panels therefore plot a partial residual: abundance change
  # with the study mean and the species mean both removed. Colour is
  # constant within a species, so what remains is how far a species
  # departs from its own average within each biome or land use, which is
  # what the colour interactions estimate. Individual residuals span
  # roughly -0.08 to +0.14 while the effect is about 0.009 per standard
  # deviation, so the point cloud is twenty times the height of the
  # signal and is replaced by binned means with their intervals.
  #
  # The biome x land-use panel needs no residualising, because both
  # factors vary within a study, so it shows the observed values.
  #
  # DI1_2004__Naidoo 1 is excluded throughout this figure. It supplies
  # 94% of the tropical open cropland cell, has a median of 25 species
  # per site against 50 elsewhere, and an interquartile range of 0.115
  # against 0.012 for the other study in that cell. Because the response
  # is a share of the site total, species-poor sites inflate it
  # mechanically, so this one study would set the vertical scale.
  load("fits/A2d_model_setup.RData")
  DROP_SS <- "DI1_2004__Naidoo 1"
  base <- dat_final |>
    filter(!is.na(malecolcooney), SS != DROP_SS) |>
    mutate(Biome4 = factor(relabel(Biome4), levels = biome_levels),
           Predominant_simple = factor(relabel(Predominant_simple), levels = lu_d)) |>
    filter(!is.na(Biome4), !is.na(Predominant_simple))

  pr <- base |>
    group_by(SS)  |> mutate(r = diff_abund - mean(diff_abund)) |> ungroup() |>
    group_by(Best_guess_binomial) |> mutate(r = r - mean(r)) |> ungroup() |>
    select(Biome4, Predominant_simple, r,
           `Male colourfulness` = malecolcooney, `Sexual dichromatism` = dichrodiff) |>
    pivot_longer(c(`Male colourfulness`, `Sexual dichromatism`),
                 names_to = "response", values_to = "x") |>
    filter(!is.na(x)) |>
    mutate(response = factor(response, levels = resp_levels))

  NBIN <- 12
  resid_panel <- function(rsp, facet_by, tag, title, note = FALSE) {
    d  <- filter(pr, response == rsp)
    xr <- quantile(d$x, c(0.01, 0.99))
    d  <- filter(d, x >= xr[1], x <= xr[2])
    bins <- d |> group_by(.data[[facet_by]]) |> mutate(b = ntile(x, NBIN)) |>
      group_by(.data[[facet_by]], b) |>
      summarise(x = mean(x), m = mean(r), se = sd(r) / sqrt(n()), .groups = "drop")
    # A note in the first facet only, since the binning is the one thing
    # about these panels a reader cannot infer from the axes.
    note_df <- tibble(f = factor(levels(droplevels(d[[facet_by]]))[1],
                                 levels = levels(droplevels(d[[facet_by]]))),
                      # short lines: the first facet of panel a is only
                      # about 20 mm wide
                      lab = paste0("point: mean of\n1/12 of records\n",
                                   "bar: 95% CI"))
    names(note_df)[1] <- facet_by

    p <- ggplot(bins, aes(x, m)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
      geom_smooth(data = d, aes(x, r), method = "lm", formula = y ~ x, se = TRUE,
                  colour = "grey25", fill = "grey80", linewidth = 0.5, alpha = 0.5) +
      geom_errorbar(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se),
                    width = 0, linewidth = 0.35, colour = resp_pal[[rsp]]) +
      geom_point(size = 1.3, shape = 21, fill = resp_pal[[rsp]],
                 colour = "grey15", stroke = 0.3) +
      facet_wrap(as.formula(paste("~", facet_by)), nrow = 1, scales = "free_x",
                 labeller = labeller(.default = function(x) str_wrap(x, 9))) +
      labs(x = paste0(rsp, " (LociUVS)"), y = "Abundance change",
           title = paste0("(", tag, ") ", title)) +
      raw_theme +
      theme(strip.text = element_text(size = 6.6),
            axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 1))
    if (note)
      p <- p + geom_text(data = note_df, aes(x = -Inf, y = -Inf, label = lab),
                         hjust = -0.06, vjust = -0.25, size = 1.6, colour = "grey25",
                         lineheight = 1.0, inherit.aes = FALSE)
    p
  }

  # The forests run the full height of the figure in a narrow right-hand
  # column. A broad forest is harder to read, since the eye has to travel
  # a long way from label to interval, and the estimates here span barely
  # 0.03 units.
  mk_forest <- function(rsp, tag) {
    short <- rsp
    forest_panel(filter(a2d, response == rsp), resp_pal[[rsp]],
                 paste0("(", tag, ") Effects"),
                 "Effect (95% CrI)") +
      scale_x_continuous(breaks = c(-0.01, 0, 0.01)) +
      theme(axis.text.y = element_text(size = 5.0),
            strip.text.y.left = element_text(angle = 0, size = 4.8, face = "bold"),
            axis.text.x = element_text(size = 5.0, angle = 90, hjust = 1, vjust = 0.5),
            axis.title.x = element_text(size = 5.8))
  }

  pa <- resid_panel("Male colourfulness", "Biome4", "a",
                    "Male colourfulness by biome: partial residuals", note = TRUE)
  pb <- resid_panel("Male colourfulness", "Predominant_simple", "b",
                    "Male colourfulness by land use: partial residuals")
  pc <- resid_panel("Sexual dichromatism", "Biome4", "c",
                    "Sexual dichromatism by biome: partial residuals")
  pd <- resid_panel("Sexual dichromatism", "Predominant_simple", "d",
                    "Sexual dichromatism by land use: partial residuals")
  pf <- mk_forest("Male colourfulness", "f")
  pg <- mk_forest("Sexual dichromatism", "g")

  # biome x land use, observed values, colour ignored
  ylim_g <- whisker_lim(base$diff_abund,
                        interaction(base$Predominant_simple, base$Biome4,
                                    drop = TRUE))
  ylim_g <- ylim_g + c(-1, 1) * diff(ylim_g) * 0.04
  pe <- ggplot(base, aes(Predominant_simple, diff_abund, fill = Biome4)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.26, colour = "grey25",
                 position = position_dodge(width = 0.8), width = 0.74) +
    scale_fill_manual(values = biome_pal, name = "Biome") +
    scale_x_discrete(labels = lu_short) +
    coord_cartesian(ylim = ylim_g) +
    labs(x = "Land use", y = "Abundance change",
         title = "(e) Biome by land use: observed values") +
    raw_theme +
    theme(legend.position = "right",
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.6),
          legend.key.size = unit(2.8, "mm"))

  # Nested rather than one design grid. A design grid sizes rows to the
  # content that spans them, so the tall forests pulled the c/d band
  # taller than the a/b band even with equal row counts. Nesting makes
  # the two scatter rows a single two-row patchwork, which patchwork
  # splits evenly, so a, b, c and d are equal by construction. The same
  # holds for f and g.
  # wrap_plots throughout: the `|` and `/` operators let a nested
  # plot_layout be absorbed by the outer composition, which silently
  # drops the widths given to the row holding panel e.
  row_ab <- wrap_plots(pa, pb, nrow = 1, widths = c(12, 14))
  row_cd <- wrap_plots(pc, pd, nrow = 1, widths = c(12, 14))
  scat   <- wrap_plots(row_ab, row_cd, ncol = 1)
  e_row  <- wrap_plots(plot_spacer(), pe, plot_spacer(),
                       nrow = 1, widths = c(0.29, 0.42, 0.29))
  left   <- wrap_plots(scat, e_row, ncol = 1, heights = c(2, 0.85))
  right  <- wrap_plots(pf, pg, ncol = 1)
  # Two forests stacked in a 184 mm column are 92 mm each, so about
  # 23 mm wide holds them near 4:1.
  fig5   <- wrap_plots(left, right, nrow = 1, widths = c(205, 23)) &
    theme(plot.title = element_text(size = 7.0, face = "bold"))
  save_fig(fig5, "Figure_5_A2d_abundance_change", mm2in(232), mm2in(184))
})

# ============================================================
# FIGURE 4 - A1 species level
# a/d  radial phylogeny with ancestral state reconstruction
# b/e  species colour against land-use association score
# c/f  every pairwise land-use contrast, from the posterior
#
# Base graphics throughout, since phytools cannot render into a ggplot.
#
# Two decisions in panels b and e need stating. Species scoring
# exactly zero for a land use are omitted, because a zero says only
# that the species was never recorded there and 66% of species score
# zero for pasture, which buries the informative points under a wall
# at the axis. The fitted lines are the MODEL's predictions, drawn
# from the phylogenetic posterior, not least squares through the
# cloud. They are solid across the range of scores actually observed
# and dotted beyond it, because no species reaches a score of 1 for
# any land use except secondary vegetation, so a coefficient read at
# a score of 1 is an extrapolation.
# ============================================================
message("Figure 3 (contMap; ~1-2 min) ...")
try_fig({
  load("fits/A1_model_setup.RData")   # spdat, tree
  load("fits/A1_all_fits.RData")      # fits_A1

  prop_cols <- c("prop_Primary_vegetation", "prop_Secondary",
                 "prop_Plantation_forest", "prop_Cropland", "prop_Pasture")
  stopifnot(all(prop_cols %in% names(spdat)))

  # Colour scales are clipped to a robust range. Both metrics have long
  # tails which would otherwise compress every branch into the middle
  # of the ramp.
  build_map <- function(v, logscale, pal, lims) {
    x <- setNames(spdat[[v]], spdat$phylo); x <- x[!is.na(x)]
    tr <- drop.tip(tree, setdiff(tree$tip.label, names(x)))
    xx <- if (logscale) log(x) else x
    xx <- pmin(pmax(xx, lims[1]), lims[2])
    obj <- setMap(contMap(tr, xx, plot = FALSE, lims = lims), pal)
    list(obj = obj, lims = if (logscale) round(exp(lims)) else round(lims))
  }
  lim_mal <- unname(log(quantile(spdat$malecolcooney, c(0.02, 0.98), na.rm = TRUE)))
  # 2%/98% spans -119 to +119, but 82% of species sit within +/-30, so
  # that range parks almost every branch at the midpoint of the ramp.
  # The 10%/90% magnitude puts the colour where the species actually are
  # and clips only the tails, which were already clipped.
  q_dich  <- unname(max(abs(quantile(spdat$dichrodiff, c(0.10, 0.90), na.rm = TRUE))))
  # Mako, not magma. Zissou 1 in panel d is warm throughout, so a warm
  # sequential ramp here would let the two trees read as one scale. Mako
  # is cool from end to end (hues 94 to 300) and shares no hue with it.
  m_mal  <- build_map("malecolcooney", TRUE, hcl.colors(20, "Mako"), lim_mal)
  # Zissou 1 rather than a white-centred diverging ramp. Its midpoint sits
  # at lightness 82 against 97 for Blue-Red 3, and its whole range spans
  # only L 52 to 82, so a near-zero branch stays visible on white instead
  # of disappearing into the page. It also shares no hue with panel a.
  m_dich <- build_map("dichrodiff", FALSE, hcl.colors(20, "Zissou 1"), c(-q_dich, q_dich))

  contr <- read.csv("results/cooney/posthoc_contrasts.csv") |> filter(model == "A1")

  # Model-implied line for land use j. The five proportions sum to 1,
  # so raising p_j to t forces the rest down. We redistribute the
  # remaining 1-t in proportion to the observed mean composition of
  # the other four land uses, giving
  #     E[y | p_j = t] = t * b_j + (1 - t) * sum_k w_k b_k
  # a straight line in t. The band is its 95% credible interval.
  model_line <- function(draws, j, wbar, ts) {
    bj   <- draws[, paste0("b_", prop_cols[j])]
    w    <- wbar[-j] / sum(wbar[-j])
    rest <- as.vector(draws[, paste0("b_", prop_cols[-j])] %*% w)
    vapply(ts, function(t) {
      v <- t * bj + (1 - t) * rest
      c(median(v), quantile(v, 0.025), quantile(v, 0.975))
    }, numeric(3))
  }

  scatter_panel <- function(resp, ylab, logscale, tag, title, legend) {
    draws <- posterior::as_draws_matrix(fits_A1[[resp]])
    d  <- spdat[!is.na(spdat[[resp]]), ]
    wbar <- colMeans(d[, prop_cols])
    yv <- if (logscale) log(d[[resp]]) else d[[resp]]

    # Bins and lines are built before the device is opened, so the axis
    # can be set from what is actually drawn. The species themselves span
    # 32 to 560 LociUVS while the decile means span barely 70 to 130, and
    # a limit taken from the raw values would leave the panel almost empty.
    NBIN <- 10
    ns <- integer(0); B <- list(); L <- list()
    ts <- seq(0, 1, length.out = 80)
    for (j in seq_along(prop_cols)) {
      l <- lu_levels[j]; p <- d[[prop_cols[j]]]
      k <- p > 0                       # omit the uninformative zeros
      ns[l] <- sum(k)
      pk <- p[k]; yk <- yv[k]
      if (length(pk) >= 2 * NBIN) {
        br <- unique(quantile(pk, seq(0, 1, length.out = NBIN + 1)))
        b  <- cut(pk, breaks = br, include.lowest = TRUE, labels = FALSE)
        B[[l]] <- list(x  = tapply(pk, b, mean),
                       y  = tapply(yk, b, mean),
                       se = tapply(yk, b, function(z) sd(z) / sqrt(length(z))))
      }
      L[[l]] <- list(ts = ts, fit = model_line(draws, j, wbar, ts),
                     pmax_obs = max(p))
    }
    span <- range(c(unlist(lapply(B, function(z) c(z$y - 1.96 * z$se,
                                                   z$y + 1.96 * z$se))),
                    unlist(lapply(L, function(z) z$fit[1, ]))), na.rm = TRUE)
    ylim <- span + c(-1, 1) * diff(span) * 0.10

    plot(NA, xlim = c(0, 1), ylim = ylim,
         xlab = "Land-use association score (proportion of a species' records)",
         ylab = ylab, las = 1, cex.lab = 0.85, cex.axis = 0.8,
         yaxt = if (logscale) "n" else "s")
    if (logscale) {
      at <- pretty(exp(ylim), 5); at <- at[at > 0 &
                                           log(at) >= ylim[1] & log(at) <= ylim[2]]
      axis(2, at = log(at), labels = at, las = 1, cex.axis = 0.8)
    }
    # Decile means rather than the raw cloud, matching Figure 5.
    for (l in names(B)) {
      z <- B[[l]]
      segments(z$x, z$y - 1.96 * z$se, z$x, z$y + 1.96 * z$se,
               col = lu_pal[[l]], lwd = 0.9)
      points(z$x, z$y, pch = 21, cex = 0.72, lwd = 0.4,
             bg = lu_pal[[l]], col = "grey20")
    }
    for (l in names(L)) {
      z <- L[[l]]; inr <- z$ts <= z$pmax_obs
      # No credible band: its width is set by uncertainty in the shared
      # level, which is common to all five lines and cancels in the
      # contrasts. Panel c carries the uncertainty that matters.
      lines(z$ts[inr], z$fit[1, inr], col = lu_pal[[l]], lwd = 2.6)
      lines(z$ts[!inr], z$fit[1, !inr], col = lu_pal[[l]], lwd = 1.5, lty = 3)
    }
    mtext(paste0("(", tag, ") ", title), side = 3, adj = 0, line = 0.5,
          font = 2, cex = 0.8)
    if (legend) {
      legend("topright", legend = sprintf("%s (n = %d)", lu_display[lu_levels], ns[lu_levels]),
             col = lu_pal[lu_levels], lwd = 2.6, bty = "n", cex = 0.6, seg.len = 1.4)
      legend("topleft", legend = c("point: decile mean", "bar: 95% CI"),
             bty = "n", cex = 0.55, text.col = "grey25")
    }
  }

  forest_base <- function(resp_lab, tag, xlab) {
    d <- contr[contr$response == resp_lab, ]
    d <- d[rev(seq_len(nrow(d))), ]
    n <- nrow(d); ys <- seq_len(n)
    xr <- range(c(d$lo, d$hi, 0)); pad <- diff(xr) * 0.06
    par(mar = c(6.0, 4.2, 2.4, 2.2), mgp = c(2.6, 0.5, 0))
    plot(NA, xlim = c(xr[1] - pad, xr[2] + pad), ylim = c(0.5, n + 0.5),
         yaxt = "n", xaxt = "n", xlab = "", ylab = "", las = 1)
    axis(1, cex.axis = 0.62, las = 2, tcl = -0.25)
    mtext(xlab, side = 1, line = 4.6, adj = 0.5, cex = 0.48)
    abline(v = 0, lty = 2, col = "grey60")
    # 24 mm of column leaves no room for "Pasture - Plantation forest",
    # so each land use is cut to its shortest unambiguous stem.
    abbr <- c("Primary vegetation" = "Prim", "Plantation forest" = "Plant",
              "Secondary" = "Sec", "Cropland" = "Crop", "Pasture" = "Past")
    lab <- d$contrast
    for (nm in names(abbr)) lab <- gsub(nm, abbr[[nm]], lab, fixed = TRUE)
    lab <- gsub(" - ", "-", lab, fixed = TRUE)
    axis(2, at = ys, labels = lab, las = 1, cex.axis = 0.52, tick = FALSE, line = -0.6)
    for (i in ys) {
      cr <- d$credible[i]
      cl <- if (cr) "#B02020" else "grey55"
      segments(d$lo[i], i, d$hi[i], i, col = cl, lwd = if (cr) 2.1 else 1.3)
      points(d$estimate[i], i, pch = if (cr) 19 else 1,
             col = cl, cex = if (cr) 1.05 else 0.85)
    }
    mtext(paste0("(", tag, ") Contrasts"), side = 3, adj = 0.5,
          line = 0.5, font = 2, cex = 0.6)
  }

  tree_panel <- function(m, title, barlab) {
    tree_mar <- c(2.8, 0.5, 2.4, 0.5)
    par(mar = tree_mar)
    # plot.contMap resets par(mar) internally, so pass the margins in
    # or the panel label lands outside the device.
    plot(m$obj, type = "fan", legend = FALSE, fsize = 1e-4, lwd = 1.1, mar = tree_mar)
    mtext(title, side = 3, adj = 0, line = 0.7, font = 2, cex = 0.8)
    u <- par("usr")
    # title = "" and a separate mtext, because add.color.bar draws its
    # title at the bar's own y and the two collide.
    add.color.bar(diff(u[1:2]) * 0.5, m$obj$cols, title = "", lims = m$lims,
                  digits = 0, x = u[1] + diff(u[1:2]) * 0.25,
                  y = u[3] + diff(u[3:4]) * 0.045,
                  lwd = 8, fsize = 0.6, prompt = FALSE, subtitle = "")
    mtext(barlab, side = 1, line = 1.2, cex = 0.6)
  }

  draw3 <- function() {
    # Measured from the rendered page rather than guessed: at a 28 mm
    # column the plotted box came out 7 mm by 58 mm, or 8.4:1, because
    # the label margin eats most of a narrow column. 36 mm puts the box
    # near 15 mm by 58 mm, which is the 4:1 asked for.
    layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(95, 106, 36))

    tree_panel(m_mal, "(a) Male colourfulness across the phylogeny", "Male colourfulness (LociUVS)")
    par(mar = c(4.2, 4.6, 2.4, 1.2))
    scatter_panel("malecolcooney", "Male colourfulness (LociUVS)", TRUE, "b",
                  "Male colourfulness vs land-use association", TRUE)
    forest_base("Male colourfulness", "c", "log difference")

    tree_panel(m_dich, "(d) Sexual dichromatism across the phylogeny",
               "Sexual dichromatism (LociUVS)")
    par(mar = c(4.2, 4.6, 2.4, 1.2))
    scatter_panel("dichrodiff", "Sexual dichromatism (LociUVS)", FALSE, "e",
                  "Sexual dichromatism vs land-use association", FALSE)
    forest_base("Sexual dichromatism", "f", "LociUVS diff.")
  }

  pdf(file.path(OUT, "Figure_3_A1_phylogeny_landuse.pdf"),
      width = mm2in(241), height = mm2in(190)); draw3(); dev.off()
  png(file.path(OUT, "Figure_3_A1_phylogeny_landuse.png"),
      width = mm2in(241), height = mm2in(190), units = "in", res = 400); draw3(); dev.off()
  message("Saved Figure_3_A1_phylogeny_landuse")
})

message("\nDone. Publication figures in ", OUT)
