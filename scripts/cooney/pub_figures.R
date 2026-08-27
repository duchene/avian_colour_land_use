# ============================================================
# pub_figures.R
# Publication figures for the avian colour x land-use manuscript.
#
# Figure 1 is a conceptual diagram, made separately.
# Figure 2  A3b community colour
# Figure 3  A2d abundance change
# Figure 4  A1 species level
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
    facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_discrete(labels = function(x) sub("^[0-9]+\u001f", "", x)) +
    labs(x = xlab, y = NULL, title = title) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 6.4),
          strip.placement = "outside",
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text.y.left = element_text(angle = 0, size = 6.2, face = "bold"),
          panel.spacing.y = unit(1.2, "pt"),
          plot.title = element_text(size = 8.5, face = "bold"))
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
message("Figure 2 ...")
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
  ylim_dich <- quantile(raw$value[raw$response == "Sexual dichromatism"],
                        c(0.005, 0.995), na.rm = TRUE)
  raw <- raw |> filter(response == "Male colourfulness" |
                         (value >= ylim_dich[1] & value <= ylim_dich[2]))

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
      raw_theme + theme(legend.position = "right")
    pa <- if (logy) pa + scale_y_continuous(trans = "log10")
          else pa + coord_cartesian(ylim = ylim_dich)

    pb <- forest_panel(filter(a1b, response == rsp), col,
                       paste0("(", tags[2], ") ", rsp, ": biome x land use"),
                       if (show_x) "Contrast or slope (95% CrI)" else NULL)
    list(pa, pb)
  }

  r1 <- panel_row("Male colourfulness", c("a", "b"), FALSE)
  r2 <- panel_row("Sexual dichromatism", c("c", "d"), TRUE)

  fig2 <- (r1[[1]] | r1[[2]]) / (r2[[1]] | r2[[2]]) +
    plot_layout(widths = c(1.0, 1.0))
  save_fig(fig2, "Figure_2_A3b_community_colour", mm2in(225), mm2in(160))
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
message("Figure 3 ...")
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
                 paste0("(", tag, ") ", short, " model"),
                 "Effect on abundance change (95% CrI)") +
      scale_x_continuous(breaks = c(-0.01, 0, 0.01)) +
      theme(axis.text.y = element_text(size = 6),
            strip.text.y.left = element_text(angle = 0, size = 5.6, face = "bold"),
            axis.text.x = element_text(size = 6),
            axis.title.x = element_text(size = 7))
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
  ylim_g <- quantile(base$diff_abund, c(0.06, 0.94))
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
  fig3   <- wrap_plots(left, right, nrow = 1, widths = c(2.75, 1)) &
    theme(plot.title = element_text(size = 7.6, face = "bold"))
  save_fig(fig3, "Figure_3_A2d_abundance_change", mm2in(280), mm2in(184))
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
message("Figure 4 (contMap; ~1-2 min) ...")
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
  q_dich  <- unname(max(abs(quantile(spdat$dichrodiff, c(0.02, 0.98), na.rm = TRUE))))
  m_mal  <- build_map("malecolcooney", TRUE, viridisLite::magma(20), lim_mal)
  # Diverging, since the difference is signed and centred near zero, and
  # deliberately far from magma so the two trees do not read as one scale.
  m_dich <- build_map("dichrodiff", FALSE, hcl.colors(20, "Tropic"), c(-q_dich, q_dich))

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

  scatter_panel <- function(resp, ylab, logscale, tag, title, ylim, legend) {
    draws <- posterior::as_draws_matrix(fits_A1[[resp]])
    d  <- spdat[!is.na(spdat[[resp]]), ]
    wbar <- colMeans(d[, prop_cols])
    yv <- if (logscale) log(d[[resp]]) else d[[resp]]
    plot(NA, xlim = c(0, 1), ylim = ylim,
         xlab = "Land-use association score (proportion of a species' records)",
         ylab = ylab, las = 1, cex.lab = 0.85, cex.axis = 0.8,
         yaxt = if (logscale) "n" else "s")
    if (logscale) {
      at <- c(30, 50, 100, 200, 400)
      axis(2, at = log(at), labels = at, las = 1, cex.axis = 0.8)
    }
    ns <- integer(0)
    for (j in seq_along(prop_cols)) {
      l <- lu_levels[j]; p <- d[[prop_cols[j]]]
      k <- p > 0                       # omit the uninformative zeros
      ns[l] <- sum(k)
      points(p[k], yv[k], pch = 16, cex = 0.62,
             col = adjustcolor(lu_pal[[l]], alpha.f = 0.28))
    }
    for (j in seq_along(prop_cols)) {
      l <- lu_levels[j]; p <- d[[prop_cols[j]]]
      pmax_obs <- max(p)
      ts  <- seq(0, 1, length.out = 80)
      fit <- model_line(draws, j, wbar, ts)
      inr <- ts <= pmax_obs
      # No credible band: its width is set by uncertainty in the shared
      # level, which is common to all five lines and cancels in the
      # contrasts. Panel c carries the uncertainty that matters.
      lines(ts[inr], fit[1, inr], col = lu_pal[[l]], lwd = 2.6)
      lines(ts[!inr], fit[1, !inr], col = lu_pal[[l]], lwd = 1.5, lty = 3)
    }
    mtext(paste0("(", tag, ") ", title), side = 3, adj = 0, line = 0.5,
          font = 2, cex = 0.8)
    if (legend)
      legend("topright", legend = sprintf("%s (n = %d)", lu_display[lu_levels], ns[lu_levels]),
             col = lu_pal[lu_levels], lwd = 2.6, bty = "n", cex = 0.6, seg.len = 1.4)
  }

  forest_base <- function(resp_lab, tag, xlab) {
    d <- contr[contr$response == resp_lab, ]
    d <- d[rev(seq_len(nrow(d))), ]
    n <- nrow(d); ys <- seq_len(n)
    xr <- range(c(d$lo, d$hi, 0)); pad <- diff(xr) * 0.06
    par(mar = c(4.2, 8.8, 2.4, 0.8))
    plot(NA, xlim = c(xr[1] - pad, xr[2] + pad), ylim = c(0.5, n + 0.5),
         yaxt = "n", xlab = xlab, ylab = "", las = 1,
         cex.lab = 0.85, cex.axis = 0.8)
    abline(v = 0, lty = 2, col = "grey60")
    lab <- gsub(" vegetation", "", gsub(" forest", "", d$contrast))
    axis(2, at = ys, labels = lab, las = 1, cex.axis = 0.66, tick = FALSE, line = -0.4)
    for (i in ys) {
      cr <- d$credible[i]
      cl <- if (cr) "#B02020" else "grey55"
      segments(d$lo[i], i, d$hi[i], i, col = cl, lwd = if (cr) 2.1 else 1.3)
      points(d$estimate[i], i, pch = if (cr) 19 else 1,
             col = cl, cex = if (cr) 1.05 else 0.85)
    }
    mtext(paste0("(", tag, ") Pairwise land-use contrasts"), side = 3, adj = 0,
          line = 0.5, font = 2, cex = 0.8)
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

  ylim_mal  <- log(c(28, 480))
  ylim_dich <- c(-75, 130)

  draw4 <- function() {
    layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(1.05, 1.35, 1.05))

    tree_panel(m_mal, "(a) Male colourfulness across the phylogeny", "Male colourfulness (LociUVS)")
    par(mar = c(4.2, 4.6, 2.4, 1.2))
    scatter_panel("malecolcooney", "Male colourfulness (LociUVS)", TRUE, "b",
                  "Male colourfulness vs land-use association", ylim_mal, TRUE)
    forest_base("Male colourfulness", "c", "Difference in log male colourfulness (95% CrI)")

    tree_panel(m_dich, "(d) Sexual dichromatism across the phylogeny",
               "Sexual dichromatism (LociUVS)")
    par(mar = c(4.2, 4.6, 2.4, 1.2))
    scatter_panel("dichrodiff", "Sexual dichromatism (LociUVS)", FALSE, "e",
                  "Sexual dichromatism vs land-use association", ylim_dich, FALSE)
    forest_base("Sexual dichromatism", "f", "Difference in LociUVS (95% CrI)")
  }

  pdf(file.path(OUT, "Figure_4_A1_phylogeny_landuse.pdf"),
      width = mm2in(265), height = mm2in(190)); draw4(); dev.off()
  png(file.path(OUT, "Figure_4_A1_phylogeny_landuse.png"),
      width = mm2in(265), height = mm2in(190), units = "in", res = 400); draw4(); dev.off()
  message("Saved Figure_4_A1_phylogeny_landuse")
})

message("\nDone. Publication figures in ", OUT)
