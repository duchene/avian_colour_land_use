# pub_figures_dale.R
# Publication figures for the Dale-colour (B) replication of the manuscript.
# Parallels scripts/cooney/pub_figures.R exactly, using the B analyses and Dale
# plumage metrics: B1b (community colour), B1a (guild/size interactions), B2c
# (relative abundance), B2d (abundance change), B3 (phylogenetic association).
# Figure numbers match their Cooney counterparts one-for-one (1-7).
# Driven from results/dale/*.csv plus the B3 tree. Run: Rscript scripts/dale/pub_figures_dale.R

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr)
  library(patchwork); library(scales); library(ape); library(phytools)
})
theme_set(theme_bw(base_size = 11))
mm2in <- function(x) x / 25.4
save_fig <- function(p, name, w, h) {
  dir.create("figures/exploratory/dale", showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0("figures/exploratory/dale/", name, ".png"), p, width = w, height = h, dpi = 300)
  ggsave(paste0("figures/exploratory/dale/", name, ".pdf"), p, width = w, height = h)
  message("Saved ", name)
}
try_fig <- function(expr) tryCatch(expr, error = function(e) message("  !! FIGURE FAILED: ", conditionMessage(e)))

# ---- palettes / labels ------------------------------------------------------
lu_pal <- c("Primary vegetation"="#1B7837","Secondary"="#A6D96A",
            "Plantation forest"="#FEC44F","Cropland"="#D95F0E","Pasture"="#993404")
lu_nonref  <- c("Secondary","Plantation forest","Cropland","Pasture")
lu_full    <- c("Primary vegetation", lu_nonref)
biome_nonref <- c("Temperate Forest","Tropical Open")
resp_lab   <- c(meancoldale="Mean colour", malecoldale="Male colour",
                dichrodale="Dichromatism ratio", dichrodiffdale="Dichromatism difference")
resp_order <- unname(resp_lab)
cred_col   <- c("TRUE"="#0072B2","FALSE"="#B0B0B0")
cred_shape <- c("TRUE"=16,"FALSE"=1)
cred_labs  <- c("TRUE"="yes","FALSE"="no")

pretty_lu    <- function(x) str_replace(x, "Plantationforest", "Plantation forest")
pretty_biome <- function(x) x |>
  str_replace("TemperateForest","Temperate Forest") |>
  str_replace("TropicalOpen","Tropical Open") |>
  str_replace("TropicalForest","Tropical Forest")

# ============================================================
# Figure 1 — B1b community colour: main effects
# ============================================================
message("Figure 1 ...")
try_fig({
  b1b <- read.csv("results/dale/B1b_fixed_effects.csv") |>
    mutate(response = resp_lab[str_remove(model, "_(base|covariate)$")],
           variant  = str_extract(model, "base|covariate"),
           credible = sign(Q2.5) == sign(Q97.5))

  main1 <- b1b |>
    filter(variant == "covariate", parameter != "Intercept", !str_detect(parameter, ":"),
           str_detect(parameter, "^Predominant_simple|^Biome4|^z_logMass|^Trophic\\.Niche")) |>
    mutate(group = case_when(str_detect(parameter,"^Predominant_simple") ~ "Land use",
                             str_detect(parameter,"^Biome4")             ~ "Biome",
                             str_detect(parameter,"^z_logMass")          ~ "Body mass",
                             TRUE                                        ~ "Trophic niche"),
           label = parameter |>
             str_replace("^Predominant_simple","") |>
             str_replace("^Biome4","") |>
             str_replace("^Trophic\\.Niche","") |>
             str_replace("^z_logMass","Body mass (log)") |>
             pretty_lu() |> pretty_biome() |>
             str_replace("Herbivoreterrestrial","Herbivore terrestrial"),
           response = factor(response, levels = resp_order))
  ord <- main1 |> distinct(group, label) |>
    mutate(group = factor(group, levels = c("Land use","Biome","Body mass","Trophic niche"))) |>
    arrange(group, label) |> pull(label) |> unique()
  main1$label <- factor(main1$label, levels = rev(ord))

  fig1 <- ggplot(main1, aes(Estimate, label, colour = credible, shape = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.3, linewidth = 0.5) +
    geom_point(size = 2) +
    facet_wrap(~response, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = cred_col, name = "95% CI excludes 0", labels = cred_labs) +
    scale_shape_manual(values = cred_shape, name = "95% CI excludes 0", labels = cred_labs) +
    labs(x = "Coefficient (reference: Primary vegetation / Tropical Forest / Aquatic predator)",
         y = NULL, title = "B1b - community colour (Dale): land use, biome, body mass, trophic niche") +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 8))
  save_fig(fig1, "Figure_1_B1b_main_effects", mm2in(220), mm2in(150))
})

# ============================================================
# Figure 2 — B1b biome x land-use interactions
# ============================================================
message("Figure 2 ...")
try_fig({
  b1b <- read.csv("results/dale/B1b_fixed_effects.csv") |>
    mutate(response = resp_lab[str_remove(model, "_(base|covariate)$")],
           variant  = str_extract(model, "base|covariate"),
           credible = sign(Q2.5) == sign(Q97.5))
  int1 <- b1b |>
    filter(variant == "covariate", str_detect(parameter, "^Biome4.*:Predominant_simple"),
           !prior_only_empty_cell) |>
    mutate(biome    = str_remove(parameter, ":Predominant_simple.*$") |> str_remove("^Biome4") |> pretty_biome(),
           land_use = str_remove(parameter, "^Biome4.*:Predominant_simple") |> pretty_lu(),
           term = paste0(biome, " x ", land_use),
           response = factor(response, levels = resp_order))
  fig2 <- ggplot(int1, aes(Estimate, term, colour = credible, shape = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.25, linewidth = 0.5) +
    geom_point(size = 1.8) +
    facet_wrap(~response, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = cred_col, name = "95% CI excludes 0", labels = cred_labs) +
    scale_shape_manual(values = cred_shape, name = "95% CI excludes 0", labels = cred_labs) +
    labs(x = "Biome x land-use interaction coefficient", y = NULL,
         title = "B1b - biome x land-use interactions on community colour (Dale)") +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
  save_fig(fig2, "Figure_2_B1b_biome_landuse", mm2in(220), mm2in(150))
})

# ============================================================
# Figure 3 — B1a: guild- and size-specific colour responses to land use
#            (land use x trophic niche and land use x body mass). B1a is
#            favoured over B1b by LOO for all four responses.
# ============================================================
message("Figure 3 ...")
try_fig({
  b1a <- read.csv("results/dale/B1a_fixed_effects.csv") |>
    mutate(response = factor(resp_lab[model], levels = resp_order),
           credible = sign(Q2.5) == sign(Q97.5))

  # (a) body mass x land use
  massint <- b1a |>
    filter(str_detect(parameter, "^Predominant_simple[A-Za-z]+:z_logMass")) |>
    mutate(land_use = parameter |> str_extract("Predominant_simple[A-Za-z]+") |>
             str_remove("Predominant_simple") |> pretty_lu(),
           land_use = factor(land_use,
             levels = rev(c("Cropland","Pasture","Plantation forest","Secondary"))))
  p_mass <- ggplot(massint, aes(Estimate, land_use, colour = credible, shape = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.25, linewidth = 0.5) +
    geom_point(size = 2) +
    facet_wrap(~response, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = cred_col, name = "95% CI excludes 0", labels = cred_labs) +
    scale_shape_manual(values = cred_shape, name = "95% CI excludes 0", labels = cred_labs) +
    labs(x = "Body mass x land-use interaction (per SD log mass)", y = NULL,
         subtitle = "(a) Body mass x land use") +
    theme(legend.position = "none", axis.text.y = element_text(size = 8))

  # (b) trophic niche x land use (Herbivore terrestrial, N=2, excluded: prior-only)
  lu_ord    <- c("Cropland","Pasture","Plantation forest","Secondary")
  guild_ord <- c("Frugivore","Granivore","Invertivore","Nectarivore","Omnivore")
  ylev <- as.vector(t(outer(lu_ord, guild_ord,
                            function(a, b) paste0(a, ": ", b))))
  trophint <- b1a |>
    filter(str_detect(parameter, "^Predominant_simple[A-Za-z]+:Trophic\\.Niche"),
           !sparse_trophic) |>
    mutate(land_use = parameter |> str_extract("Predominant_simple[A-Za-z]+") |>
             str_remove("Predominant_simple") |> pretty_lu(),
           guild = parameter |> str_extract("Trophic\\.Niche[A-Za-z]+") |>
             str_remove("Trophic.Niche"),
           ylab = factor(paste0(land_use, ": ", guild), levels = rev(ylev)))
  p_troph <- ggplot(trophint, aes(Estimate, ylab, colour = credible, shape = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.25, linewidth = 0.45) +
    geom_point(size = 1.8) +
    facet_wrap(~response, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = cred_col, name = "95% CI excludes 0", labels = cred_labs) +
    scale_shape_manual(values = cred_shape, name = "95% CI excludes 0", labels = cred_labs) +
    labs(x = "Trophic niche x land-use interaction", y = NULL,
         subtitle = "(b) Trophic niche x land use") +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 7))

  fig3 <- (p_mass / p_troph) +
    plot_layout(heights = c(1, 4)) +
    plot_annotation(
      title = "B1a - guild- and size-specific colour responses to land use (Dale)")
  save_fig(fig3, "Figure_3_B1a_trophic_mass", mm2in(220), mm2in(200))
})

# ============================================================
# Figure 4 — Colour effects on abundance (B2c relative abundance, B2d change)
# ============================================================
message("Figure 4 ...")
try_fig({
  term_levels <- c("main effect","x Temperate Forest","x Tropical Open",
                   "x Cropland","x Pasture","x Plantation forest","x Secondary")
  read_ct <- function(path, lab) read.csv(path) |> mutate(analysis = lab)
  c4 <- bind_rows(read_ct("results/dale/B2c_fixed_effects.csv","B2c (relative abundance)"),
                  read_ct("results/dale/B2d_fixed_effects.csv","B2d (abundance change)")) |>
    filter(str_detect(parameter, "^z_")) |>
    mutate(metric = ifelse(str_detect(parameter,"malecoldale"),"Male colour","Dichromatism difference"),
           term = case_when(
             parameter %in% c("z_malecoldale","z_dichrodiffdale") ~ "main effect",
             str_detect(parameter,":Biome4") ~ paste0("x ", pretty_biome(str_remove(parameter,"^z_[a-z]+:Biome4"))),
             TRUE ~ paste0("x ", pretty_lu(str_remove(parameter,"^z_[a-z]+:Predominant_simple")))),
           credible = sign(Q2.5) == sign(Q97.5),
           term = factor(term, levels = rev(term_levels)),
           metric = factor(metric, levels = c("Male colour","Dichromatism difference")))
  fig5 <- ggplot(c4, aes(Estimate, term, colour = credible, shape = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.25, linewidth = 0.5) +
    geom_point(size = 2) +
    facet_grid(metric ~ analysis, scales = "free_x") +
    scale_colour_manual(values = cred_col, name = "95% CI excludes 0", labels = cred_labs) +
    scale_shape_manual(values = cred_shape, name = "95% CI excludes 0", labels = cred_labs) +
    labs(x = "Colour effect on abundance / abundance change", y = NULL,
         title = "Dale colour effects on abundance: credible terms persist in B2c") +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 8))
  save_fig(fig5, "Figure_4_B2c_B2d_colour_terms", mm2in(200), mm2in(140))
})

# ============================================================
# Figure 5 — Abundance variance components (phylogeny/study dominate)
# ============================================================
message("Figure 5 ...")
try_fig({
  vc <- bind_rows(read.csv("results/dale/B2c_variance_components.csv") |> mutate(analysis="B2c (relative abundance)"),
                  read.csv("results/dale/B2d_variance_components.csv") |> mutate(analysis="B2d (abundance change)")) |>
    filter(group != "residual__") |>
    mutate(component = recode(group, phylo="Phylogeny", SS="Study", SSB="Block",
                              SSBS="Site", sigma="Residual"),
           metric = ifelse(model == "z_malecoldale","Male colour","Dichromatism difference"),
           component = factor(component, levels = c("Phylogeny","Study","Block","Site","Residual")))
  fig6 <- ggplot(vc, aes(component, sd, fill = component)) +
    geom_col() +
    facet_grid(metric ~ analysis, scales = "free_y") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = NULL, y = "Random-effect / residual SD",
         title = "Abundance variance is dominated by phylogeny and study, not colour (Dale)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(fig6, "Figure_5_variance_components", mm2in(180), mm2in(130))
})

# ============================================================
# Figure 6 — B3 radial phylogenies with ancestral colour states (Dale)
# ============================================================
message("Figure 6 (contMap; ~1-2 min) ...")
try_fig({
  load("fits/B3_model_setup.RData")   # tree, spdat
  x_mal <- setNames(spdat$malecoldale, spdat$phylo); x_mal <- x_mal[!is.na(x_mal)]
  tree_mal <- drop.tip(tree, setdiff(tree$tip.label, names(x_mal)))
  obj_a <- contMap(tree_mal, log(x_mal), plot = FALSE)
  obj_a <- setMap(obj_a, viridisLite::magma(20))
  lims_mal <- round(exp(obj_a$lims), 1)

  x_dich <- setNames(spdat$dichrodiffdale, spdat$phylo); x_dich <- x_dich[!is.na(x_dich)]
  tree_dich <- drop.tip(tree, setdiff(tree$tip.label, names(x_dich)))
  lim_d <- max(abs(x_dich))
  obj_b <- contMap(tree_dich, x_dich, plot = FALSE, lims = c(-lim_d, lim_d))
  obj_b <- setMap(obj_b, hcl.colors(20, "Cork"))

  draw6 <- function() {
    layout(matrix(c(1,2,3,4), nrow = 2, byrow = TRUE), heights = c(5,1))
    par(mar = c(0,1,2,1))
    plot(obj_a, type = "fan", legend = FALSE, fsize = 1e-4, lwd = 1.4, main = "(a) Male colour")
    par(mar = c(0,1,2,1))
    plot(obj_b, type = "fan", legend = FALSE, fsize = 1e-4, lwd = 1.4, main = "(b) Dichromatism difference")
    par(mar = c(3,2,0,2)); plot.new()
    add.color.bar(0.9, obj_a$cols, title = "Male colour (Dale)", lims = lims_mal,
                  digits = 1, x = 0.05, y = 0.6, lwd = 10, fsize = 0.9, prompt = FALSE, subtitle = "")
    par(mar = c(3,2,0,2)); plot.new()
    add.color.bar(0.9, obj_b$cols, title = "Dichromatism difference (Dale)", lims = obj_b$lims,
                  digits = 1, x = 0.05, y = 0.6, lwd = 10, fsize = 0.9, prompt = FALSE, subtitle = "")
  }
  dir.create("figures/exploratory/dale", showWarnings = FALSE, recursive = TRUE)
  pdf("figures/exploratory/dale/Figure_6_B3_phylogeny.pdf", width = mm2in(200), height = mm2in(130)); draw6(); dev.off()
  png("figures/exploratory/dale/Figure_6_B3_phylogeny.png", width = mm2in(200), height = mm2in(130), units = "in", res = 300); draw6(); dev.off()
  message("Saved Figure_6_B3_phylogeny")
})

# ============================================================
# Figure 7 — B3 expected colour by land-use association
# ============================================================
message("Figure 7 ...")
try_fig({
  b3 <- read.csv("results/dale/B3_fixed_effects.csv") |>
    mutate(land_use = str_remove(parameter, "^prop_") |> str_replace_all("_", " "),
           land_use = factor(land_use, levels = rev(lu_full)),
           response = factor(resp_lab[response], levels = resp_order))
  fig8 <- ggplot(b3, aes(Estimate, land_use, colour = land_use)) +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.2, linewidth = 0.7) +
    geom_point(size = 2.6) +
    facet_wrap(~response, scales = "free_x") +
    scale_colour_manual(values = lu_pal, guide = "none") +
    labs(x = "Expected colour for a species found exclusively in that land use", y = NULL,
         title = "B3 - phylogenetic association of colour with land use (Dale)") +
    theme(axis.text.y = element_text(size = 8))
  save_fig(fig8, "Figure_7_B3_landuse_association", mm2in(200), mm2in(140))
})

message("\nDone. Figures in figures/exploratory/dale/")
