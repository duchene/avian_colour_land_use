# pub_figures.R
# Publication figures for the avian colour x land-use manuscript.
# Focus: A1b (community colour), A2c (relative abundance), A2d (abundance
# change), A3 (phylogenetic association). Driven from results/*.csv (no need
# to reload the large brms fits) plus the A3 tree. Run: Rscript scripts/pub_figures.R

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr)
  library(patchwork); library(scales); library(ape); library(phytools)
})
theme_set(theme_bw(base_size = 11))
mm2in <- function(x) x / 25.4
save_fig <- function(p, name, w, h) {
  dir.create("figures/pub/cooney", showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0("figures/pub/cooney/", name, ".png"), p, width = w, height = h, dpi = 300)
  ggsave(paste0("figures/pub/cooney/", name, ".pdf"), p, width = w, height = h)
  message("Saved ", name)
}
try_fig <- function(expr) tryCatch(expr, error = function(e) message("  !! FIGURE FAILED: ", conditionMessage(e)))

# ---- palettes / labels ------------------------------------------------------
lu_pal <- c("Primary vegetation"="#1B7837","Secondary"="#A6D96A",
            "Plantation forest"="#FEC44F","Cropland"="#D95F0E","Pasture"="#993404")
lu_nonref  <- c("Secondary","Plantation forest","Cropland","Pasture")
lu_full    <- c("Primary vegetation", lu_nonref)
biome_nonref <- c("Temperate Forest","Temperate Open","Tropical Open")
resp_lab   <- c(meancolcooney="Mean colour", malecolcooney="Male colour",
                dichrocooney="Dichromatism ratio", dichrodiff="Dichromatism difference")
resp_order <- unname(resp_lab)
cred_col   <- c("TRUE"="#0072B2","FALSE"="#B0B0B0")
cred_shape <- c("TRUE"=16,"FALSE"=1)
cred_labs  <- c("TRUE"="yes","FALSE"="no")

pretty_lu    <- function(x) str_replace(x, "Plantationforest", "Plantation forest")
pretty_biome <- function(x) x |>
  str_replace("TemperateForest","Temperate Forest") |>
  str_replace("TemperateOpen","Temperate Open") |>
  str_replace("TropicalOpen","Tropical Open") |>
  str_replace("TropicalForest","Tropical Forest")

# ============================================================
# Figure 1 — A1b community colour: main effects
# ============================================================
message("Figure 1 ...")
try_fig({
  a1b <- read.csv("results/cooney/A1b_fixed_effects.csv") |>
    mutate(response = resp_lab[str_remove(model, "_(base|covariate)$")],
           variant  = str_extract(model, "base|covariate"),
           credible = sign(Q2.5) == sign(Q97.5))

  main1 <- a1b |>
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
         y = NULL, title = "A1b - community colour: land use, biome, body mass, trophic niche") +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 8))
  save_fig(fig1, "Figure_1_A1b_main_effects", mm2in(220), mm2in(150))
})

# ============================================================
# Figure 2 — A1b biome x land-use interactions
# ============================================================
message("Figure 2 ...")
try_fig({
  a1b <- read.csv("results/cooney/A1b_fixed_effects.csv") |>
    mutate(response = resp_lab[str_remove(model, "_(base|covariate)$")],
           variant  = str_extract(model, "base|covariate"),
           credible = sign(Q2.5) == sign(Q97.5))
  int1 <- a1b |>
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
         title = "A1b - biome x land-use interactions on community colour") +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
  save_fig(fig2, "Figure_2_A1b_biome_landuse", mm2in(220), mm2in(150))
})

# ============================================================
# Figure 3 — A1a: guild- and size-specific colour responses to land use
#            (land use x trophic niche and land use x body mass). A1a is
#            favoured over A1b by LOO for all four responses; these are the
#            interactions A1b's main-effects-only structure cannot capture.
# ============================================================
message("Figure 3 ...")
try_fig({
  a1a <- read.csv("results/cooney/A1a_fixed_effects.csv") |>
    mutate(response = factor(resp_lab[model], levels = resp_order),
           credible = sign(Q2.5) == sign(Q97.5))

  # (a) body mass x land use
  massint <- a1a |>
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
  trophint <- a1a |>
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
      title = "A1a - guild- and size-specific colour responses to land use")
  save_fig(fig3, "Figure_3_A1a_trophic_mass", mm2in(220), mm2in(200))
})

# ============================================================
# Figure 4 — Colour x land-use effects vanish once biome + phylogeny added
# ============================================================
message("Figure 4 ...")
try_fig({
  get_cxlu <- function(path, cpcol, label) {
    read.csv(path) |> rename(cp = all_of(cpcol)) |>
      filter(cp %in% c("z_malecolcooney","z_dichrodiff"),
             str_detect(parameter, "^z_[a-z]+:Predominant_simple")) |>
      mutate(land_use = str_remove(parameter, "^z_[a-z]+:Predominant_simple") |> pretty_lu(),
             metric   = ifelse(cp == "z_malecolcooney", "Male colour", "Dichromatism difference"),
             model    = label, credible = sign(Q2.5) == sign(Q97.5))
  }
  d3 <- bind_rows(
    get_cxlu("results/cooney/A2_fixed_effects.csv",  "colour_predictor", "A2 (colour x land use only)"),
    get_cxlu("results/cooney/A2c_fixed_effects.csv", "model",            "A2c (+ biome + phylogeny)")) |>
    mutate(land_use = factor(land_use, levels = rev(lu_nonref)),
           metric   = factor(metric, levels = c("Male colour","Dichromatism difference")),
           model    = factor(model, levels = c("A2 (colour x land use only)","A2c (+ biome + phylogeny)")))
  fig3 <- ggplot(d3, aes(Estimate, land_use, colour = model, shape = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.25, linewidth = 0.6,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.5)) +
    facet_wrap(~metric) +
    scale_colour_manual(values = c("A2 (colour x land use only)"="#D55E00",
                                   "A2c (+ biome + phylogeny)"="#0072B2"), name = NULL) +
    scale_shape_manual(values = cred_shape, name = "95% CI excludes 0", labels = cred_labs) +
    labs(x = "Colour x land-use interaction on relative abundance", y = NULL,
         title = "Colour x land-use effects disappear once biome and phylogeny are included") +
    theme(legend.position = "bottom")
  save_fig(fig3, "Figure_4_A2_vs_A2c_vanishing", mm2in(200), mm2in(110))
})

# ============================================================
# Figure 5 — Colour effects on abundance: null (A2c) / negligible (A2d)
# ============================================================
message("Figure 5 ...")
try_fig({
  term_levels <- c("main effect","x Temperate Forest","x Temperate Open","x Tropical Open",
                   "x Cropland","x Pasture","x Plantation forest","x Secondary")
  read_ct <- function(path, lab) read.csv(path) |> mutate(analysis = lab)
  c4 <- bind_rows(read_ct("results/cooney/A2c_fixed_effects.csv","A2c (relative abundance)"),
                  read_ct("results/cooney/A2d_fixed_effects.csv","A2d (abundance change)")) |>
    filter(str_detect(parameter, "^z_")) |>
    mutate(metric = ifelse(str_detect(parameter,"malecolcooney"),"Male colour","Dichromatism difference"),
           term = case_when(
             parameter %in% c("z_malecolcooney","z_dichrodiff") ~ "main effect",
             str_detect(parameter,":Biome4") ~ paste0("x ", pretty_biome(str_remove(parameter,"^z_[a-z]+:Biome4"))),
             TRUE ~ paste0("x ", pretty_lu(str_remove(parameter,"^z_[a-z]+:Predominant_simple")))),
           credible = sign(Q2.5) == sign(Q97.5),
           term = factor(term, levels = rev(term_levels)),
           metric = factor(metric, levels = c("Male colour","Dichromatism difference")))
  fig4 <- ggplot(c4, aes(Estimate, term, colour = credible, shape = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.25, linewidth = 0.5) +
    geom_point(size = 2) +
    facet_grid(metric ~ analysis, scales = "free_x") +
    scale_colour_manual(values = cred_col, name = "95% CI excludes 0", labels = cred_labs) +
    scale_shape_manual(values = cred_shape, name = "95% CI excludes 0", labels = cred_labs) +
    labs(x = "Colour effect on abundance / abundance change", y = NULL,
         title = "Colour effects on abundance are null (A2c) or negligible (A2d)") +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 8))
  save_fig(fig4, "Figure_5_A2c_A2d_colour_terms", mm2in(200), mm2in(140))
})

# ============================================================
# Figure 6 — Abundance variance components (phylogeny/study dominate)
# ============================================================
message("Figure 6 ...")
try_fig({
  vc <- bind_rows(read.csv("results/cooney/A2c_variance_components.csv") |> mutate(analysis="A2c (relative abundance)"),
                  read.csv("results/cooney/A2d_variance_components.csv") |> mutate(analysis="A2d (abundance change)")) |>
    filter(group != "residual__") |>
    mutate(component = recode(group, phylo="Phylogeny", SS="Study", SSB="Block",
                              SSBS="Site", sigma="Residual"),
           metric = ifelse(model == "z_malecolcooney","Male colour","Dichromatism difference"),
           component = factor(component, levels = c("Phylogeny","Study","Block","Site","Residual")))
  fig5 <- ggplot(vc, aes(component, sd, fill = component)) +
    geom_col() +
    facet_grid(metric ~ analysis, scales = "free_y") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = NULL, y = "Random-effect / residual SD",
         title = "Abundance variance is dominated by phylogeny and study, not colour") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(fig5, "Figure_6_variance_components", mm2in(180), mm2in(130))
})

# ============================================================
# Figure 7 — A3 radial phylogenies with ancestral colour states
# ============================================================
message("Figure 7 (contMap; ~1-2 min) ...")
try_fig({
  load("fits/A3_model_setup.RData")   # tree, spdat
  x_mal <- setNames(spdat$malecolcooney, spdat$phylo); x_mal <- x_mal[!is.na(x_mal)]
  tree_mal <- drop.tip(tree, setdiff(tree$tip.label, names(x_mal)))
  obj_a <- contMap(tree_mal, log(x_mal), plot = FALSE)
  obj_a <- setMap(obj_a, viridisLite::magma(20))
  lims_mal <- round(exp(obj_a$lims), 1)

  x_dich <- setNames(spdat$dichrodiff, spdat$phylo); x_dich <- x_dich[!is.na(x_dich)]
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
    add.color.bar(0.9, obj_a$cols, title = "Male colour (LociUVS)", lims = lims_mal,
                  digits = 1, x = 0.05, y = 0.6, lwd = 10, fsize = 0.9, prompt = FALSE, subtitle = "")
    par(mar = c(3,2,0,2)); plot.new()
    add.color.bar(0.9, obj_b$cols, title = "Dichromatism difference (LociUVS)", lims = obj_b$lims,
                  digits = 1, x = 0.05, y = 0.6, lwd = 10, fsize = 0.9, prompt = FALSE, subtitle = "")
  }
  dir.create("figures/pub/cooney", showWarnings = FALSE, recursive = TRUE)
  pdf("figures/pub/cooney/Figure_7_A3_phylogeny.pdf", width = mm2in(200), height = mm2in(130)); draw6(); dev.off()
  png("figures/pub/cooney/Figure_7_A3_phylogeny.png", width = mm2in(200), height = mm2in(130), units = "in", res = 300); draw6(); dev.off()
  message("Saved Figure_7_A3_phylogeny")
})

# ============================================================
# Figure 8 — A3 expected colour by land-use association
# ============================================================
message("Figure 8 ...")
try_fig({
  a3 <- read.csv("results/cooney/A3_fixed_effects.csv") |>
    mutate(land_use = str_remove(parameter, "^prop_") |> str_replace_all("_", " "),
           land_use = factor(land_use, levels = rev(lu_full)),
           response = factor(resp_lab[response], levels = resp_order))
  fig7 <- ggplot(a3, aes(Estimate, land_use, colour = land_use)) +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.2, linewidth = 0.7) +
    geom_point(size = 2.6) +
    facet_wrap(~response, scales = "free_x") +
    scale_colour_manual(values = lu_pal, guide = "none") +
    labs(x = "Expected colour for a species found exclusively in that land use", y = NULL,
         title = "A3 - phylogenetic association of colour with land use") +
    theme(axis.text.y = element_text(size = 8))
  save_fig(fig7, "Figure_8_A3_landuse_association", mm2in(200), mm2in(140))
})

message("\nDone. Figures in figures/pub/")
