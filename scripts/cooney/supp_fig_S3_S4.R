# ============================================================
# supp_fig_S3_S4.R
# Supplementary Figures S3 and S4, from saved result tables only.
#
# Figure S3 (MR-PMM, A4_colour_mrpmm.R):
#   a  correlations among the eight colours and four land-use scores,
#      phylogenetic level below the diagonal, residual level above it;
#      cells whose 95% credible interval excludes zero are outlined
#   b  phylogenetic signal of each response, S_phy / (S_phy + S_res)
# Figure S4 (colour and land use, each colour a separate variable):
#   a  MR-PMM: conditional slope of each colour on each land-use score
#      at the total (species) level, in SD units, primary vegetation as
#      the implicit reference
#   b  abundance change (A2d fitted once per colour): slope of abundance
#      change on the standardised colour proportion within each land use
#
# Inputs:  results/cooney/A4_{correlations,phylogenetic_signal,
#          colour_landuse_slopes}.csv, results/cooney/A2d_colours_slopes.csv
# Outputs: figures/pub/cooney/Figure_S3_mrpmm_correlations.{png,pdf}
#          figures/pub/cooney/Figure_S4_colour_landuse.{png,pdf}
# ============================================================

source("scripts/cooney/delhey_colours.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

LUS  <- c("Secondary", "Plantation", "Cropland", "Pasture")
COLS <- setdiff(DELHEY_COLOURS, c("grey", "purple"))
LAB  <- c(COLOUR_LABEL[COLS], setNames(LUS, LUS))
ORD  <- unname(LAB)
hue  <- c(Black = "#111111", White = "#9a9a9a", Grey = "#6e6e6e", Brown = "#7b4a1e",
          Rufous = "#c0561b", Yellow = "#d4af00", Red = "#d01c1c", Green = "#2e8b2e",
          Blue = "#1f5fd1", Purple = "#7b2fa8")
out <- "figures/pub/cooney"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---------- Figure S3 ---------------------------------------------------
cr <- read.csv("results/cooney/A4_correlations.csv", stringsAsFactors = FALSE)
cr$t1 <- LAB[cr$trait_1]; cr$t2 <- LAB[cr$trait_2]
# Phylogenetic below the diagonal (row later in ORD than column), residual above
pos <- function(a, b) match(a, ORD) - match(b, ORD)
ph  <- cr[cr$level == "phylogenetic", ]; rs <- cr[cr$level == "residual", ]
sw  <- function(d, below) {
  flip <- (pos(d$t1, d$t2) > 0) != below
  d[flip, c("t1", "t2")] <- d[flip, c("t2", "t1")]; d
}
hm <- rbind(sw(ph, TRUE), sw(rs, FALSE))
hm$x <- factor(hm$t2, ORD); hm$y <- factor(hm$t1, rev(ORD))
hm$lab <- sub("^(-?)0\\.", "\\1.", sprintf("%.2f", hm$median))
pa <- ggplot(hm, aes(x, y, fill = median)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = hm[hm$excludes_zero, ], colour = "black", linewidth = 0.45, fill = NA) +
  geom_text(aes(label = lab), size = 1.9) +
  geom_abline(intercept = length(ORD) + 1, slope = -1, colour = "grey60", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1),
                       name = "Correlation") +
  annotate("text", x = 1, y = 0.2, label = "Below diagonal: phylogenetic", hjust = 0, size = 2.3) +
  annotate("text", x = length(ORD), y = length(ORD) + 0.8, label = "Above diagonal: residual",
           hjust = 1, size = 2.3) +
  coord_equal(clip = "off") + labs(x = NULL, y = NULL, title = "a") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))

ps <- read.csv("results/cooney/A4_phylogenetic_signal.csv", stringsAsFactors = FALSE)
ps$trait <- factor(LAB[ps$trait], rev(ORD))
ps$type <- ifelse(ps$trait %in% LUS, "Land-use score", "Colour")
pb <- ggplot(ps, aes(median, trait, colour = type)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.25, linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_colour_manual(values = c(Colour = "#2a6f97", "Land-use score" = "#c46a00"), name = NULL) +
  labs(x = "Phylogenetic signal", y = NULL, title = "b") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
s3 <- pa + pb + plot_layout(widths = c(3, 1))
ggsave(file.path(out, "Figure_S3_mrpmm_correlations.png"), s3, width = 9, height = 5.6, dpi = 400, bg = "white")
ggsave(file.path(out, "Figure_S3_mrpmm_correlations.pdf"), s3, width = 9, height = 5.6)

# ---------- Figure S4 ---------------------------------------------------
sl <- read.csv("results/cooney/A4_colour_landuse_slopes.csv", stringsAsFactors = FALSE)
sl <- sl[sl$level == "total", ]
sl$colour <- factor(COLOUR_LABEL[sl$colour], rev(COLOUR_LABEL[COLS]))
sl$land_use <- factor(sl$land_use, LUS)
# lim: optional symmetric axis limit; intervals beyond it are clipped and
# drawn with an arrow; estimates beyond it are not drawn.
forest <- function(d, x, lo, hi, xlab, title, lim = NULL) {
  d$.x <- d[[x]]; d$.lo <- d[[lo]]; d$.hi <- d[[hi]]
  arr <- NULL
  if (!is.null(lim)) {
    cl <- d[d$.lo < -lim | d$.hi > lim, ]
    arr <- rbind(
      transform(cl[cl$.lo < -lim, ], a = -lim * 0.8, b = -lim),
      transform(cl[cl$.hi > lim, ], a = lim * 0.8, b = lim))
    d$.lo <- pmax(d$.lo, -lim); d$.hi <- pmin(d$.hi, lim)
    d$.x[abs(d$.x) > lim] <- NA   # estimate off the axis: interval only
  }
  g <- ggplot(d, aes(y = colour, colour = colour)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
    geom_linerange(aes(xmin = .lo, xmax = .hi), linewidth = 0.5) +
    geom_point(data = d[!is.na(d$.x), ], aes(x = .x, shape = credible), size = 1.2,
               stroke = 0.8, fill = "white")
  if (!is.null(arr) && nrow(arr))
    g <- g + geom_segment(data = arr, aes(x = a, xend = b, y = colour, yend = colour),
                          arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.5)
  g +
    facet_wrap(~ land_use, nrow = 1) +
    scale_colour_manual(values = hue, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
    labs(x = xlab, y = NULL, title = title) +
    theme_classic(base_size = 8) +
    theme(strip.background = element_blank(), plot.title = element_text(face = "bold"))
}
sl$credible <- sl$excludes_zero
s4a <- forest(sl, "sd_median", "sd_lower", "sd_upper",
              "Conditional slope of colour on land-use score (SD units; reference: primary vegetation)", "a")

ab <- read.csv("results/cooney/A2d_colours_slopes.csv", stringsAsFactors = FALSE)
ab <- ab[ab$slope_within %in% c("Overall", "Land use"), ]
ab$land_use <- factor(ifelse(ab$slope_within == "Overall", "All land uses", ab$level),
                      c("All land uses", LUS))
ab$colour <- factor(ab$colour, rev(COLOUR_LABEL))
s4b <- forest(ab, "estimate", "lo", "hi",
              "Slope of abundance change on colour proportion (per SD of the colour; arrows: interval truncated)",
              "b", lim = 0.015)
# Green's land-use slopes average over biomes and inherit its temperate
# closed slope (-0.18, 95% CI -0.37 to 0.02; few green species), so its
# intervals run beyond the axis and are truncated.
s4 <- s4a / s4b
ggsave(file.path(out, "Figure_S4_colour_landuse.png"), s4, width = 8.5, height = 6.5, dpi = 400, bg = "white")
ggsave(file.path(out, "Figure_S4_colour_landuse.pdf"), s4, width = 8.5, height = 6.5)
cat("Saved", file.path(out, "Figure_S3_mrpmm_correlations.{png,pdf}"), "and",
    file.path(out, "Figure_S4_colour_landuse.{png,pdf}"), "\n")
