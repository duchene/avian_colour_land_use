# ============================================================
# supp_fig_S1_tree.R
# Figure S1: the PREDICTS passerines across the passerine phylogeny.
#
# A circular BBtree2 phylogeny of every passerine with Delhey et al. 2023
# colour data, surrounded by rings (inner to outer):
#   1  species present in the PREDICTS data (dark tick)
#   2  male colourfulness (Cooney et al. 2022, LociUVS, log scale)
#   3  sexual dichromatism (Cooney et al. 2022 plate score)
#   4-13  proportion of male plumage in each of the ten Delhey colours,
#         each ring shaded from white (none) to that colour, capped at the
#         97.5th percentile of the colour across passerines
# Grey cells mark species without that measurement.
# The statistical comparison of these distributions is in Figure S2
# (C1_sample_bias.R).
#
# Output: figures/pub/cooney/Figure_S1_passerine_tree.{png,pdf}
# ============================================================

source("scripts/cooney/delhey_colours.R")
suppressPackageStartupMessages(library(ape))

tree <- read.tree("data/BBtree2.tre")
dl   <- delhey_male()
dl   <- dl[dl$phylo %in% tree$tip.label, ]
tr   <- drop.tip(tree, setdiff(tree$tip.label, dl$phylo))
tr   <- ladderize(tr)
tips <- tr$tip.label

raw  <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"),
                 stringsAsFactors = FALSE)
ours <- unique(c(raw$jetz_sp, gsub(" ", "_", raw$Best_guess_binomial)))
co   <- read.csv("data/external/cooney2022_passerines.csv", stringsAsFactors = FALSE)
dlm  <- dl[match(tips, dl$phylo), ]
com  <- co[match(tips, co$binomial), ]
cat("Figure S1:", length(tips), "passerine tips,", sum(tips %in% ours), "in the PREDICTS data\n")

ramp <- function(x, cols, lim) {
  pal <- colorRampPalette(cols)(100)
  i <- pmin(pmax(round(99 * (x - lim[1]) / diff(lim)) + 1, 1), 100)
  out <- pal[i]; out[is.na(x)] <- "grey85"; out
}
hue <- c(black = "#111111", white = "#9a9a9a", grey = "#6e6e6e", brown = "#7b4a1e",
         rufous = "#c0561b", yellow = "#e8c000", red = "#d01c1c", green = "#2e8b2e",
         blue = "#1f5fd1", purple = "#7b2fa8")
lc   <- log(com$LociUVS_male)
rings <- c(list(
  list(lab = "In PREDICTS data", col = ifelse(tips %in% ours, "#1b1b1b", "white")),
  list(lab = "Male colourfulness", col = ramp(lc, hcl.colors(9, "Mako"),
                                              quantile(lc, c(0.025, 0.975), na.rm = TRUE))),
  list(lab = "Sexual dichromatism", col = ramp(com$Dichro, hcl.colors(9, "Zissou 1"), c(-1, 1)))),
  lapply(DELHEY_COLOURS, function(cc)
    list(lab = COLOUR_LABEL[[cc]],
         col = ramp(dlm[[cc]], c("white", hue[[cc]]),
                    c(0, quantile(dl[[cc]], 0.975))))))

draw <- function() {
  par(mar = c(0, 0, 0, 0))
  w <- 0.035; gap <- 0.006; R0 <- 1.02
  H   <- max(node.depth.edgelength(tr))   # fan coordinates are in branch-length units
  lim <- H * (R0 + length(rings) * (w + gap) + 0.02)
  plot(tr, type = "fan", show.tip.label = FALSE, edge.width = 0.15, edge.color = "grey30",
       x.lim = c(-lim, lim), y.lim = c(-lim, lim), no.margin = TRUE)
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  xx <- pp$xx[seq_along(tips)]; yy <- pp$yy[seq_along(tips)]
  rmax <- max(sqrt(xx^2 + yy^2))
  # plot.phylo scales the tree to the plot; express ring radii relative to it
  sc <- rmax
  th <- atan2(yy, xx)
  dth <- 2 * pi / length(tips) / 2
  for (k in seq_along(rings)) {
    r1 <- sc * (R0 + (k - 1) * (w + gap)); r2 <- r1 + sc * w
    cl <- rings[[k]]$col
    for (i in seq_along(tips)) {
      a <- c(th[i] - dth, th[i] + dth)
      polygon(c(r1 * cos(a), rev(r2 * cos(a))), c(r1 * sin(a), rev(r2 * sin(a))),
              col = cl[i], border = NA)
    }
    # ring label at the gap just above the 3 o'clock position
    text(r1 + sc * w / 2, 0, labels = k, cex = 0.35, col = "black", srt = 0)
  }
  legend("bottomleft", inset = 0.005, bty = "n", cex = 0.55, title = "Rings (inner to outer)",
         legend = paste(seq_along(rings), vapply(rings, `[[`, "", "lab")), title.adj = 0)
  # Colour keys for rings 2 and 3 (bottom right); colour rings 4-13 run
  # from white (none) to the full colour at the 97.5th percentile.
  key <- function(y0, cols, lims, labs, title) {
    x0 <- 0.66 * lim; x1 <- 0.93 * lim; h <- 0.025 * lim
    pal <- colorRampPalette(cols)(100); xs <- seq(x0, x1, length.out = 101)
    rect(xs[-101], y0, xs[-1], y0 + h, col = pal, border = NA)
    rect(x0, y0, x1, y0 + h, border = "grey30", lwd = 0.4)
    text(c(x0, x1), y0 - 0.5 * h, labs, cex = 0.45, adj = c(0.5, 1))
    text(x0, y0 + 1.6 * h, title, cex = 0.5, adj = c(0, 0.5))
  }
  ql <- quantile(lc, c(0.025, 0.975), na.rm = TRUE)
  key(-0.80 * lim, hcl.colors(9, "Mako"), ql, sprintf("%.0f", exp(ql)),
      "2 Male colourfulness (LociUVS, log scale)")
  key(-0.93 * lim, hcl.colors(9, "Zissou 1"), c(-1, 1),
      c("-1 or below (female)", "1 or above (male)"), "3 Sexual dichromatism (plate score)")
}
out <- "figures/pub/cooney"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
png(file.path(out, "Figure_S1_passerine_tree.png"), width = 9, height = 9, units = "in", res = 400)
draw(); dev.off()
pdf(file.path(out, "Figure_S1_passerine_tree.pdf"), width = 9, height = 9)
draw(); dev.off()
cat("Saved", file.path(out, "Figure_S1_passerine_tree.{png,pdf}"), "\n")
