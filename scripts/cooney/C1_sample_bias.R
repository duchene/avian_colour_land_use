# ============================================================
# C1_sample_bias.R
# Are the PREDICTS passerines a biased sample of passerine colour?
#
# Compares the colour of the species in our data with all passerines
# measured by the source datasets:
#   male colourfulness and plate-scored sexual dichromatism: Cooney et
#     al. 2022 (4,527 passerines)
#   proportion of male plumage in each of ten colours: Delhey et al.
#     2023 (5,447 passerines), each colour separately (no ratios)
# for two samples: every passerine in the PREDICTS extract, and the
# species entering the species-level model (A1b).
#
# Two tests per metric and sample:
#   1. Distributions: medians, means and a two-sample Kolmogorov-Smirnov
#      test, sampled species versus the remaining passerines.
#   2. Phylogenetic regression (PGLS, Pagel's lambda, phylolm) on all
#      passerines on the BBtree2 tips: metric ~ in_sample. The coefficient
#      is the colour difference of sampled species after accounting for
#      the fact that sampled species cluster in clades.
# Colourfulness is log-transformed before the PGLS; dichromatism and the
# colour proportions are used as they are.
#
# Outputs: results/cooney/C1_sample_bias.csv and the supplementary
# figure figures/pub/cooney/Figure_S2_sample_bias.{png,pdf}.
# Run from the repo root after A1b_01_setup.R:
#   Rscript scripts/cooney/C1_sample_bias.R
# ============================================================

source("scripts/00_config.R")
source("scripts/cooney/delhey_colours.R")
suppressPackageStartupMessages({ library(ape); library(phylolm); library(ggplot2) })

tree <- read_phylo()   # 2026-11: PHYLO selects the tree (00_config.R)
bb_tips <- read.tree("data/BBtree2.tre")$tip.label
raw  <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"),
                 stringsAsFactors = FALSE)
ours <- unique(raw[, c("Best_guess_binomial", "jetz_sp")])
ours_names <- unique(c(ours$jetz_sp, gsub(" ", "_", ours$Best_guess_binomial)))
# Load into its own environment: the setup file also holds `tree` (pruned
# to the A1b species), which must not overwrite the full tree above.
a1b <- new.env(); load("fits/A1b_model_setup.RData", envir = a1b)   # spdat: the A1b species
a1b_names <- unique(c(a1b$spdat$jetz_sp, gsub(" ", "_", a1b$spdat$Best_guess_binomial)))

co <- read.csv("data/external/cooney2022_passerines.csv", stringsAsFactors = FALSE)
dl <- delhey_male()

ref <- c(list(
  list(metric = "Male colourfulness", key = co$binomial, x = co$LociUVS_male,
       tx = function(v) log(v), scale = "log", source = "Cooney et al. 2022"),
  list(metric = "Sexual dichromatism", key = co$binomial, x = co$Dichro,
       tx = identity, scale = "raw", source = "Cooney et al. 2022")),
  lapply(DELHEY_COLOURS, function(cc)
    list(metric = COLOUR_LABEL[[cc]], key = dl$phylo, x = dl[[cc]], tx = identity,
         scale = "proportion", source = "Delhey et al. 2023")))
samples <- list("PREDICTS passerines" = ours_names, "Species-level model (A1b)" = a1b_names)

res <- list(); plot_df <- list()
for (r in ref) for (sname in names(samples)) {
  ok  <- !is.na(r$x)
  key <- r$key[ok]; x <- r$x[ok]
  ins <- key %in% samples[[sname]]
  kt  <- suppressWarnings(ks.test(x[ins], x[!ins]))
  # PGLS on the passerines that sit on BBtree2 tips
  # the same species under both trees: those on the BBtree2 tips
  on <- key %in% tree$tip.label & key %in% bb_tips
  dd <- data.frame(y = r$tx(x[on]), in_sample = as.numeric(ins[on]), row.names = key[on])
  tr <- drop.tip(tree, setdiff(tree$tip.label, rownames(dd)))
  pg <- phylolm(y ~ in_sample, data = dd, phy = tr, model = "lambda")
  cf <- summary(pg)$coefficients["in_sample", ]
  res[[length(res) + 1]] <- data.frame(
    metric = r$metric, source = r$source, sample = sname,
    n_sample = sum(ins), n_rest = sum(!ins),
    median_sample = median(x[ins]), median_rest = median(x[!ins]),
    mean_sample = mean(x[ins]), mean_rest = mean(x[!ins]),
    ks_D = unname(kt$statistic), ks_p = kt$p.value,
    pgls_n = nrow(dd), pgls_lambda = pg$optpar,
    pgls_diff = cf[["Estimate"]], pgls_se = cf[["StdErr"]], pgls_p = cf[["p.value"]],
    pgls_scale = r$scale)
  if (sname == "PREDICTS passerines")
    plot_df[[length(plot_df) + 1]] <- data.frame(metric = r$metric, value = x,
      group = ifelse(ins, "In PREDICTS data", "Other passerines"))
}
res <- do.call(rbind, res)
dir.create("results/cooney", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "results/cooney/C1_sample_bias.csv", row.names = FALSE)
print(res[, c("metric", "sample", "n_sample", "n_rest", "mean_sample", "mean_rest",
              "ks_D", "ks_p", "pgls_lambda", "pgls_diff", "pgls_se", "pgls_p")], digits = 3)

# ---------- Figure S2 ---------------------------------------------------
# Cumulative distributions (the quantity the Kolmogorov-Smirnov test
# compares) for the PREDICTS sample against the other passerines, with
# the phylogenetic (PGLS) difference written in each panel. Cumulative
# curves stay readable for the rare colours, where most species carry
# none of the colour. Colourfulness is on a log axis, its PGLS scale.
lv <- vapply(ref, `[[`, "", "metric")
pd <- do.call(rbind, plot_df); pd$metric <- factor(pd$metric, lv)
pd$value[pd$metric == "Male colourfulness"] <- log10(pd$value[pd$metric == "Male colourfulness"])
levels(pd$metric)[1] <- "Male colourfulness (log10)"
an <- res[res$sample == "PREDICTS passerines", ]
an$metric <- factor(an$metric, lv); levels(an$metric)[1] <- "Male colourfulness (log10)"
an$label <- sprintf("PGLS diff. %s\np %s", formatC(an$pgls_diff, digits = 2, format = "g"),
                    ifelse(an$pgls_p < 0.001, "< 0.001", paste("=", formatC(an$pgls_p, digits = 2, format = "f"))))
p <- ggplot(pd, aes(value, colour = group)) +
  stat_ecdf(linewidth = 0.45, pad = FALSE) +
  geom_text(data = an, aes(x = Inf, y = -Inf, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = -0.4, size = 2.3, lineheight = 0.9) +
  facet_wrap(~ metric, scales = "free_x", ncol = 4) +
  scale_colour_manual(values = c("In PREDICTS data" = "#2a6f97", "Other passerines" = "#999999")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x = "Value (proportion of male plumage for the ten colours)",
       y = "Cumulative proportion of species", colour = NULL) +
  theme_classic(base_size = 9) + theme(legend.position = "bottom",
                                       strip.background = element_blank())
out <- "figures/pub/cooney"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out, "Figure_S2_sample_bias.png"), p, width = 7.5, height = 6, dpi = 400)
ggsave(file.path(out, "Figure_S2_sample_bias.pdf"), p, width = 7.5, height = 6)
unlink(c("figures/cooney/C1_sample_bias.png", "figures/cooney/C1_sample_bias.pdf"))
cat("Saved results/cooney/C1_sample_bias.csv and", file.path(out, "Figure_S2_sample_bias.{png,pdf}"), "\n")
