# ============================================================
# T1_tree_comparison.R
# Do the results depend on the phylogeny? Every phylogenetic quantity
# reported from the Claramunt et al. (2025) tree (results/cooney/) is set
# against the same quantity from the Jetz et al. (2012) MCC tree
# (runs/jetz/results/cooney/, written by hpc/run_jetz.sh).
#
# Quantities: A1b land-use contrasts and R2; A2d post-hoc contrasts,
# interaction contrasts and colour slopes; A2d per-colour slopes; A4
# MR-PMM colour on land-use slopes, correlations and phylogenetic signal;
# C1 sample-bias PGLS. A3a/A3b have no phylogeny and are not compared.
#
# same_conclusion: the 95% interval excludes zero under both trees with the
# same sign, or under neither.
# Outputs: results/cooney/T1_tree_comparison.csv (every quantity) and
#          results/cooney/T1_tree_comparison_summary.csv (per table)
# Run from the repo root.
# ============================================================

suppressPackageStartupMessages(library(dplyr))
A <- "results/cooney/"; B <- "runs/jetz/results/cooney/"

cmp <- function(file, keys, est, lo, hi, cred = NULL, filter = NULL) {
  a <- read.csv(paste0(A, file), stringsAsFactors = FALSE)
  b <- read.csv(paste0(B, file), stringsAsFactors = FALSE)
  if (!is.null(filter)) { a <- filter(a); b <- filter(b) }
  std <- function(d) {
    o <- data.frame(key = do.call(paste, c(d[keys], sep = " | ")),
                    est = d[[est]], lo = d[[lo]], hi = d[[hi]])
    o$cred <- if (is.null(cred)) sign(o$lo) == sign(o$hi) else as.logical(d[[cred]])
    o
  }
  m <- merge(std(a), std(b), by = "key", suffixes = c("_claramunt", "_jetz"))
  stopifnot(nrow(m) == nrow(a), nrow(m) == nrow(b))
  m$same_conclusion <- m$cred_claramunt == m$cred_jetz &
    (!m$cred_claramunt | sign(m$est_claramunt) == sign(m$est_jetz))
  cbind(table = sub("\\.csv$", "", file), m)
}
not_a3 <- function(d) d[d$model != "A3b" & d$model != "A3a", ]

res <- bind_rows(
  cmp("A1b_landuse_contrasts.csv", c("response", "contrast"), "Estimate", "Q2.5", "Q97.5", "credible"),
  cmp("A1b_r2_summary.csv", "response", "R2_estimate", "R2_lower", "R2_upper"),
  cmp("posthoc_contrasts.csv", c("model", "response", "factor", "contrast"), "estimate", "lo", "hi",
      "credible", filter = not_a3),
  cmp("posthoc_interactions.csv", c("model", "response", "source_fit", "quantity", "factor_compared",
                                    "held_fixed", "contrast"), "estimate", "lo", "hi", "credible",
      filter = not_a3),
  cmp("colour_simple_slopes.csv", c("model", "colour_predictor", "slope_within", "level"),
      "estimate", "lo", "hi", "credible", filter = not_a3),
  cmp("A2d_colours_slopes.csv", c("colour", "slope_within", "level"), "estimate", "lo", "hi", "credible"),
  cmp("A4_colour_landuse_slopes.csv", c("level", "colour", "land_use"), "sd_median", "sd_lower",
      "sd_upper", "excludes_zero"),
  cmp("A4_correlations.csv", c("level", "trait_1", "trait_2"), "median", "lower", "upper", "excludes_zero"),
  cmp("A4_phylogenetic_signal.csv", "trait", "median", "lower", "upper"))
# C1: PGLS difference with its p-value (no interval): significant at 0.05
c1 <- function(d) data.frame(key = paste(d$metric, d$sample, sep = " | "), est = d$pgls_diff,
                             lo = NA_real_, hi = NA_real_, cred = d$pgls_p < 0.05)
ca <- c1(read.csv(paste0(A, "C1_sample_bias.csv"))); cb <- c1(read.csv(paste0(B, "C1_sample_bias.csv")))
m <- merge(ca, cb, by = "key", suffixes = c("_claramunt", "_jetz"))
m$same_conclusion <- m$cred_claramunt == m$cred_jetz
res <- bind_rows(res, cbind(table = "C1_sample_bias", m))

write.csv(res, "results/cooney/T1_tree_comparison.csv", row.names = FALSE)
sm <- res %>% group_by(table) %>%
  summarise(n = n(), same_conclusion = sum(same_conclusion),
            credible_claramunt = sum(cred_claramunt), credible_jetz = sum(cred_jetz),
            cor_estimates = cor(est_claramunt, est_jetz),
            max_abs_diff = max(abs(est_claramunt - est_jetz)), .groups = "drop")
write.csv(sm, "results/cooney/T1_tree_comparison_summary.csv", row.names = FALSE)
print(as.data.frame(sm), digits = 3)
cat("\nQuantities whose conclusion differs between trees:\n")
print(res[!res$same_conclusion, c("table", "key", "est_claramunt", "lo_claramunt", "hi_claramunt",
                                  "est_jetz", "lo_jetz", "hi_jetz")], digits = 3, row.names = FALSE)
cat("Saved results/cooney/T1_tree_comparison{,_summary}.csv\n")
