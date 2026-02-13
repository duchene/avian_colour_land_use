# ============================================================
# A1_03_model_comparison.R
# Analysis 1: Compare models using LOO-CV
# ============================================================

library(tidyverse)
library(brms)
library(loo)

# Load fits
load("fits/A1_model_setup.RData")
load("fits/A1_all_fits.RData")

cat("Loaded", length(fits_reduced), "reduced models and",
    length(fits_minimal), "minimal models\n")

# ============================================================
# LOO-CV COMPARISON
# ============================================================

loo_comparisons <- list()

for (resp in responses) {
  cat("\n", rep("=", 50), "\n")
  cat("LOO comparison for:", resp, "\n")
  cat(rep("=", 50), "\n")

  loo_red <- loo(fits_reduced[[resp]])
  loo_min <- loo(fits_minimal[[resp]])

  comp <- loo_compare(loo_red, loo_min)

  cat("\nReduced model LOO:\n")
  print(loo_red)

  cat("\nMinimal model LOO:\n")
  print(loo_min)

  cat("\nComparison (positive elpd_diff favors first model):\n")
  print(comp)

  loo_comparisons[[resp]] <- list(
    reduced = loo_red,
    minimal = loo_min,
    comparison = comp
  )
}

# ============================================================
# SUMMARY TABLE
# ============================================================

comparison_summary <- map_dfr(responses, function(resp) {
  comp <- loo_comparisons[[resp]]$comparison
  tibble(
    response = resp,
    best_model = rownames(comp)[1],
    elpd_diff = comp[2, "elpd_diff"],
    se_diff = comp[2, "se_diff"]
  )
})

cat("\n", rep("=", 50), "\n")
cat("SUMMARY: Model Comparison Results\n")
cat(rep("=", 50), "\n")
print(comparison_summary)

# ============================================================
# SAVE RESULTS
# ============================================================

dir.create("results", showWarnings = FALSE)
save(loo_comparisons, comparison_summary, file = "fits/A1_loo_comparisons.RData")
write_csv(comparison_summary, "results/A1_model_comparison_summary.csv")

cat("\nResults saved to results/\n")
