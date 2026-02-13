# ============================================================
# A1_04_diagnostics.R
# Analysis 1: Convergence diagnostics and posterior predictive checks
# ============================================================

library(tidyverse)
library(brms)
library(bayesplot)
library(posterior)

# Load fits
load("fits/A1_model_setup.RData")
load("fits/A1_all_fits.RData")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ============================================================
# CONVERGENCE DIAGNOSTICS FUNCTION
# ============================================================

check_convergence <- function(fit, model_name = "") {
  cat("\n=== Convergence Diagnostics:", model_name, "===\n")

  fit_summary <- summary(fit)$fixed
  rhat_vals <- fit_summary$Rhat
  ess_bulk <- fit_summary$Bulk_ESS
  ess_tail <- fit_summary$Tail_ESS

  # Random effects
  re_summary <- summary(fit)$random
  if (length(re_summary) > 0) {
    for (re in re_summary) {
      rhat_vals <- c(rhat_vals, re$Rhat)
      ess_bulk <- c(ess_bulk, re$Bulk_ESS)
      ess_tail <- c(ess_tail, re$Tail_ESS)
    }
  }

  # Divergences
  np <- nuts_params(fit)
  n_divergent <- sum(np$Value[np$Parameter == "divergent__"])
  max_td <- max(np$Value[np$Parameter == "treedepth__"])

  cat("Max Rhat:", round(max(rhat_vals, na.rm = TRUE), 4), "(target: < 1.01)\n")
  cat("Min ESS bulk:", round(min(ess_bulk, na.rm = TRUE), 0), "(target: > 400)\n")
  cat("Min ESS tail:", round(min(ess_tail, na.rm = TRUE), 0), "(target: > 400)\n")
  cat("Divergent transitions:", n_divergent, "(target: 0)\n")
  cat("Max treedepth:", max_td, "\n")

  converged <- max(rhat_vals, na.rm = TRUE) < 1.01 && n_divergent == 0
  cat("Converged:", ifelse(converged, "YES", "NO - check model"), "\n")

  tibble(
    model = model_name,
    max_rhat = max(rhat_vals, na.rm = TRUE),
    min_ess_bulk = min(ess_bulk, na.rm = TRUE),
    min_ess_tail = min(ess_tail, na.rm = TRUE),
    n_divergent = n_divergent,
    converged = converged
  )
}

# ============================================================
# RUN DIAGNOSTICS
# ============================================================

diagnostics <- list()

for (resp in responses) {
  diagnostics[[paste0(resp, "_reduced")]] <-
    check_convergence(fits_reduced[[resp]], paste(resp, "- reduced"))

  diagnostics[[paste0(resp, "_minimal")]] <-
    check_convergence(fits_minimal[[resp]], paste(resp, "- minimal"))
}

diagnostics_df <- bind_rows(diagnostics)

cat("\n", rep("=", 60), "\n")
cat("CONVERGENCE SUMMARY\n")
cat(rep("=", 60), "\n")
print(diagnostics_df)

write_csv(diagnostics_df, "results/A1_convergence_diagnostics.csv")

# ============================================================
# POSTERIOR PREDICTIVE CHECKS
# ============================================================

cat("\n", rep("=", 60), "\n")
cat("GENERATING POSTERIOR PREDICTIVE PLOTS\n")
cat(rep("=", 60), "\n")

for (resp in responses) {
  cat("\nPP check for:", resp, "\n")

  p <- pp_check(fits_reduced[[resp]], ndraws = 100) +
    ggtitle(paste("Posterior Predictive Check:", resp))

  ggsave(paste0("figures/A1_pp_check_", resp, ".png"), p,
         width = 8, height = 5, dpi = 150)
}

cat("\nDiagnostics complete. Results in results/ and figures/\n")
