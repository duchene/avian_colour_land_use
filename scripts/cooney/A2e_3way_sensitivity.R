# ============================================================
# A2e_3way_sensitivity.R
# Sensitivity check for A2c: FULL three-way interaction
#   relabund ~ z_malecolcooney * Biome4 * Predominant_simple
#              + (1|SS)+(1|SSB)+(1|SSBS) + (1|gr(phylo,cov=A))
# Compared against the reported 2-way A2c model via LOO.
# Reuses the A2c setup + sampler settings so the comparison is fair.
# ============================================================
source("scripts/00_config.R")
library(tidyverse); library(brms); library(cmdstanr); library(posterior)
options(brms.backend = "cmdstanr")
dir.create("results/cooney", showWarnings = FALSE)

load("fits/A2c_model_setup.RData")   # dat, A, colour_vars, formulas, priors_A2c, nthreads
cv <- "z_malecolcooney"

d <- dat %>% filter(!is.na(.data[[cv]]))
d$phylo <- droplevels(d$phylo)
tips <- levels(d$phylo); A_cv <- A[tips, tips]

form3 <- bf(relabund ~ z_malecolcooney * Biome4 * Predominant_simple +
              (1 | SS) + (1 | SSB) + (1 | SSBS) + (1 | gr(phylo, cov = A)))

f3 <- "fits/A2e_fit_z_malecolcooney.RData"
if (file.exists(f3)) { load(f3) } else {
  cat("=== fitting A2e (3-way) | n =", nrow(d), "===\n")
  fit3 <- brm(form3, data = d, data2 = list(A = A_cv), family = lognormal(),
              prior = priors_A2c, chains = 4, iter = 2000, warmup = 1000,
              cores = 4, threads = threading(nthreads),
              control = list(adapt_delta = 0.9, max_treedepth = 10),
              seed = 1, refresh = 200)
  save(fit3, file = f3)
}

# ---- LOO comparison vs reported 2-way A2c model ----
fit3 <- add_criterion(fit3, "loo")
save(fit3, file = f3)
load("fits/A2c_fit_z_malecolcooney.RData")     # `fit` = 2-way A2c
fit2 <- add_criterion(fit, "loo")
cmp <- loo_compare(fit2, fit3)
cat("\n=== LOO compare: fit2 = A2c 2-way, fit3 = A2e 3-way ===\n"); print(cmp)
write.csv(as.data.frame(cmp), "results/cooney/A2e_loo_compare.csv")

# ---- convergence ----
s  <- summary(fit3)
rh <- c(s$fixed[,"Rhat"], s$spec_pars[,"Rhat"], unlist(lapply(s$random, function(x) x[,"Rhat"])))
eb <- c(s$fixed[,"Bulk_ESS"], s$spec_pars[,"Bulk_ESS"], unlist(lapply(s$random, function(x) x[,"Bulk_ESS"])))
np <- nuts_params(fit3)
conv <- tibble(model = "A2e_3way_malecolcooney",
               max_rhat = round(max(rh, na.rm = TRUE), 4),
               min_ess_bulk = round(min(eb, na.rm = TRUE)),
               n_divergent = sum(np$Value[np$Parameter == "divergent__"]))
write_csv(conv, "results/cooney/A2e_convergence.csv"); print(conv)

# ---- fixed effects with empty-cell flag (any term touching an empty biome x LU cell) ----
emptns <- empty_cell_params(d)
fx <- as.data.frame(fixef(fit3)); fx$parameter <- rownames(fx)
fx$prior_only_empty_cell <- vapply(fx$parameter,
  function(p) any(vapply(emptns, function(e) grepl(e, p, fixed = TRUE), logical(1))), logical(1))
write_csv(fx %>% select(parameter, Estimate, Est.Error, Q2.5, Q97.5, prior_only_empty_cell),
          "results/cooney/A2e_fixed_effects.csv")

r2 <- bayes_R2(fit3)
write_csv(tibble(model = "A2e_3way_malecolcooney", R2 = r2[,"Estimate"], Q2.5 = r2[,"Q2.5"], Q97.5 = r2[,"Q97.5"]),
          "results/cooney/A2e_r2_summary.csv")
cat("\nA2e complete.\n")
