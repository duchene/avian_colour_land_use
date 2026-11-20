# Tiny fits of every (model, new metric) pair. Called by hpc/01_smoke.sh.
source("scripts/00_config.R")
suppressMessages({library(brms); library(cmdstanr); library(ape)})
quick <- function(...) brm(..., chains = 1, iter = 120, warmup = 60, refresh = 0,
                           seed = 1, silent = 2, backend = "cmdstanr")
ok <- function(tag, expr) {
  r <- tryCatch({ fit <- expr; sprintf("OK   n=%d", nobs(fit)) },
                error = function(e) paste("FAIL", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", tag, r))
}
# A1b: species level, phylogeny, 300-species subset
load("fits/A1b_model_setup.RData")
s <- spdat[spdat$phylo %in% sample(spdat$phylo, 300), ]
for (resp in "dichroplate") {
  d <- s[!is.na(s[[resp]]), ]; fp <- family_prior(resp, "species")
  ok(paste("A1b", resp), quick(bf(as.formula(paste0(resp, " ~ 0 + ", paste(lu_predictors, collapse = " + "),
        " + (1 | gr(phylo, cov = A))"))), data = d, data2 = list(A = A[d$phylo, d$phylo]),
        family = fp$family, prior = fp$prior))
}
# A2d: colour predictor, 4 studies
load("fits/A2d_model_setup.RData")
d0 <- dat_final[dat_final$SS %in% head(levels(droplevels(dat_final$SS)), 4), ]
for (cv in "z_dichroplate") {
  d <- droplevels(d0[!is.na(d0[[cv]]), ]); tips <- levels(d$phylo)
  ok(paste("A2d", cv), quick(formulas[[cv]], data = d, data2 = list(A = A[tips, tips]),
                             family = gaussian(), prior = priors_A2d))
}
# A3b/A3a: record level, 3000-record subset
for (m in c("A3b", "A3a")) {
  load(sprintf("fits/%s_model_setup.RData", m))
  d <- dat[sample(nrow(dat), 3000), ]
  for (resp in "dichroplate") {
    f <- if (m == "A3b") formulas[[resp]]$covariate else formulas[[resp]]
    fp <- family_prior(resp, "record")
    ok(paste(m, resp), quick(f, data = d, family = fp$family, prior = fp$prior))
  }
}
