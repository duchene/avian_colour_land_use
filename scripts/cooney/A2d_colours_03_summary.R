# ============================================================
# A2d_colours_03_summary.R
# Convergence and colour slopes for the per-colour A2d fits
# (A2d_colours_02_fit.R). For each colour, the slope of abundance change
# on the standardised colour proportion (per 1 SD), overall and within
# each land use and each biome, with random effects excluded; slopes are
# computed exactly as for the reported A2d (posthoc_helpers.R).
#
# Outputs: results/cooney/A2d_colours_convergence.csv
#          results/cooney/A2d_colours_slopes.csv
# ============================================================

source("scripts/cooney/posthoc_helpers.R")   # attaches rstan and brms
source("scripts/cooney/delhey_colours.R")
suppressPackageStartupMessages(library(posterior))

conv <- list(); rows <- list()
for (cc in DELHEY_COLOURS) {
  cv <- paste0("z_col_", cc)
  f  <- sprintf("fits/A2d_fit_%s.RData", cv)
  if (!file.exists(f)) { message("missing ", f, ", skipped"); next }
  load(f)
  sm <- summarise_draws(as_draws_array(fit, variable = "^(b_|sd_|sigma)", regex = TRUE),
                        "rhat", "ess_bulk", "ess_tail")
  np <- nuts_params(fit)
  conv[[cc]] <- data.frame(colour = COLOUR_LABEL[[cc]], n_records = nobs(fit),
    n_species = length(unique(fit$data$phylo)),
    max_rhat = max(sm$rhat), min_ess_bulk = min(sm$ess_bulk), min_ess_tail = min(sm$ess_tail),
    n_divergent = sum(np$Value[np$Parameter == "divergent__"]))

  g <- make_grid(Predominant_simple = LU4, Biome4 = BIOME3)
  g[[cv]] <- 0
  g$overall <- "All land uses and biomes"
  sdl <- link_sd(fit)
  for (wi in c("overall", "Predominant_simple", "Biome4"))
    rows[[paste(cc, wi)]] <- columnwise(slope_draws(fit, g, wi, cv), sdl, FALSE,
      colour = COLOUR_LABEL[[cc]],
      slope_within = c(overall = "Overall", Predominant_simple = "Land use", Biome4 = "Biome")[[wi]])
  rm(fit); gc(verbose = FALSE)
}
conv <- do.call(rbind, conv)
print(conv, row.names = FALSE, digits = 3)
write.csv(conv, "results/cooney/A2d_colours_convergence.csv", row.names = FALSE)

slopes <- bind_rows(rows) %>%
  rename(level = contrast) %>%
  mutate(level = relabel_levels(level)) %>%
  select(colour, slope_within, level, estimate, lo, hi, pd, credible, effect_size)
write.csv(slopes, "results/cooney/A2d_colours_slopes.csv", row.names = FALSE)
print(as.data.frame(slopes %>% filter(credible)), digits = 3)
cat("Saved results/cooney/A2d_colours_{convergence,slopes}.csv\n")
