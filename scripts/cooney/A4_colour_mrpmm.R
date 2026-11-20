# ============================================================
# A4_colour_mrpmm.R
# Analysis 4: multi-response phylogenetic mixed model (MR-PMM) of male
# plumage colours and land-use association, at the species level.
#
# Following Halliwell, Holland & Yates (2025, Biol Rev 100:1294-1316,
# doi:10.1111/brv.70001), every trait is a RESPONSE, including the
# land-use association scores. The covariance between traits is
# partitioned into a phylogenetic and a residual (non-phylogenetic)
# matrix, and the relationships of interest are derived from those
# matrices rather than fixed in advance as predictor -> response.
#
# Responses (12), for the A1b species that have Delhey et al. 2023 data:
#   eight colours: the proportion of male plumage in black, white, brown,
#     rufous, yellow, red, green and blue. Each colour is its own
#     variable; no ratios are formed (co-author decision, 2026-11-13).
#     Purple is left out: 96% of species carry less than 0.5% purple, and
#     it failed the posterior predictive check on every scale tried
#     (pilot of 2026-11-13; decision 2026-11-15).
#   four land-use association scores: the proportion of a species'
#     detection records in secondary vegetation, plantation forest,
#     cropland and pasture (as in A1b).
# Grey and primary vegetation are also left out. The ten colours and the five
# land-use scores each sum to one, so including all of either set makes
# the covariance matrices singular. Grey is the commonest colour and is
# not part of any hypothesis. Primary vegetation becomes the implicit
# reference: a slope on a land-use score is the change expected when a
# species' records shift into that land use away from primary
# vegetation, the other three scores held constant.
#
# Every response is standardised (mean 0, SD 1) so that the identity
# prior scale suits all of them; correlations are unaffected, and
# slopes are also reported per unit of the fitted (logit) scale.
#
# Model (MCMCglmm, Hadfield 2010), as in the Halliwell et al. tutorial:
#   cbind(12 responses) ~ trait - 1
#   random = ~ us(trait):phylo  (BBtree2, inverse relationship matrix)
#   rcov   = ~ us(trait):units
#   priors: inverse-Wishart, V = I, nu = k + 1, for both matrices.
#
# Derived per posterior draw, at the phylogenetic, residual and total
# (phylogenetic + residual) levels:
#   correlations among all responses (are the colours orthogonal?)
#   conditional slopes of each colour on the four land-use scores,
#     B = S[colour, LU] %*% solve(S[LU, LU])
#   phylogenetic signal of each response, S_phy / (S_phy + S_res)
# The total level is the species-level association, the analogue of A1b.
#
# Posterior predictive check: the colours are skewed, bounded
# proportions with many species lacking a colour, fitted as Gaussian.
# For each response, data are simulated from the fitted model (with
# phylogenetically correlated effects) and compared with the observed
# distribution: skewness, and the share of predictions below zero, an
# impossible proportion.
#
# Response scale (MRPMM_SCALE): "logit" (default) fits each variable's
# own logit, after the
# Smithson & Verkuilen (2006) squeeze (p (n - 1) + 0.5) / n that keeps 0
# and 1 finite. Each variable is transformed on its own; nothing is
# combined across colours. "raw" fits the proportions as they are; in the
# pilot it predicted 13-48% negative proportions and reproduced none of
# the skew, whereas the logit scale reproduced both for every response
# kept here (decision 2026-11-15).
#
# Run from the repo root after A1b_01_setup.R (hpc/15_mrpmm.sh).
# Environment overrides: MRPMM_NITT, MRPMM_BURN, MRPMM_THIN, MRPMM_SCALE,
# MRPMM_TAG (suffix for the fit, result and figure files).
# ============================================================

source("scripts/00_config.R")
source("scripts/cooney/delhey_colours.R")
suppressPackageStartupMessages({ library(MCMCglmm); library(ape); library(posterior)
                                  library(ggplot2) })

a1b <- new.env(); load("fits/A1b_model_setup.RData", envir = a1b)
spdat <- add_delhey_colours(a1b$spdat, a1b$spdat$phylo)
tree  <- a1b$tree

COLS <- setdiff(DELHEY_COLOURS, c("grey", "purple"))
LU   <- c(Secondary = "prop_Secondary", Plantation = "prop_Plantation_forest",
          Cropland = "prop_Cropland", Pasture = "prop_Pasture")
raw  <- spdat[, c("phylo", COLS, unname(LU))]
names(raw) <- c("phylo", COLS, names(LU))
raw  <- raw[complete.cases(raw) & raw$phylo %in% tree$tip.label, ]
TR   <- c(COLS, names(LU)); k <- length(TR)
SCALE <- Sys.getenv("MRPMM_SCALE", "logit"); stopifnot(SCALE %in% c("raw", "logit"))
nsq  <- nrow(raw)
fwd  <- if (SCALE == "logit") function(p) qlogis((p * (nsq - 1) + 0.5) / nsq) else identity
bwd  <- if (SCALE == "logit") function(y) (plogis(y) * nsq - 0.5) / (nsq - 1) else identity
mod  <- as.data.frame(lapply(raw[TR], fwd))
ctr  <- colMeans(mod); scl <- apply(mod, 2, sd)
d    <- data.frame(phylo = raw$phylo, scale(mod))
cat("A4 MR-PMM:", nrow(d), "species,", k, "responses, scale:", SCALE, "\n")

tr <- drop.tip(tree, setdiff(tree$tip.label, d$phylo))
tr$node.label <- NULL
tr <- makeNodeLabel(tr)
Ainv <- inverseA(tr, nodes = "ALL", scale = TRUE)$Ainv

fixed <- as.formula(paste0("cbind(", paste(TR, collapse = ", "), ") ~ trait - 1"))
prior <- list(R = list(V = diag(k), nu = k + 1),
              G = list(G1 = list(V = diag(k), nu = k + 1)))

env_int <- function(v, def) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else def }
CHAINS <- 4
NITT <- env_int("MRPMM_NITT", 105000); BURN <- env_int("MRPMM_BURN", 5000)
THIN <- env_int("MRPMM_THIN", 50)
TAG  <- Sys.getenv("MRPMM_TAG", "")
OUT  <- paste0("fits/A4_mrpmm", TAG, ".RData")
RES  <- function(x) sprintf("results/cooney/A4%s_%s", TAG, x)

t0 <- Sys.time()
fits <- parallel::mclapply(seq_len(CHAINS), function(i) {
  set.seed(2026 + i)
  MCMCglmm(fixed, random = ~ us(trait):phylo, rcov = ~ us(trait):units,
           family = rep("gaussian", k), ginverse = list(phylo = Ainv),
           prior = prior, data = d, nitt = NITT, burnin = BURN, thin = THIN,
           verbose = FALSE, pr = FALSE)
}, mc.cores = CHAINS)
stopifnot(all(vapply(fits, inherits, TRUE, "MCMCglmm")))
cat("Sampling took", format(Sys.time() - t0), "for nitt =", NITT, "\n")
dir.create("fits", showWarnings = FALSE)
save(fits, d, raw, TR, COLS, LU, ctr, scl, SCALE, file = OUT)

dir.create("results/cooney", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/cooney", showWarnings = FALSE, recursive = TRUE)
q <- function(x) c(median = median(x), lower = unname(quantile(x, 0.025)),
                   upper = unname(quantile(x, 0.975)))

# ---------- convergence -------------------------------------------------
arr <- function(get) {
  m <- lapply(fits, function(f) as.matrix(get(f)))
  a <- array(unlist(m), dim = c(nrow(m[[1]]), ncol(m[[1]]), CHAINS))
  a <- aperm(a, c(1, 3, 2)); dimnames(a) <- list(NULL, NULL, colnames(m[[1]])); a
}
sm <- summarise_draws(as_draws_array(arr(function(f) f$Sol)), "rhat", "ess_bulk", "ess_tail")
sv <- summarise_draws(as_draws_array(arr(function(f) f$VCV)), "rhat", "ess_bulk", "ess_tail")
sv <- sv[!is.na(sv$rhat), ]
conv <- data.frame(n_species = nrow(d), n_responses = k, chains = CHAINS,
                   draws_per_chain = nrow(fits[[1]]$VCV),
                   max_rhat = max(c(sm$rhat, sv$rhat)),
                   min_ess_bulk = min(c(sm$ess_bulk, sv$ess_bulk)),
                   min_ess_tail = min(c(sm$ess_tail, sv$ess_tail)))
print(conv)
write.csv(conv, RES("convergence.csv"), row.names = FALSE)

# ---------- covariance matrices per draw --------------------------------
VCV <- do.call(rbind, lapply(fits, function(f) as.matrix(f$VCV)))
SOL <- do.call(rbind, lapply(fits, function(f) as.matrix(f$Sol)))
S <- function(level) {
  a <- array(NA_real_, c(nrow(VCV), k, k), dimnames = list(NULL, TR, TR))
  for (i in TR) for (j in TR) a[, i, j] <- VCV[, sprintf("trait%s:trait%s.%s", i, j, level)]
  a
}
LEV <- list(phylogenetic = S("phylo"), residual = S("units"))
LEV$total <- LEV$phylogenetic + LEV$residual

# ---------- correlations among all responses ----------------------------
cors <- list()
for (lv in names(LEV)) {
  a <- LEV[[lv]]
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    r <- a[, i, j] / sqrt(a[, i, i] * a[, j, j])
    cors[[length(cors) + 1]] <- data.frame(level = lv, trait_1 = TR[i], trait_2 = TR[j],
      t(q(r)), excludes_zero = prod(quantile(r, c(0.025, 0.975))) > 0)
  }
}
cors <- do.call(rbind, cors)
rawcor <- cor(raw[TR], method = "spearman")
cors$raw_spearman <- mapply(function(a, b) rawcor[a, b], cors$trait_1, cors$trait_2)
cors$pair_type <- ifelse(cors$trait_1 %in% COLS & cors$trait_2 %in% COLS, "colour-colour",
                  ifelse(cors$trait_1 %in% COLS | cors$trait_2 %in% COLS, "colour-land use",
                         "land use-land use"))
write.csv(cors, RES("correlations.csv"), row.names = FALSE)

# ---------- phylogenetic signal -----------------------------------------
h2 <- do.call(rbind, lapply(TR, function(t) {
  p <- LEV$phylogenetic[, t, t]; u <- LEV$residual[, t, t]
  data.frame(trait = t, t(q(p / (p + u))))
}))
write.csv(h2, RES("phylogenetic_signal.csv"), row.names = FALSE)

# ---------- conditional slopes of each colour on land use ---------------
L <- names(LU)
slopes <- list()
for (lv in names(LEV)) {
  a <- LEV[[lv]]
  B <- array(NA_real_, c(nrow(VCV), length(COLS), length(L)), dimnames = list(NULL, COLS, L))
  for (s in seq_len(nrow(VCV)))
    B[s, , ] <- a[s, COLS, L, drop = FALSE][1, , ] %*% solve(a[s, L, L])
  for (cc in COLS) for (l in L) {
    b <- B[, cc, l]
    slopes[[length(slopes) + 1]] <- data.frame(level = lv, colour = cc, land_use = l,
      setNames(as.data.frame(t(q(b))), c("sd_median", "sd_lower", "sd_upper")),
      # model scale (raw proportion or logit): change in the colour per unit
      # (0 -> 1) shift of a species' records into this land use away from
      # primary vegetation
      setNames(as.data.frame(t(q(b * scl[[cc]] / scl[[l]]))), c("unit_median", "unit_lower", "unit_upper")),
      prob_positive = mean(b > 0),
      excludes_zero = prod(quantile(b, c(0.025, 0.975))) > 0)
  }
}
slopes <- do.call(rbind, slopes)
write.csv(slopes, RES("colour_landuse_slopes.csv"), row.names = FALSE)

# ---------- posterior predictive check ----------------------------------
set.seed(1)
NREP <- 300
idx  <- sample(nrow(VCV), NREP)
Lc   <- t(chol(vcv(tr, corr = TRUE)[d$phylo, d$phylo]))
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
ppc <- list(); dens <- list()
for (t in TR) {
  obs <- raw[[t]]
  rep_stats <- t(vapply(idx, function(s) {
    y <- SOL[s, paste0("trait", t)] +
         sqrt(VCV[s, sprintf("trait%s:trait%s.phylo", t, t)]) * as.vector(Lc %*% rnorm(nrow(d))) +
         sqrt(VCV[s, sprintf("trait%s:trait%s.units", t, t)]) * rnorm(nrow(d))
    y <- bwd(y * scl[[t]] + ctr[[t]])   # back on the proportion scale
    c(skew = skew(y), low = mean(y < 0.005), below0 = mean(y < 0))
  }, numeric(3)))
  # "low": share of species with < 0.5% of plumage in this colour (or
  # association score), observed against replicated.
  ppc[[t]] <- data.frame(trait = t, scale = SCALE,
    obs_low_share = mean(obs < 0.005), rep_low_median = median(rep_stats[, "low"]),
    obs_skew = skew(obs), rep_skew_median = median(rep_stats[, "skew"]),
    ppp_skew = mean(rep_stats[, "skew"] >= skew(obs)),
    rep_below0_median = median(rep_stats[, "below0"]))
  # a few replicated data sets for the overlay figure
  for (s in idx[1:20]) {
    y <- SOL[s, paste0("trait", t)] +
         sqrt(VCV[s, sprintf("trait%s:trait%s.phylo", t, t)]) * as.vector(Lc %*% rnorm(nrow(d))) +
         sqrt(VCV[s, sprintf("trait%s:trait%s.units", t, t)]) * rnorm(nrow(d))
    dens[[length(dens) + 1]] <- data.frame(trait = t, draw = s, value = bwd(y * scl[[t]] + ctr[[t]]))
  }
}
ppc <- do.call(rbind, ppc)
write.csv(ppc, RES("ppc.csv"), row.names = FALSE)
print(ppc, digits = 3, row.names = FALSE)

lab <- c(COLOUR_LABEL, Secondary = "Secondary", Plantation = "Plantation",
         Cropland = "Cropland", Pasture = "Pasture")
dd <- do.call(rbind, dens); dd$trait <- factor(lab[dd$trait], lab[TR])
ob <- data.frame(trait = factor(rep(lab[TR], each = nrow(raw)), lab[TR]),
                 value = unlist(raw[TR], use.names = FALSE))
p <- ggplot() +
  geom_density(data = dd, aes(value, group = draw), colour = "#2a6f97", alpha = 0.3, linewidth = 0.2) +
  geom_density(data = ob, aes(value), colour = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = 3) +
  facet_wrap(~ trait, scales = "free", ncol = 5) +
  labs(x = "Proportion (raw scale)", y = "Density") +
  theme_classic(base_size = 9)
ggsave(sprintf("figures/cooney/A4%s_ppc.png", TAG), p, width = 9, height = 5.5, dpi = 300)
cat("\nA4 complete:", RES("*.csv"), "\n")
