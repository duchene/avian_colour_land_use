# ============================================================
# 00_data_preparation.R
# Build the analytical dataset from raw PREDICTS + trait data
# ============================================================
#
# Input:  data/aves_predicts_passerines_dale_cooney_avonet.csv.gz
# Output: data/present.csv  (analytical dataset)
#
# Run first, from the repository root. Every shared decision
# (biome map, biome drop, land-use reference) lives in
# scripts/00_config.R, not here.
# ============================================================

source("scripts/00_config.R")

fulldat <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"))

# ------------------------------------------------------------
# Colour metrics
# ------------------------------------------------------------
# Cooney UVS (A family): LociUVS loci counts under the
# ultraviolet-sensitive visual model.
fulldat$meancolcooney <- rowMeans(fulldat[, c("LociUVS_male_cooney",
                                              "LociUVS_female_cooney")],
                                  na.rm = TRUE)
fulldat$malecolcooney <- fulldat$LociUVS_male_cooney
fulldat$dichrocooney  <- fulldat$LociUVS_male_cooney / fulldat$LociUVS_female_cooney
fulldat$dichrodiff    <- fulldat$LociUVS_male_cooney - fulldat$LociUVS_female_cooney

# Dale (B family): Dale et al. 2015 plumage scores, the parallel set.
fulldat$meancoldale    <- rowMeans(fulldat[, c("Male_plumage_score_dale",
                                               "Female_plumage_score_dale")],
                                   na.rm = TRUE)
fulldat$malecoldale    <- fulldat$Male_plumage_score_dale
fulldat$dichrodale     <- fulldat$Male_plumage_score_dale /
                          fulldat$Female_plumage_score_dale
fulldat$dichrodiffdale <- fulldat$Male_plumage_score_dale -
                          fulldat$Female_plumage_score_dale

# ------------------------------------------------------------
# Step 1: drop absences (zero-abundance records)
# ------------------------------------------------------------
# PREDICTS stores a near-complete site x species matrix per study,
# so ~84% of raw rows are structural zeros for species not detected
# at a site. Only presences carry colour/abundance information here.
presdat <- fulldat[fulldat$Effort_corrected_measurement != 0, ]
cat("Step 1  drop zero-abundance rows:", nrow(fulldat), "->", nrow(presdat), "\n")

# ------------------------------------------------------------
# Step 2: apply the biome policy
# ------------------------------------------------------------
# Assigns Biome4 and drops Temperate Open. Runs BEFORE the
# primary-vegetation filter, so a study is judged on the records
# that actually enter the models.
presdat <- apply_biome_policy(presdat)
cat("Step 2  biome policy applied ->", nrow(presdat), "rows\n")

# ------------------------------------------------------------
# Step 3: keep only units contrasting Primary vegetation
#         with at least one other land use
# ------------------------------------------------------------
# The analyses model land use against a Primary vegetation
# reference, so a unit only contributes a within-unit contrast if
# it holds Primary vegetation AND at least one other land use.
#
# FILTER_UNIT: "Reference" (the publication) is the original and
# more conservative choice, since one publication can hold several
# PREDICTS studies. "SS" (the study carrying the (1|SS) intercept)
# differs by only a few dozen records.
FILTER_UNIT <- "Reference"   # or "SS"

unit        <- presdat[[FILTER_UNIT]]
has_primary <- tapply(presdat$Predominant_simple, unit,
                      function(x) any(x == LU_REF))
n_landuse   <- tapply(presdat$Predominant_simple, unit,
                      function(x) length(unique(x)))
keep_units  <- names(has_primary)[has_primary & n_landuse >= 2]

cat("Step 3  filter by", FILTER_UNIT, ":", length(has_primary), "units ->",
    length(keep_units), "kept\n")
cat("        dropped, no Primary vegetation:", sum(!has_primary), "\n")
cat("        dropped, Primary vegetation but only one land use:",
    sum(has_primary & n_landuse < 2), "\n")

presdat <- presdat[unit %in% keep_units, ]

# The saved file must satisfy the stated criterion. This assertion
# is the reason the 2026-08 filter bug cannot recur silently.
stopifnot(all(tapply(presdat$Predominant_simple, presdat[[FILTER_UNIT]],
                     function(x) any(x == LU_REF) && length(unique(x)) >= 2)))

# ------------------------------------------------------------
# Step 4: factor levels and save
# ------------------------------------------------------------
presdat <- set_landuse(presdat)

write.csv(presdat, file = "data/present.csv", row.names = FALSE)

cat("\nAnalytical dataset saved to data/present.csv\n")
cat("Rows:", nrow(presdat),
    "| species:", length(unique(presdat$Best_guess_binomial)),
    "| SS:", length(unique(presdat$SS)),
    "| SSB:", length(unique(presdat$SSB)),
    "| SSBS:", length(unique(presdat$SSBS)), "\n\n")
cat("Records per land use:\n")
print(table(presdat$Predominant_simple))
check_cells(presdat)
