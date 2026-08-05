# ============================================================
# 00_data_preparation.R
# Prepare analytical dataset from raw PREDICTS + trait data
# ============================================================
#
# Input:  data/aves_predicts_passerines_dale_cooney_avonet.csv.gz
# Output: data/present.csv  (analytical dataset)
# ============================================================

fulldat <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"))

# ----------------------------------------------------------
# Compute colour response variables (Cooney UVS metrics)
# ----------------------------------------------------------

# Mean plumage colour (average of male + female)
fulldat$meancolcooney <- rowMeans(fulldat[, c("LociUVS_male_cooney",
                                               "LociUVS_female_cooney")],
                                   na.rm = TRUE)

# Sexual dichromatism — ratio (male / female)
fulldat$dichrocooney <- fulldat$LociUVS_male_cooney /
                        fulldat$LociUVS_female_cooney

# Sexual dichromatism — absolute difference (male - female)
fulldat$dichrodiff <- fulldat$LociUVS_male_cooney -
                      fulldat$LociUVS_female_cooney

# Male plumage colour only
fulldat$malecolcooney <- fulldat$LociUVS_male_cooney

# ----------------------------------------------------------
# Compute colour response variables (Dale metrics)
# ----------------------------------------------------------
# Dale male/female plumage scores (Dale et al. 2015), the parallel of
# the Cooney UVS metrics above. Used by the B-family analyses (B1a, B1b,
# B2c, B2d, B3), which mirror A1a/A1b/A2c/A2d/A3 with Dale colour.

# Mean plumage colour (average of male + female)
fulldat$meancoldale <- rowMeans(fulldat[, c("Male_plumage_score_dale",
                                            "Female_plumage_score_dale")],
                                na.rm = TRUE)

# Sexual dichromatism, ratio (male / female)
fulldat$dichrodale <- fulldat$Male_plumage_score_dale /
                      fulldat$Female_plumage_score_dale

# Sexual dichromatism, absolute difference (male - female)
fulldat$dichrodiffdale <- fulldat$Male_plumage_score_dale -
                          fulldat$Female_plumage_score_dale

# Male plumage colour only
fulldat$malecoldale <- fulldat$Male_plumage_score_dale

# ----------------------------------------------------------
# Simplified biome classification (11 -> 4 levels)
# ----------------------------------------------------------
# Reduces parameter count and enables Biome x Land-use interaction
fulldat$Biome4 <- ifelse(
  fulldat$Biome %in% c("Tropical & Subtropical Moist Broadleaf Forests",
                        "Tropical & Subtropical Dry Broadleaf Forests",
                        "Tropical & Subtropical Coniferous Forests"),
  "Tropical Forest",
  ifelse(
    fulldat$Biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands",
                          "Mangroves"),
    "Tropical Open",
    ifelse(
      fulldat$Biome %in% c("Temperate Broadleaf & Mixed Forests",
                            "Temperate Conifer Forests",
                            "Mediterranean Forests, Woodlands & Scrub"),
      "Temperate Forest",
      "Temperate Open"  # Temperate Grasslands, Montane Grasslands, Tundra
    )
  )
)

# ----------------------------------------------------------
# Step 1: remove absences (zero-abundance records)
# ----------------------------------------------------------
# PREDICTS stores a near-complete site x species matrix per study, so
# ~84% of raw rows are structural zeros for species not detected at a
# site. Only presences carry colour/abundance information here.

presdat <- fulldat[fulldat$Effort_corrected_measurement != 0, ]

cat("Step 1  drop zero-abundance rows:", nrow(fulldat), "->", nrow(presdat), "\n")

# ----------------------------------------------------------
# Step 2: keep only studies that contrast Primary vegetation
#         with at least one other land use
# ----------------------------------------------------------
# The analyses model land use against a Primary vegetation reference, so
# a study only contributes a within-study contrast if it holds Primary
# vegetation AND at least one other land-use type.
#
# Applied AFTER step 1, on the records that actually enter the models.
# Testing on the raw file (zeros included) would only ask whether
# primary-vegetation sites were sampled, which a study can satisfy while
# contributing no non-zero primary-vegetation record.
#
# FILTER_UNIT: "Reference" (publication) is the original choice and the
# more conservative, since a publication can hold several PREDICTS
# studies. "SS" (the study carrying the (1|SS) intercept) differs by
# only 32 records here.

FILTER_UNIT <- "Reference"   # or "SS"

unit <- presdat[[FILTER_UNIT]]

has_primary <- tapply(presdat$Predominant_simple, unit,
                      function(x) any(x == "Primary vegetation"))
n_landuse   <- tapply(presdat$Predominant_simple, unit,
                      function(x) length(unique(x)))

keep_units <- names(has_primary)[has_primary & n_landuse >= 2]

cat("Step 2  filter by", FILTER_UNIT, ":", length(has_primary), "units ->",
    length(keep_units), "kept\n")
cat("        dropped, no Primary vegetation:", sum(!has_primary), "\n")
cat("        dropped, Primary vegetation but only one land use:",
    sum(has_primary & n_landuse < 2), "\n")
cat("        dropped units:",
    paste(setdiff(names(has_primary), keep_units), collapse = ", "), "\n")

presdat <- presdat[unit %in% keep_units, ]

# Every retained unit must now satisfy the stated criterion.
stopifnot(all(tapply(presdat$Predominant_simple, presdat[[FILTER_UNIT]],
                     function(x) any(x == "Primary vegetation") &&
                                 length(unique(x)) >= 2)))

# ----------------------------------------------------------
# Step 3: derived factors and save
# ----------------------------------------------------------

# Set reference level for land use
presdat$Predominant_simple <- relevel(
  as.factor(presdat$Predominant_simple),
  ref = "Primary vegetation"
)

write.csv(presdat, file = "data/present.csv", row.names = FALSE)

cat("\nAnalytical dataset saved to data/present.csv\n")
cat("Rows:", nrow(presdat),
    "| species:", length(unique(presdat$Best_guess_binomial)),
    "| SS:", length(unique(presdat$SS)),
    "| SSB:", length(unique(presdat$SSB)),
    "| SSBS:", length(unique(presdat$SSBS)), "\n\n")
cat("Records per land use:\n")
print(table(presdat$Predominant_simple))
cat("\nBiome4 x land use, records (check for sparse cells before fitting):\n")
print(table(presdat$Biome4, presdat$Predominant_simple))
cat("\nBiome4 x land use, species:\n")
print(tapply(presdat$Best_guess_binomial,
             list(presdat$Biome4, presdat$Predominant_simple),
             function(x) length(unique(x))))
