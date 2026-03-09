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
# Filter: keep studies that include Primary vegetation
# ----------------------------------------------------------

primpresence <- by(fulldat, fulldat$Reference,
                   function(x) any(x$Predominant_simple == "Primary vegetation"))

fulldat_primary <- fulldat[fulldat$Reference %in%
                             names(primpresence)[which(primpresence)], ]

# Confirm each study has >1 land-use type
npred <- by(fulldat_primary, fulldat_primary$Reference,
            function(x) length(unique(x$Predominant_simple)))
table(npred)

# Random factor: site nested within study
fulldat_primary$site_study <- paste(fulldat_primary$Reference, "_",
                                    fulldat_primary$Site_number)

# Set reference level for land use
fulldat_primary$Predominant_simple <- relevel(
  as.factor(fulldat_primary$Predominant_simple),
  ref = "Primary vegetation"
)

# ----------------------------------------------------------
# Remove absences (zero-effort records) and save
# ----------------------------------------------------------

presdat <- fulldat[fulldat$Effort_corrected_measurement != 0, ]

write.csv(presdat, file = "data/present.csv", row.names = FALSE)

cat("Analytical dataset saved to data/present.csv\n")
cat("Rows:", nrow(presdat), "\n")
