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

# Sexual dichromatism (male / female ratio)
fulldat$dichrocooney <- fulldat$LociUVS_male_cooney /
                        fulldat$LociUVS_female_cooney

# Male plumage colour only
fulldat$malecolcooney <- fulldat$LociUVS_male_cooney

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
