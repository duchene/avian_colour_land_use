# ============================================================
# 00b_external_colour.R
# External colour data: download the published supplements and
# write the derived species-level tables the analyses use.
#
#   Cooney et al. 2022 Nat Ecol Evol 6:622 (doi:10.1038/s41559-022-01714-1)
#     Supplementary Data 1: male and female LociUVS and the plate-scored
#     sexual dichromatism (Dichro) for 4,527 passerines. Used for the
#     sample-bias test (C1) against all passerines.
#   Delhey et al. 2023 PNAS 120:e2217692120 (doi:10.1073/pnas.2217692120)
#     Dataset S7: proportion of plumage in each colour category, per sex,
#     for ~9,200 bird species named on the Jetz et al. 2012 tips.
#
# Raw downloads go to data/external/raw/ (gitignored). Only the derived
# tables in data/external/ are committed. Rerun this script to rebuild
# them from the original sources.
#
# Run from the repo root:  Rscript scripts/00b_external_colour.R
# ============================================================

library(readxl)

raw_dir <- "data/external/raw"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

options(timeout = 900)   # the Europe PMC supplement zip is slow to serve
fetch <- function(url, dest) {
  if (!file.exists(dest)) {
    cat("downloading", basename(dest), "\n")
    ok <- tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE) == 0,
                   error = function(e) FALSE)
    if (!ok) { unlink(dest); stop("download failed: ", url) }
  }
  dest
}

# ------------------------------------------------------------
# Cooney et al. 2022: species-level LociUVS and Dichro
# ------------------------------------------------------------
cooney_xlsx <- fetch(paste0("https://static-content.springer.com/esm/",
                            "art%3A10.1038%2Fs41559-022-01714-1/MediaObjects/",
                            "41559_2022_1714_MOESM3_ESM.xlsx"),
                     file.path(raw_dir, "cooney2022_SD1.xlsx"))
co <- as.data.frame(read_excel(cooney_xlsx, sheet = 1))
stopifnot(all(c("Binomial", "Sex", "LociUVS", "Dichro") %in% names(co)))

m <- co[co$Sex == "M", c("Binomial", "LociUVS", "Dichro")]
f <- co[co$Sex == "F", c("Binomial", "LociUVS")]
names(m) <- c("binomial", "LociUVS_male", "Dichro")
names(f) <- c("binomial", "LociUVS_female")
cooney <- merge(m, f, by = "binomial", all = TRUE)
cooney <- cooney[, c("binomial", "LociUVS_male", "LociUVS_female", "Dichro")]
cat("Cooney 2022:", nrow(cooney), "passerine species\n")
write.csv(cooney, "data/external/cooney2022_passerines.csv", row.names = FALSE)

# ------------------------------------------------------------
# Delhey et al. 2023: colour-category proportions, grouped
# ------------------------------------------------------------
# Supplement served as a zip by Europe PMC (PMC10450850).
delhey_zip <- fetch("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10450850/supplementaryFiles",
                    file.path(raw_dir, "delhey2023_supplement.zip"))
xl <- "pnas.2217692120.sd01.xlsx"
if (!file.exists(file.path(raw_dir, xl))) unzip(delhey_zip, files = xl, exdir = raw_dir)
del <- as.data.frame(read_excel(file.path(raw_dir, xl), sheet = "DatasetS7"))

# The 17 base categories partition the plumage (they sum to 1 per row);
# the .ex/.res/.alt/.sat/.unsat columns are alternative splits and are
# ignored. Column names carry the RGB centroid as a suffix.
base <- c("black", "white", "grey.dark", "grey.light", "brown.dark", "brown.light",
          "rufous.dark", "rufous.light", "yellow", "red.dark", "red.light",
          "green.dark", "green.light", "blue.dark", "blue.light",
          "purple.dark", "purple.light")
col_of <- function(b) {
  hit <- grep(paste0("^", gsub(".", "\\.", b, fixed = TRUE), "_[0-9]"), names(del), value = TRUE)
  stopifnot(length(hit) == 1)
  hit
}
P <- as.matrix(del[, vapply(base, col_of, "")]); colnames(P) <- base
stopifnot(abs(rowSums(P) - 1) < 1e-6)

groups <- list(black  = "black", white = "white",
               grey   = c("grey.dark", "grey.light"),
               brown  = c("brown.dark", "brown.light"),
               rufous = c("rufous.dark", "rufous.light"),
               yellow = "yellow",
               red    = c("red.dark", "red.light"),
               green  = c("green.dark", "green.light"),
               blue   = c("blue.dark", "blue.light"),
               purple = c("purple.dark", "purple.light"))
G <- sapply(groups, function(g) rowSums(P[, g, drop = FALSE]))

out <- data.frame(phylo = del$phylo, sex = del$sex, order = del$order_,
                  family = del$family, round(G, 6))
out <- out[out$order == "passeriformes", ]
cat("Delhey 2023:", length(unique(out$phylo)), "passerine species\n")
write.csv(out, "data/external/delhey2023_passerine_colour_groups.csv", row.names = FALSE)
cat("Wrote data/external/cooney2022_passerines.csv and",
    "data/external/delhey2023_passerine_colour_groups.csv\n")
