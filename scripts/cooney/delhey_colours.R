# ============================================================
# delhey_colours.R
# Male plumage colour proportions from Delhey et al. 2023 (PNAS,
# Dataset S7), one variable per colour.
#
# Decision of 2026-11-13 (co-authors): each colour is analysed as a
# separate variable, as the proportion of the plumage in that colour.
# No ratios, log-ratios or composite scores are formed. The ten colours
# are Delhey's 17 base categories with their dark and light shades
# combined (scripts/00b_external_colour.R); they partition the plumage,
# so the ten proportions of a species sum to one.
#
# Delhey names follow the BirdTree (Jetz) tips, the same labels as the
# `phylo` column of every model table, so species are matched on those.
# ============================================================

DELHEY_COLOURS <- c("black", "white", "grey", "brown", "rufous",
                    "yellow", "red", "green", "blue", "purple")
COLOUR_LABEL <- c(black = "Black", white = "White", grey = "Grey", brown = "Brown",
                  rufous = "Rufous", yellow = "Yellow", red = "Red", green = "Green",
                  blue = "Blue", purple = "Purple")

delhey_male <- function(path = "data/external/delhey2023_passerine_colour_groups.csv") {
  dl <- read.csv(path, stringsAsFactors = FALSE)
  dl[dl$sex == "male", c("phylo", DELHEY_COLOURS)]
}

# Adds one column per colour to `df`, matched on the tip label in `key`.
add_delhey_colours <- function(df, key) {
  dl <- delhey_male()
  m  <- match(as.character(key), dl$phylo)
  for (cc in DELHEY_COLOURS) df[[cc]] <- dl[[cc]][m]
  df
}
