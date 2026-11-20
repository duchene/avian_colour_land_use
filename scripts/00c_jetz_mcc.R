# ============================================================
# 00c_jetz_mcc.R
# Second phylogeny for the tree-sensitivity runs: the maximum clade
# credibility (MCC) tree of 1,000 posterior trees of Jetz et al. (2012,
# Nature 491:444-448), Stage 2 (all species), Hackett et al. backbone.
#
# The reported analyses use the single maximum-likelihood time tree of
# Claramunt et al. (2025), data/BBtree2.tre. Every phylogenetic analysis is
# repeated on this tree, so the study rests on two independent trees.
#
# Input (gitignored download, rebuild with):
#   mkdir -p data/external/raw/jetz2012 && cd data/external/raw/jetz2012
#   curl -LO https://data.vertlife.org/birdtree/Stage2/HackettStage2_0001_1000.zip
#   unzip -j HackettStage2_0001_1000.zip      # AllBirdsHackett1.tre
# The 1,000 trees are pruned to every passerine any analysis can use (the
# Cooney et al. 2022 and Delhey et al. 2023 passerines, the PREDICTS
# passerines and the tree synonyms in 00_config.R) before the MCC tree is
# chosen (phangorn::maxCladeCred, rooted). Branch lengths are those of the
# selected posterior tree, in millions of years.
#
# Output: data/external/jetz2012_hackett_mcc_passerines.tre
# Run on a compute node (about 30 GB of memory).
# ============================================================

source("scripts/00_config.R")
suppressPackageStartupMessages({ library(ape); library(phangorn) })

co  <- read.csv("data/external/cooney2022_passerines.csv", stringsAsFactors = FALSE)
dl  <- read.csv("data/external/delhey2023_passerine_colour_groups.csv", stringsAsFactors = FALSE)
raw <- read.csv(gzfile("data/aves_predicts_passerines_dale_cooney_avonet.csv.gz"),
                stringsAsFactors = FALSE)
want <- unique(c(co$binomial, dl$phylo, raw$jetz_sp, gsub(" ", "_", raw$Best_guess_binomial),
                 unname(TREE_SYNONYMS)))

f  <- "data/external/raw/jetz2012/AllBirdsHackett1.tre"
tr <- read.tree(f)
stopifnot(length(tr) == 1000)
keep <- intersect(want, tr[[1]]$tip.label)
cat("Jetz trees:", length(tr), "| tips per tree:", Ntip(tr[[1]]),
    "| passerine tips kept:", length(keep), "of", length(want), "wanted\n")
tr <- lapply(tr, keep.tip, tip = keep)
class(tr) <- "multiPhylo"
tr <- .compressTipLabel(tr)
mcc <- maxCladeCred(tr, rooted = TRUE)
stopifnot(is.ultrametric(mcc, tol = 1e-6), Ntip(mcc) == length(keep))

out <- "data/external/jetz2012_hackett_mcc_passerines.tre"
write.tree(mcc, out)
bb <- read.tree("data/BBtree2.tre")
shared <- intersect(mcc$tip.label, bb$tip.label)
cat("Saved", out, "|", Ntip(mcc), "tips | root age", round(max(node.depth.edgelength(mcc)), 1),
    "Myr | tips shared with BBtree2:", length(shared), "\n")
