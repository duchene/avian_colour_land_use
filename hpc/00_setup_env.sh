#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=colour_env
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#
# One-off: build the R environment for this repo on GenomeDK.
# R 4.4.1 comes from the shared rocker container; packages missing from it
# (cmdstanr + CmdStan, tidyverse, phylolm) go into a separate user library so the
# avian_flight_genomics library is not touched.
# Usage, from the repo root:  sbatch hpc/00_setup_env.sh
set -euo pipefail
source hpc/env.sh
mkdir -p "$RLIB_COLOUR" "$CMDSTAN_DIR"
$RSCRIPT -e '
lib <- Sys.getenv("RLIB_COLOUR")
repos <- c(stan = "https://stan-dev.r-universe.dev", CRAN = getOption("repos")[["CRAN"]])
need <- c("cmdstanr", "tidyverse", "phylolm", "phytools")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, lib = lib, repos = repos, Ncpus = 8)
library(cmdstanr)
if (is.null(tryCatch(cmdstan_path(), error = function(e) NULL)))
  install_cmdstan(dir = Sys.getenv("CMDSTAN_DIR"), cores = 8, quiet = TRUE)
cat("cmdstan:", cmdstan_path(), cmdstan_version(), "\n")
for (p in c("brms", need)) cat(p, as.character(packageVersion(p)), "\n")'
