#!/usr/bin/env Rscript

cat("=== Starting calc_alpha_diversity.R script ===\n")

library(argparse)
library(tidyverse)
library(mia)
library(TreeSummarizedExperiment)

# Source shared utilities
source("r_utils.R")

# Argument parsing
parser <- ArgumentParser(description = "Calculate alpha diversity indices for a TSE object")
parser$add_argument("--tse", type = "character", required = TRUE, help = "Path to input TSE RDS file")
parser$add_argument("--tse_type", type = "character", required = FALSE, default = NULL, help = "Type of TSE (for output naming, auto-detected if not provided)")
parser$add_argument("--output_dir", type = "character", default = ".", help = "Output directory for modified TSE files")
parser$add_argument("--seed", type = "integer", default = 14100715, help = "Random seed")
parser$add_argument("--indices", type = "character", default = "shannon,simpson_lambda,faith_diversity", help = "Comma-separated list of alpha diversity indices to calculate")
parser$add_argument("--names", type = "character", default = "Shannon,Simpson,Faith_PD", help = "Comma-separated list of names for the calculated indices")
parser$add_argument("--niter", type = "integer", default = 100, help = "Number of rarefaction iterations")
args <- parser$parse_args()

log_args(args)
set.seed(args$seed)

# Load and setup
tse <- load_tse(args$tse)
args$tse_type <- detect_tse_type(args$tse, args$tse_type)

# Parse indices and names
alpha_indices <- trimws(strsplit(args$indices, ",")[[1]])
alpha_names <- trimws(strsplit(args$names, ",")[[1]])

if (length(alpha_indices) != length(alpha_names)) {
  stop("Number of indices (", length(alpha_indices), ") must match number of names (", length(alpha_names), ")")
}

cat("Indices:", paste(alpha_indices, collapse = ", "), "\n")

# Check if counts assay exists
if (!"counts" %in% assayNames(tse)) {
  stop("Counts assay not found in TSE object. Required for alpha diversity calculation with rarefaction.")
}

# Calculate minimum library size for rarefaction
min_lib_size <- min(colSums(assay(tse, "counts")))
cat("Rarefaction depth:", min_lib_size, "\n")

cat("=== Calculating alpha diversity indices ===\n")
tse <- mia::addAlpha(
  tse,
  index = alpha_indices,
  name = alpha_names,
  assay.type = "counts",
  sample = min_lib_size,
  niter = args$niter
)

# Verify indices were added
for (name in alpha_names) {
  if (!name %in% colnames(colData(tse))) {
    warning("Index ", name, " not found in colData after calculation")
  }
}

# Save the modified TSE object
output_file <- file.path(args$output_dir, paste0("tse_", args$tse_type, "_with_alpha.rds"))
saveRDS(tse, file = output_file)
cat("Saved modified TSE to:", output_file, "\n")

# Save alpha diversity results separately
alpha_results <- list(
  alpha_diversity = colData(tse)[alpha_names],
  metadata = build_metadata(args, tse, params = list(
    indices_calculated = alpha_names,
    niter = args$niter,
    rarefaction_depth = min_lib_size
  ))
)

results_file <- file.path(args$output_dir, paste0("alpha_diversity_", args$tse_type, ".rds"))
saveRDS(alpha_results, file = results_file)
cat("Saved alpha diversity results to:", results_file, "\n")

cat("=== SUCCESS: Alpha diversity calculation completed for", args$tse_type, "TSE! ===\n")
