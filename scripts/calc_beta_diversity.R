#!/usr/bin/env Rscript

cat("=== Starting calc_beta_diversity.R script ===\n")

library(argparse)
library(tidyverse)
library(mia)
library(TreeSummarizedExperiment)

source("r_utils.R")

# Argument parsing
parser <- ArgumentParser(description = "Calculate beta diversity distance matrices and MDS ordination for a TSE object")
parser$add_argument("--tse", type = "character", required = TRUE, help = "Path to input TSE RDS file")
parser$add_argument("--tse_type", type = "character", required = FALSE, default = NULL, help = "Type of TSE (for output naming, auto-detected if not provided)")
parser$add_argument("--output_dir", type = "character", default = ".", help = "Output directory for modified TSE files")
parser$add_argument("--seed", type = "integer", default = 14100715, help = "Random seed")
parser$add_argument("--distance_methods", type = "character", default = "jaccard,bray,unifrac_unweighted,unifrac_weighted", help = "Comma-separated list of distance methods to calculate")
parser$add_argument("--ncomponents", type = "integer", default = 3, help = "Number of MDS components to calculate")
parser$add_argument("--niter", type = "integer", default = 100, help = "Number of rarefaction iterations for distance calculations")
args <- parser$parse_args()

log_args(args)
set.seed(args$seed)

# Load and setup
tse <- load_tse(args$tse)
args$tse_type <- detect_tse_type(args$tse, args$tse_type)

# Parse distance methods
distance_methods <- trimws(strsplit(args$distance_methods, ",")[[1]])
cat("Methods:", paste(distance_methods, collapse = ", "), "\n")

# Check required assays
required_assays <- c("relabundance", "counts")
missing_assays <- setdiff(required_assays, assayNames(tse))
if (length(missing_assays) > 0) {
  stop("Missing required assays: ", paste(missing_assays, collapse = ", "))
}

# Check if tree exists for UniFrac calculations
has_tree <- !is.null(rowTree(tse))
if (!has_tree && any(grepl("unifrac", distance_methods))) {
  warning("Row tree not found - UniFrac calculations will be skipped")
  distance_methods <- distance_methods[!grepl("unifrac", distance_methods)]
}

niter <- args$niter

cat("=== Calculating beta diversity distance matrices ===\n")

for (method in distance_methods) {
  cat("Calculating", method, "...\n")

  if (method == "jaccard") {
    tse <- mia::transformAssay(tse, method = "pa")
    tse <- mia::addDissimilarity(tse, assay.type = "pa", method = "jaccard", niter = niter)
  } else if (method == "bray") {
    tse <- mia::addDissimilarity(tse, method = "bray", name = "bray", niter = niter)
  } else if (method == "unifrac_unweighted") {
    if (has_tree) {
      tse <- mia::addDissimilarity(tse, method = "unifrac", tree = rowTree(tse), weighted = FALSE, niter = niter)
    } else {
      cat("Skipping unweighted UniFrac - no tree available\n")
    }
  } else if (method == "unifrac_weighted") {
    if (has_tree) {
      tse <- mia::addDissimilarity(tse, assay.type = "counts", method = "unifrac",
        tree = rowTree(tse), weighted = TRUE, name = "unifrac_weighted", niter = niter)
    } else {
      cat("Skipping weighted UniFrac - no tree available\n")
    }
  } else {
    stop("Unknown distance method: ", method)
  }
}

cat("=== Calculating MDS ordinations ===\n")

for (method in distance_methods) {
  cat("Calculating MDS for", method, "...\n")

  # Map method name to metadata key (unifrac_unweighted is stored as "unifrac")
  meta_key <- if (method == "unifrac_unweighted") "unifrac" else method
  dimred_name <- paste0("MDS_", if (method == "unifrac_unweighted") "unifrac" else method)

  if (meta_key %in% names(metadata(tse))) {
    mds <- cmdscale(as.dist(metadata(tse)[[meta_key]]), k = args$ncomponents, eig = TRUE)
    reducedDim(tse, dimred_name) <- mds$points
    metadata(tse)[[paste0(meta_key, "_eig")]] <- mds$eig
  } else {
    cat("Skipping MDS for", method, "- dissimilarity matrix not found in metadata\n")
  }
}

cat("Distance matrices:", paste(names(metadata(tse)), collapse = ", "), "\n")
cat("Reduced dimensions:", paste(reducedDimNames(tse), collapse = ", "), "\n")

# Save the modified TSE object
output_file <- file.path(args$output_dir, paste0("tse_", args$tse_type, "_with_beta.rds"))
saveRDS(tse, file = output_file)
cat("Saved modified TSE to:", output_file, "\n")

# Save beta diversity results separately
beta_results <- list(
  distance_matrices = metadata(tse)[names(metadata(tse))],
  reduced_dimensions = reducedDims(tse),
  metadata = build_metadata(args, tse, params = list(
    methods_used = distance_methods,
    niter = niter,
    ncomponents = args$ncomponents
  ))
)

results_file <- file.path(args$output_dir, paste0("beta_diversity_", args$tse_type, ".rds"))
saveRDS(beta_results, file = results_file)
cat("Saved beta diversity results to:", results_file, "\n")

cat("=== SUCCESS: Beta diversity calculation completed for", args$tse_type, "TSE! ===\n")
