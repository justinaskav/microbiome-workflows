#! /usr/bin/env Rscript
library(maaslin3)
library(qiime2R)
library(argparse)
library(TreeSummarizedExperiment)

# Create a parser
parser <- ArgumentParser(description = "Run MaAsLin3 analysis using qiime2R::qza_to_tse with optional taxa agglomeration")

# Add arguments
parser$add_argument("--features", type = "character", required = TRUE, help = "Path to feature table QZA file (e.g., table.qza)")
parser$add_argument("--tree", type = "character", required = FALSE, default = NULL, help = "Path to rooted tree QZA file (e.g., rooted_tree.qza, optional)")
parser$add_argument("--taxonomy", type = "character", required = TRUE, help = "Path to taxonomy QZA file (e.g., taxonomy.qza)")
parser$add_argument("--metadata", type = "character", required = TRUE, help = "Path to metadata TSV file")
parser$add_argument("--formula", type = "character", required = TRUE, help = "MaAsLin3 formula string (e.g., 'Description' or '~Var1+Var2')")
parser$add_argument("--reference_levels", type = "character", required = FALSE, default = NULL, help = "Reference levels for factors, comma-separated, e.g., 'Factor1,RefLevel1;Factor2,RefLevel2'")
parser$add_argument("--taxa_level",
  type = "character", required = FALSE, default = "all",
  help = "Taxonomic level for agglomeration (e.g., 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', or 'all' for no agglomeration)"
)
parser$add_argument("--seed",
  type = "integer", required = FALSE, default = 42,
  help = "Random seed for reproducibility (default: 42)"
)

# Parse arguments
args <- parser$parse_args()
set.seed(args$seed) # Set seed for reproducibility

# Prepare arguments for qza_to_tse
tse_args <- list(
  features = args$features,
  taxonomy = args$taxonomy,
  metadata = args$metadata
)

if (!is.null(args$tree) && file.exists(args$tree)) {
  tse_args$tree <- args$tree
} else if (!is.null(args$tree)) {
  print(paste("Warning: Tree file not found:", args$tree, "- proceeding without tree."))
}

# Load data using qza_to_tse
tse <- do.call(qiime2R::qza_to_tse, tse_args)

# Prepare MaAsLin3 parameters
maaslin3_params <- list(
  input_data = tse,
  output = tolower(args$taxa_level), # MaAsLin3 creates this directory
  formula = args$formula
)

# Add reference levels if provided
if (!is.null(args$reference_levels) && nzchar(args$reference_levels)) {
  # Maaslin3 expects reference as a comma-separated string like "Factor1,RefLevel1,Factor2,RefLevel2"
  # The input format is "Factor1,RefLevel1;Factor2,RefLevel2" for easier parsing if multiple.
  # We need to convert ";" to "," if multiple are provided this way.
  # However, maaslin3 documentation for `reference` says:
  # "character string, e.g. 'X1,X1_ref;X2,X2_ref'. This defines the reference
  # levels for variables in the formula. Multiple variables are separated by
  #   ';'. Each variable and its reference level are separated by ','."
  # So the input format is already correct for MaAsLin3.
  maaslin3_params$reference <- args$reference_levels
  print(paste("Using reference levels:", args$reference_levels))
}

# Function to agglomerate TSE at a specific taxonomic level
agglomerate_tse <- function(tse, taxa_level) {
  if (taxa_level == "all") {
    return(tse) # Return original TSE without agglomeration
  }

  tse_agglomerated <- mia::agglomerateByRank(tse, rank = taxa_level)

  return(tse_agglomerated)
}

# Agglomerate TSE at specified taxonomic level
tse_processed <- agglomerate_tse(tse, args$taxa_level)
maaslin3_params$input_data <- tse_processed

# Run MaAsLin3
maaslin3_results <- do.call(maaslin3, maaslin3_params)
