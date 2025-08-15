#!/usr/bin/env Rscript

# Format DADA2 Taxonomy for QIIME2 Compatibility
# This script reformats the DADA2 taxonomy output to be compatible with QIIME2 visualization
#
# Arguments:
# 1. Input taxonomy TSV file (from DADA2 addSpecies)
# 2. Output taxonomy TSV file

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Function to generate a MD5 hash for a string
hash_string <- function(x) {
  # Create MD5 hash
  md5 <- digest::digest(x, algo = "md5", serialize = FALSE)
  return(md5)
}

# Read the input taxonomy file
message("Reading DADA2 taxonomy file...")
tax_df <- read.table(input_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# Add taxonomy rank prefixes and convert separators
add_prefixes <- function(taxon_string) {
  # Split the taxonomy string by semicolon and space
  taxa <- unlist(strsplit(taxon_string, "; "))
  
  # Assign rank prefixes based on position
  rank_prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
  
  # Format each taxon with its prefix
  formatted_taxa <- character(length = length(taxa))
  for (i in 1:length(taxa)) {
    if (i <= length(rank_prefixes)) {
      formatted_taxa[i] <- paste0(rank_prefixes[i], taxa[i])
    } else {
      # For any additional taxonomy levels beyond the standard 7
      formatted_taxa[i] <- paste0("u__", taxa[i]) # Use "u__" for unspecified rank
    }
  }
  
  # If less than 7 levels, fill with empty ranks
  if (length(formatted_taxa) < 7) {
    for (i in (length(formatted_taxa) + 1):7) {
      formatted_taxa <- c(formatted_taxa, rank_prefixes[i])
    }
  }
  
  # Join with semicolons (no spaces)
  return(paste(formatted_taxa, collapse = ";"))
}

message("Formatting taxonomy entries...")
# Apply the prefix function to each taxonomy entry
tax_df$Taxon <- sapply(tax_df$Taxon, add_prefixes)

# Hash feature IDs if requested
message("Hashing feature IDs...")
suppressPackageStartupMessages(library(digest))
tax_df$Feature_ID <- sapply(tax_df$Feature_ID, hash_string)

# Rename columns to match expected format
names(tax_df)[names(tax_df) == "Feature_ID"] <- "Feature ID"

# Make sure the column order matches exactly what QIIME2 expects
tax_df <- tax_df[, c("Feature ID", "Taxon")]

# Write the formatted taxonomy file
message("Writing formatted taxonomy file...")
write.table(tax_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

message("Taxonomy formatting completed.")
