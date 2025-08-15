#!/usr/bin/env Rscript

# DADA2 Add Species Script
# This script adds species-level assignments to taxonomy results using exact matching
# 
# Arguments:
# 1. Taxonomy RDS file path (output from assignTaxonomy)
# 2. Representative sequences file path
# 3. Species assignment reference file path
# 4. allowMultiple - Whether to allow multiple species assignments
# 5. verbose - Whether to show verbose output
# 6. seed - Random seed for reproducibility

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
taxonomy_rds_file <- args[1]
rep_seqs_file <- args[2]
species_file <- args[3]
allowMultiple <- as.integer(args[4])
verbose <- as.logical(args[5])
seed <- as.numeric(args[6])

# Set seed for reproducibility
set.seed(seed)

# Load required packages
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(Biostrings))
if (grepl("\\.qza$", rep_seqs_file)) {
  suppressPackageStartupMessages(library(qiime2R))
}

# Function to handle QIIME2 artifacts if qiime2R is not available
extract_sequences_from_qza <- function(qza_file) {
  # Create a temporary directory
  temp_dir <- tempdir()
  
  # Unzip the QZA file (which is a zip archive)
  message("Extracting QIIME2 artifact...")
  system(paste("unzip -q -o", qza_file, "-d", temp_dir))
  
  # The data in a .qza file is typically in the data directory
  data_dir <- file.path(temp_dir, "data")
  
  # Find the FASTA file (should be in the data directory)
  fasta_files <- list.files(path = data_dir, pattern = "\\.fa$|\\.fasta$|\\.fna$", recursive = TRUE, full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No FASTA files found in the QIIME2 artifact")
  }
  
  # Read the first FASTA file found
  seqs <- Biostrings::readDNAStringSet(fasta_files[1])
  
  return(seqs)
}

# Read the taxonomy from assignTaxonomy
message("Loading taxonomy assignments...")
tax <- readRDS(taxonomy_rds_file)

# Check if the species file is compressed and uncompress if needed
if (grepl("\\.gz$", species_file)) {
  message("Decompressing species file...")
  uncompressed_species <- gsub("\\.gz$", "", species_file)
  system(paste("gunzip -c", species_file, ">", uncompressed_species))
  species_file <- uncompressed_species
}

# Add species assignment
message("Running addSpecies with Silva v138.2 species reference...")
tax_with_species <- dada2::addSpecies(tax, 
                               species_file, 
                               allowMultiple=allowMultiple, 
                               verbose=verbose)

# Save the results as RDS
message("Saving taxonomy with species assignments...")
saveRDS(tax_with_species, "taxonomy_with_species.rds")

# Create a TSV output format
write_taxonomy_tsv <- function(tax_table) {
  # Combine taxonomic levels into strings
  tax_strings <- apply(tax_table, 1, function(x) {
    paste(na.omit(x), collapse="; ")
  })
  
  # Create a data frame for output
  tax_df <- data.frame(
    Feature_ID = names(tax_strings),
    Taxon = tax_strings,
    stringsAsFactors = FALSE
  )
  
  # Write to TSV
  write.table(tax_df, "taxonomy_with_species.tsv", sep="\t", quote=FALSE, row.names=FALSE)
}

write_taxonomy_tsv(tax_with_species)

message("DADA2 species assignment completed.")
