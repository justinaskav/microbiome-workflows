#!/usr/bin/env Rscript

# DADA2 Taxonomy Assignment Script
# This script assigns taxonomy to ASV sequences using the DADA2 naive Bayesian classifier
#
# Arguments:
# 1. Representative sequences file path
# 2. Silva trainset file path
# 3. minBoot - Minimum bootstrap confidence threshold
# 4. cpus - Number of CPU cores to use
# 5. seed - Random seed for reproducibility

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
rep_seqs_file <- args[1]
trainset_file <- args[2]
minBoot <- as.numeric(args[3])
cpus <- as.numeric(args[4])
seed <- as.numeric(args[5])

# Set seed for reproducibility
set.seed(seed)

# Load required packages
suppressPackageStartupMessages(library(dada2))
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

# Read representative sequences
message("Reading representative sequences...")
if (grepl("\\.qza$", rep_seqs_file)) {
  # Try to use qiime2R if available
  if (requireNamespace("qiime2R", quietly = TRUE)) {
    message("Using qiime2R to read QIIME2 artifact...")
    seqs <- qiime2R::read_qza(rep_seqs_file)$data
  } else {
    # Fall back to manual extraction
    seqs <- extract_sequences_from_qza(rep_seqs_file)
  }
} else {
  # Handle regular FASTA file
  seqs <- Biostrings::readDNAStringSet(rep_seqs_file)
}

# Check if the trainset file is compressed and uncompress if needed
if (grepl("\\.gz$", trainset_file)) {
  message("Decompressing trainset file...")
  uncompressed_trainset <- gsub("\\.gz$", "", trainset_file)
  system(paste("gunzip -c", trainset_file, ">", uncompressed_trainset))
  trainset_file <- uncompressed_trainset
}

# Run assignTaxonomy using the Silva trainset
message("Running assignTaxonomy with Silva v138.2 reference...")
tax <- dada2::assignTaxonomy(seqs,
  trainset_file,
  multithread = cpus,
  minBoot = minBoot,
  tryRC = TRUE,
  verbose = TRUE
)

# Write the results
message("Saving taxonomy assignment results...")
saveRDS(tax, "taxonomy.rds")

# Create a TSV output format
write_taxonomy_tsv <- function(tax_table) {
  # Combine taxonomic levels into strings
  tax_strings <- apply(tax_table, 1, function(x) {
    paste(na.omit(x), collapse = "; ")
  })

  # Create a data frame for output
  tax_df <- data.frame(
    Feature_ID = names(tax_strings),
    Taxon = tax_strings,
    stringsAsFactors = FALSE
  )

  # Write to TSV
  write.table(tax_df, "taxonomy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
}

write_taxonomy_tsv(tax)

message("DADA2 taxonomy assignment completed.")
