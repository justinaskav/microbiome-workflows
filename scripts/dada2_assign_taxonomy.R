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

# Read the FASTA payload of a QIIME 2 artifact (.qza is a zip: <uuid>/data/dna-sequences.fasta)
read_qza_fasta <- function(qza_file) {
  members <- utils::unzip(qza_file, list = TRUE)$Name
  fasta <- grep("^[^/]+/data/.*\\.(fa|fasta|fna)$", members, value = TRUE)
  if (length(fasta) != 1) {
    stop("Expected exactly one FASTA file in ", qza_file, ", found ", length(fasta))
  }
  exdir <- tempfile("qza_")
  utils::unzip(qza_file, files = fasta, exdir = exdir, junkpaths = TRUE)
  Biostrings::readDNAStringSet(file.path(exdir, basename(fasta)))
}

# Read representative sequences
message("Reading representative sequences...")
if (grepl("\\.qza$", rep_seqs_file)) {
  seqs <- read_qza_fasta(rep_seqs_file)
} else {
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
