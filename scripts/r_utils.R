# Shared CLI utilities for calc_*.R Nextflow scripts
# Source this at the top of each calc_ script.

INCERTAE_SEDIS_PATTERN <- "Incertae.Sedis"

#' Check if taxa names match Incertae Sedis pattern
#'
#' @param taxa_names Character vector of taxa names
#' @param pattern Regex pattern to match (dot matches both underscore and space)
#' @return Logical vector
is_incertae_sedis <- function(taxa_names, pattern = INCERTAE_SEDIS_PATTERN) {
  grepl(pattern, taxa_names, ignore.case = TRUE)
}

#' Filter Incertae Sedis taxa from a TSE object
#'
#' @param tse TreeSummarizedExperiment object (already aggregated to desired rank)
#' @param rank Taxonomic rank column name in rowData
#' @return Filtered TSE object
filter_incertae_sedis <- function(tse, rank) {
  row_names <- rownames(tse)
  keep <- !is_incertae_sedis(row_names)
  if (sum(!keep) > 0) {
    cat("Filtering", sum(!keep), "Incertae Sedis taxa\n")
  }
  tse[keep, ]
}

#' Log parsed arguments compactly
#'
#' @param args Parsed argparse list
#' @param show Character vector of argument names to display (default: all)
log_args <- function(args, show = NULL) {
  if (is.null(show)) show <- names(args)
  cat("Arguments:\n")
  for (nm in show) {
    val <- args[[nm]]
    if (is.null(val)) val <- "NULL"
    else if (is.logical(val)) val <- ifelse(val, "TRUE", "FALSE")
    cat("  ", nm, ": ", as.character(val), "\n", sep = "")
  }
}

#' Auto-detect TSE type from filename
#'
#' Expects filenames like "tse_<type>.rds"
#' @param tse_path Path to the TSE RDS file
#' @param tse_type Explicit type (returned as-is if not NULL)
#' @return Character string with TSE type
detect_tse_type <- function(tse_path, tse_type = NULL) {
  if (!is.null(tse_type)) return(tse_type)
  tse_filename <- basename(tse_path)
  if (grepl("^tse_(.+)\\.rds$", tse_filename)) {
    detected <- gsub("^tse_(.+)\\.rds$", "\\1", tse_filename)
    cat("Auto-detected TSE type:", detected, "\n")
    return(detected)
  }
  cat("WARNING: Could not auto-detect TSE type, using 'unknown'\n")
  "unknown"
}

#' Load and validate a TSE object
#'
#' @param tse_path Path to the TSE RDS file
#' @return TreeSummarizedExperiment object
load_tse <- function(tse_path) {
  if (!file.exists(tse_path)) stop("Input file does not exist: ", tse_path)
  tse <- readRDS(tse_path)
  cat("Loaded TSE:", nrow(tse), "features x", ncol(tse), "samples\n")
  tse
}

#' Filter samples with NA values in formula variables
#'
#' @param tse TreeSummarizedExperiment object
#' @param formula_str Formula string (with or without ~)
#' @return Filtered TSE
filter_na_samples <- function(tse, formula_str) {
  formula_check <- if (grepl("^~", formula_str)) formula_str else paste0("~", formula_str)
  formula_vars <- all.vars(as.formula(formula_check))

  metadata_df <- as.data.frame(colData(tse))
  existing_vars <- formula_vars[formula_vars %in% colnames(metadata_df)]
  missing_vars <- formula_vars[!formula_vars %in% colnames(metadata_df)]

  if (length(missing_vars) > 0) {
    cat("WARNING: Variables not found in metadata:", paste(missing_vars, collapse = ", "), "\n")
  }

  if (length(existing_vars) > 0) {
    na_mask <- complete.cases(metadata_df[, existing_vars, drop = FALSE])
    n_removed <- sum(!na_mask)
    if (n_removed > 0) {
      cat("Removing", n_removed, "samples with NA values in formula variables\n")
      tse <- tse[, na_mask]
      cat("TSE after NA filtering:", nrow(tse), "features x", ncol(tse), "samples\n")
    }
  }
  tse
}

#' Filter samples by valid group values
#'
#' @param tse TreeSummarizedExperiment object
#' @param valid_groups_csv Comma-separated string of valid group values, or NULL
#' @return Filtered TSE
filter_by_group <- function(tse, valid_groups_csv) {
  if (is.null(valid_groups_csv)) return(tse)

  valid_groups <- trimws(strsplit(valid_groups_csv, ",")[[1]])
  cat("Valid groups:", paste(valid_groups, collapse = ", "), "\n")

  if (!"Group" %in% colnames(colData(tse))) {
    warning("Group column not found in metadata, skipping group filtering")
    return(tse)
  }

  current_groups <- as.character(colData(tse)$Group)
  group_mask <- current_groups %in% valid_groups
  n_filtered <- sum(!group_mask, na.rm = TRUE)

  if (n_filtered > 0) {
    removed_groups <- unique(current_groups[!group_mask])
    removed_groups <- removed_groups[!is.na(removed_groups)]
    cat("Filtering out", n_filtered, "samples (removed groups:", paste(removed_groups, collapse = ", "), ")\n")
    tse <- tse[, group_mask]
  }
  tse
}

#' Drop unused factor levels from TSE colData
#'
#' @param tse TreeSummarizedExperiment object
#' @return TSE with cleaned factor levels
drop_empty_levels <- function(tse) {
  metadata_df <- as.data.frame(colData(tse))
  factor_cols <- names(metadata_df)[sapply(metadata_df, is.factor)]
  changed <- character(0)
  for (col in factor_cols) {
    old_n <- length(levels(metadata_df[[col]]))
    metadata_df[[col]] <- droplevels(metadata_df[[col]])
    new_n <- length(levels(metadata_df[[col]]))
    if (old_n != new_n) changed <- c(changed, paste0(col, ": ", old_n, "->", new_n))
  }
  colData(tse) <- S4Vectors::DataFrame(metadata_df)
  if (length(changed) > 0) cat("Dropped levels:", paste(changed, collapse = ", "), "\n")
  tse
}

#' Build standardized result metadata
#'
#' @param args Parsed argparse list
#' @param tse TreeSummarizedExperiment object
#' @param params Named list of analysis-specific parameters
#' @return Named list with standardized metadata structure
build_metadata <- function(args, tse, params = list()) {
  list(
    tse_type = args$tse_type,
    timestamp = Sys.time(),
    seed = args$seed,
    n_features = nrow(tse),
    n_samples = ncol(tse),
    analysis_params = params
  )
}

#' Create an NA analysis report for metadata variables
#'
#' @param tse A TreeSummarizedExperiment object
#' @param variables Character vector of column names to check.
#'   Only columns present in colData(tse) are checked.
#' @param threshold Percentage threshold for flagging high-NA variables (default: 20).
#' @param verbose Print report to console (default: TRUE).
#' @return A list: na_report (data.frame), high_na_vars (character),
#'   na_threshold (numeric), analysis_date (Date).
create_na_report <- function(tse, variables, threshold = 20, verbose = TRUE) {
  vars_present <- intersect(variables, colnames(SummarizedExperiment::colData(tse)))
  metadata_df <- as.data.frame(SummarizedExperiment::colData(tse)[, vars_present, drop = FALSE])

  na_counts <- vapply(metadata_df, function(x) sum(is.na(x)), integer(1))
  na_pct <- (na_counts / nrow(metadata_df)) * 100

  na_report <- data.frame(
    Variable = names(na_counts),
    NACount = na_counts,
    NAPercentage = round(na_pct, 2),
    row.names = NULL
  )
  na_report <- na_report[order(-na_report$NAPercentage), ]

  high_na_vars <- na_report$Variable[na_report$NAPercentage > threshold]

  if (verbose) {
    cat("\n--- Missing Value (NA) Report ---\n")
    print(na_report)
    cat("--------------------------------\n")
    if (length(high_na_vars) > 0) {
      warning(
        "High proportion (>", threshold, "%) of NA values in: ",
        paste(high_na_vars, collapse = ", ")
      )
    }
  }

  list(
    na_report = na_report,
    high_na_vars = high_na_vars,
    na_threshold = threshold,
    analysis_date = Sys.Date()
  )
}
