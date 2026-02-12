#' Input Validation Functions for GWASTargetChase
#'
#' These functions provide input validation and error checking for the package.
#'
#' @name validation
#' @keywords internal
NULL

#' Validate GWAS input file
#'
#' Checks that the GWAS summary statistics file exists and contains required columns.
#'
#' @param gwas_file Path to GWAS summary statistics file
#' @param required_cols Character vector of required column names
#' @return TRUE if valid, throws error otherwise
#' @export
validate_gwas_input <- function(gwas_file,
                                 required_cols = c("chr", "ps", "p_wald")) {
  # Check file exists
  validate_file_exists(gwas_file)

  # Read header to check columns
  header <- tryCatch({
    data.table::fread(gwas_file, nrows = 1)
  }, error = function(e) {
    stop("Cannot read GWAS file: ", e$message, call. = FALSE)
  })

  # Check for empty file
  if (nrow(header) == 0) {
    full_data <- data.table::fread(gwas_file)
    if (nrow(full_data) == 0) {
      stop("GWAS file is empty or contains no data rows.", call. = FALSE)
    }
  }

  # Check required columns
  missing_cols <- setdiff(required_cols, names(header))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s) in GWAS file: ",
         paste(missing_cols, collapse = ", "),
         "\nRequired columns: ", paste(required_cols, collapse = ", "),
         call. = FALSE)
  }

  invisible(TRUE)
}

#' Check for significant SNPs
#'
#' Verifies that the GWAS data contains SNPs below the p-value threshold.
#'
#' @param gwas_file Path to GWAS summary statistics file
#' @param pval P-value threshold
#' @return Number of significant SNPs, or throws error/warning if none
#' @export
check_significant_snps <- function(gwas_file, pval = 5e-8) {
  gwas <- data.table::fread(gwas_file)

  if (!"p_wald" %in% names(gwas)) {
    stop("Column 'p_wald' not found in GWAS file.", call. = FALSE)
  }

  n_sig <- sum(gwas$p_wald <= pval, na.rm = TRUE)

  if (n_sig == 0) {
    stop("No significant SNPs found with p-value <= ", pval,
         ". Consider increasing the p-value threshold.",
         call. = FALSE)
  }

  message("Found ", n_sig, " SNPs with p-value <= ", pval)
  return(n_sig)
}

#' Validate file exists
#'
#' @param file_path Path to check
#' @return TRUE if exists, throws error otherwise
#' @export
validate_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path, call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate p-value threshold
#'
#' @param pval P-value threshold to validate
#' @return TRUE if valid, throws error otherwise
#' @export
validate_pvalue_threshold <- function(pval) {
  if (!is.numeric(pval)) {
    stop("P-value threshold must be numeric.", call. = FALSE)
  }

  if (length(pval) != 1) {
    stop("P-value threshold must be a single value.", call. = FALSE)
  }

  if (pval <= 0) {
    stop("P-value threshold must be greater than 0.", call. = FALSE)
  }

  if (pval > 1) {
    stop("P-value threshold must be less than or equal to 1.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Validate gene name format
#'
#' @param gene_name Gene name to validate
#' @return TRUE if valid
#' @export
validate_gene_name <- function(gene_name) {
  if (is.na(gene_name) || gene_name == "") {
    stop("Gene name cannot be empty or NA.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Normalize chromosome names
#'
#' Converts chromosome names to a standard format (1-22, X, Y).
#'
#' @param gwas_file Path to GWAS file or data.frame
#' @return Data frame with normalized chromosome names
#' @export
normalize_chromosomes <- function(gwas_file) {
  if (is.character(gwas_file)) {
    gwas <- data.table::fread(gwas_file)
  } else {
    gwas <- gwas_file
  }

  # Remove 'chr' prefix if present
  gwas$chr <- gsub("^chr", "", gwas$chr, ignore.case = TRUE)

  # Convert chromosome 23 to X if present
  gwas$chr <- gsub("^23$", "X", gwas$chr)

  # Convert chromosome 24 to Y if present
  gwas$chr <- gsub("^24$", "Y", gwas$chr)

  return(gwas)
}

#' Check for missing values in data
#'
#' @param file_path Path to data file
#' @param critical_cols Columns that must not have missing values
#' @return List with has_missing and columns_with_na
#' @export
check_missing_values <- function(file_path, critical_cols = c("chr", "ps", "p_wald")) {
  data <- data.table::fread(file_path)

  na_counts <- sapply(data, function(x) sum(is.na(x)))
  cols_with_na <- names(na_counts[na_counts > 0])

  result <- list(
    has_missing = length(cols_with_na) > 0,
    columns_with_na = cols_with_na,
    na_counts = na_counts[na_counts > 0]
  )

  # Check critical columns
  critical_with_na <- intersect(critical_cols, cols_with_na)
  if (length(critical_with_na) > 0) {
    warning("Critical columns have missing values: ",
            paste(critical_with_na, collapse = ", "),
            call. = FALSE)
  }

  return(result)
}

#' Validate numeric columns
#'
#' @param file_path Path to data file
#' @param cols Columns to check
#' @return TRUE if valid, throws warning for coercion
#' @export
validate_numeric_columns <- function(file_path, cols = c("ps", "p_wald")) {
  data <- data.table::fread(file_path)

  for (col in cols) {
    if (col %in% names(data)) {
      # Check if column can be converted to numeric
      test_numeric <- suppressWarnings(as.numeric(data[[col]]))
      n_na_before <- sum(is.na(data[[col]]))
      n_na_after <- sum(is.na(test_numeric))

      if (n_na_after > n_na_before) {
        warning("Column '", col, "' contains non-numeric values that will be coerced to NA.",
                call. = FALSE)
      }
    }
  }

  invisible(TRUE)
}

#' Validate output path
#'
#' @param output_path Path to output directory
#' @return TRUE if valid, throws error otherwise
#' @export
validate_output_path <- function(output_path) {
  # Normalize path
  output_path <- normalizePath(output_path, mustWork = FALSE)

  # Check directory exists
  if (!dir.exists(output_path)) {
    stop("Output directory does not exist: ", output_path, call. = FALSE)
  }

  # Check writability by attempting to create a temp file
  test_file <- file.path(output_path, paste0(".test_", Sys.getpid()))
  tryCatch({
    writeLines("test", test_file)
    unlink(test_file)
  }, error = function(e) {
    stop("Cannot write to output directory: ", output_path, call. = FALSE)
  })

  invisible(TRUE)
}

#' Validate GTF file format
#'
#' @param gtf_file Path to GTF file
#' @return TRUE if valid, throws error otherwise
#' @export
validate_gtf_file <- function(gtf_file) {
  validate_file_exists(gtf_file)

  # Try to read first few lines
  lines <- readLines(gtf_file, n = 20, warn = FALSE)

  # Remove comment lines
  data_lines <- lines[!grepl("^#", lines)]

  if (length(data_lines) == 0) {
    stop("GTF file contains no data (only comments).", call. = FALSE)
  }

  # Check GTF format (should have 9 tab-separated fields)
  first_data <- strsplit(data_lines[1], "\t")[[1]]

  if (length(first_data) < 9) {
    stop("Invalid GTF format: expected 9 tab-separated columns, found ",
         length(first_data), call. = FALSE)
  }

  # Check for gene_id attribute
  if (!grepl("gene_id", data_lines[1])) {
    stop("Invalid GTF format: 'gene_id' attribute not found.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Comprehensive input validation for gwasFollowup functions
#'
#' @param sumStats Path to GWAS summary statistics
#' @param felGTF Path to GTF file
#' @param pval P-value threshold
#' @param ResultsPath Output directory path
#' @return TRUE if all validations pass
#' @export
validate_gwas_followup_inputs <- function(sumStats, felGTF, pval, ResultsPath) {
  # Validate all inputs
  validate_file_exists(sumStats)
  validate_gwas_input(sumStats)
  validate_file_exists(felGTF)
  validate_gtf_file(felGTF)
  validate_pvalue_threshold(pval)
  validate_output_path(ResultsPath)

  # Check for significant SNPs
  check_significant_snps(sumStats, pval)

  # Check for missing values
  mv_result <- check_missing_values(sumStats)
  if (mv_result$has_missing) {
    message("Note: Input data contains missing values in columns: ",
            paste(mv_result$columns_with_na, collapse = ", "))
  }

  invisible(TRUE)
}
