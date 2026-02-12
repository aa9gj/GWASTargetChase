# Test input validation and edge cases for GWASTargetChase

context("Input Validation Tests")

# Helper function to create temporary test files
create_temp_gwas <- function(n_snps = 10, p_values = NULL, include_cols = TRUE) {
  if (is.null(p_values)) {
    p_values <- runif(n_snps, 1e-10, 0.5)
  }

  df <- data.frame(
    chr = sample(1:22, n_snps, replace = TRUE),
    rs = paste0("rs", sample(1000000:9999999, n_snps)),
    ps = sample(1000000:100000000, n_snps),
    p_wald = p_values,
    af = runif(n_snps, 0.01, 0.5),
    beta = rnorm(n_snps, 0, 0.3),
    se = abs(rnorm(n_snps, 0.05, 0.01)),
    stringsAsFactors = FALSE
  )

  if (!include_cols) {
    # Remove required columns for testing
    df$p_wald <- NULL
  }

  temp_file <- tempfile(fileext = ".tsv")
  write.table(df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
  return(temp_file)
}

create_temp_gtf <- function(n_genes = 5) {
  genes <- paste0("GENE", 1:n_genes)

  lines <- c(
    '##description: Test GTF file',
    '##provider: Test'
  )

  for (i in 1:n_genes) {
    start <- sample(1000000:50000000, 1)
    end <- start + sample(10000:100000, 1)
    chr <- sample(1:22, 1)
    lines <- c(lines, sprintf(
      '%s\tensembl\tgene\t%d\t%d\t.\t+\t.\tgene_id "ENSG%08d"; gene_name "%s"; gene_biotype "protein_coding";',
      chr, start, end, i, genes[i]
    ))
  }

  temp_file <- tempfile(fileext = ".gtf")
  writeLines(lines, temp_file)
  return(temp_file)
}

# ==============================================================================
# Test: validate_gwas_input function
# ==============================================================================

test_that("validate_gwas_input detects missing required columns", {
  # Create GWAS without p_wald column
  temp_gwas <- create_temp_gwas(n_snps = 10, include_cols = FALSE)

  expect_error(
    validate_gwas_input(temp_gwas),
    regexp = "Missing required column"
  )

  unlink(temp_gwas)
})

test_that("validate_gwas_input accepts valid GWAS file", {
  temp_gwas <- create_temp_gwas(n_snps = 10)

  expect_silent(validate_gwas_input(temp_gwas))

  unlink(temp_gwas)
})

# ==============================================================================
# Test: Edge case - Empty GWAS file
# ==============================================================================

test_that("Empty GWAS file is handled gracefully", {
  temp_file <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(chr = character(), rs = character(), ps = numeric(),
               p_wald = numeric(), stringsAsFactors = FALSE),
    temp_file, sep = "\t", row.names = FALSE, quote = FALSE
  )

  expect_error(
    validate_gwas_input(temp_file),
    regexp = "empty|no data|zero rows"
  )

  unlink(temp_file)
})

# ==============================================================================
# Test: Edge case - No significant SNPs
# ==============================================================================

test_that("No significant SNPs produces appropriate message", {
  # Create GWAS with all p-values > 0.05
  temp_gwas <- create_temp_gwas(n_snps = 100, p_values = runif(100, 0.1, 1))
  temp_gtf <- create_temp_gtf(n_genes = 10)

  # The function should handle this gracefully
  expect_message(
    tryCatch(
      check_significant_snps(temp_gwas, pval = 1e-10),
      error = function(e) message(e$message)
    ),
    regexp = "no significant|p-value|consider"
  )

  unlink(c(temp_gwas, temp_gtf))
})

# ==============================================================================
# Test: Edge case - Invalid file paths
# ==============================================================================

test_that("Invalid file paths are caught", {
  expect_error(
    validate_file_exists("/nonexistent/path/to/file.tsv"),
    regexp = "does not exist|not found"
  )
})

# ==============================================================================
# Test: Edge case - Extreme p-value thresholds
# ==============================================================================

test_that("Extreme p-value thresholds are validated", {
  # p-value of 0 should be invalid

expect_error(
    validate_pvalue_threshold(0),
    regexp = "p-value|must be|positive|greater than"
  )

  # p-value > 1 should be invalid
  expect_error(
    validate_pvalue_threshold(1.5),
    regexp = "p-value|must be|less than|between"
  )

  # Negative p-value should be invalid
  expect_error(
    validate_pvalue_threshold(-0.05),
    regexp = "p-value|positive|negative"
  )
})

test_that("Valid p-value thresholds are accepted", {
  expect_silent(validate_pvalue_threshold(0.05))
  expect_silent(validate_pvalue_threshold(5e-8))
  expect_silent(validate_pvalue_threshold(1e-300))
})

# ==============================================================================
# Test: Edge case - Special characters in data
# ==============================================================================

test_that("Special characters in gene names are handled", {
  # Some gene names have special characters like hyphens
  gene_names <- c("HLA-A", "HLA-B", "TP53", "BRCA1/BRCA2", "Gene.Name")

  for (gene in gene_names) {
    expect_silent(validate_gene_name(gene))
  }
})

# ==============================================================================
# Test: Edge case - Chromosome naming conventions
# ==============================================================================

test_that("Different chromosome naming conventions are handled", {
  # Create GWAS with different chr formats
  formats <- list(
    numeric = c(1, 2, 22),
    chr_prefix = c("chr1", "chr2", "chrX"),
    mixed = c("1", "chr2", "X")
  )

  for (fmt_name in names(formats)) {
    df <- data.frame(
      chr = formats[[fmt_name]],
      rs = paste0("rs", 1:3),
      ps = c(1000, 2000, 3000),
      p_wald = c(1e-10, 1e-8, 1e-6)
    )
    temp_file <- tempfile(fileext = ".tsv")
    write.table(df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

    result <- normalize_chromosomes(temp_file)
    expect_true(all(result$chr %in% c(1:22, "X", "Y")))

    unlink(temp_file)
  }
})

# ==============================================================================
# Test: Edge case - Missing values (NA) in data
# ==============================================================================

test_that("Missing values in GWAS data are detected", {
  df <- data.frame(
    chr = c(1, NA, 3),
    rs = c("rs1", "rs2", NA),
    ps = c(1000, 2000, 3000),
    p_wald = c(1e-10, NA, 1e-6)
  )
  temp_file <- tempfile(fileext = ".tsv")
  write.table(df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- check_missing_values(temp_file)
  expect_true(result$has_missing)
  expect_true("p_wald" %in% result$columns_with_na)

  unlink(temp_file)
})

# ==============================================================================
# Test: Stress test - Large dataset handling
# ==============================================================================

test_that("Large GWAS files can be processed", {
  skip_on_cran()  # Skip on CRAN due to time/memory

  # Create a moderately large GWAS file (10,000 SNPs)
  n_snps <- 10000
  temp_gwas <- create_temp_gwas(n_snps = n_snps)

  # Should complete without error
  expect_silent({
    result <- validate_gwas_input(temp_gwas)
  })

  unlink(temp_gwas)
})

# ==============================================================================
# Test: Data format validation
# ==============================================================================

test_that("Numeric columns contain valid numbers", {
  df <- data.frame(
    chr = c(1, 2, 3),
    rs = c("rs1", "rs2", "rs3"),
    ps = c("1000", "invalid", "3000"),  # Invalid position
    p_wald = c(1e-10, 1e-8, 1e-6)
  )
  temp_file <- tempfile(fileext = ".tsv")
  write.table(df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  expect_warning(
    validate_numeric_columns(temp_file, cols = c("ps", "p_wald")),
    regexp = "non-numeric|invalid|coerced"
  )

  unlink(temp_file)
})

# ==============================================================================
# Test: Output directory validation
# ==============================================================================

test_that("Output directory must exist and be writable", {
  expect_error(
    validate_output_path("/nonexistent/directory/"),
    regexp = "does not exist|not found|cannot write"
  )

  # Valid temp directory should work
  temp_dir <- tempdir()
  expect_silent(validate_output_path(temp_dir))
})

# ==============================================================================
# Test: GTF file validation
# ==============================================================================

test_that("Invalid GTF files are detected", {
  # Create invalid GTF file
  temp_file <- tempfile(fileext = ".gtf")
  writeLines(c("not", "a", "valid", "gtf"), temp_file)

  expect_error(
    validate_gtf_file(temp_file),
    regexp = "invalid|GTF|format|parse"
  )

  unlink(temp_file)
})

test_that("Valid GTF files are accepted", {
  temp_gtf <- create_temp_gtf(n_genes = 5)

  expect_silent(validate_gtf_file(temp_gtf))

  unlink(temp_gtf)
})
