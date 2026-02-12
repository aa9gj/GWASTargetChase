# Stress tests for GWASTargetChase
# These tests evaluate performance with large datasets and edge conditions

context("Stress Tests")

# ==============================================================================
# Helper functions for generating large test datasets
# ==============================================================================

generate_large_gwas <- function(n_snps = 100000, n_significant = 100) {
  # Generate a large GWAS dataset
  set.seed(42)

  # Generate p-values with some significant ones
  p_values <- c(
    runif(n_significant, 1e-15, 5e-8),  # Significant
    runif(n_snps - n_significant, 5e-8, 1)  # Non-significant
  )
  p_values <- sample(p_values)  # Shuffle

  df <- data.frame(
    chr = sample(c(1:22, "X"), n_snps, replace = TRUE),
    rs = paste0("rs", 1:n_snps),
    ps = sample(1:250000000, n_snps),
    p_wald = p_values,
    af = runif(n_snps, 0.01, 0.99),
    beta = rnorm(n_snps, 0, 0.5),
    se = abs(rnorm(n_snps, 0.05, 0.02)),
    stringsAsFactors = FALSE
  )

  temp_file <- tempfile(fileext = ".tsv")
  data.table::fwrite(df, temp_file, sep = "\t")
  return(list(file = temp_file, n_snps = n_snps, n_sig = n_significant))
}

generate_large_gtf <- function(n_genes = 20000) {
  set.seed(42)

  # Common gene name prefixes for realistic names
  prefixes <- c("GENE", "LOC", "LINC", "MIR", "SNORD", "RPL", "RPS", "MT-",
                "HLA-", "OR", "ZNF", "FAM", "KIAA", "C", "TMEM", "SLC")

  lines <- c(
    '##description: Stress test GTF file',
    '##provider: GWASTargetChase'
  )

  for (i in 1:n_genes) {
    chr <- sample(c(1:22, "X", "Y"), 1)
    start <- sample(1:200000000, 1)
    end <- start + sample(1000:500000, 1)
    strand <- sample(c("+", "-"), 1)
    gene_name <- paste0(sample(prefixes, 1), i)

    # Gene entry
    lines <- c(lines, sprintf(
      '%s\tensembl\tgene\t%d\t%d\t.\t%s\t.\tgene_id "ENSG%011d"; gene_name "%s"; gene_biotype "protein_coding";',
      chr, start, end, strand, i, gene_name
    ))

    # Add transcript and exon entries for some genes
    if (i %% 10 == 0) {
      lines <- c(lines, sprintf(
        '%s\tensembl\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "ENSG%011d"; gene_name "%s"; transcript_id "ENST%011d"; gene_biotype "protein_coding";',
        chr, start, end, strand, i, gene_name, i
      ))
    }
  }

  temp_file <- tempfile(fileext = ".gtf")
  writeLines(lines, temp_file)
  return(list(file = temp_file, n_genes = n_genes))
}

# ==============================================================================
# Memory and Performance Tests
# ==============================================================================

test_that("Large GWAS file validation completes in reasonable time", {
  skip_on_cran()
  skip_if_not(interactive(), "Skipping long-running stress test in non-interactive mode")

  gwas_data <- generate_large_gwas(n_snps = 50000, n_significant = 50)

  start_time <- Sys.time()

  expect_silent({
    validate_gwas_input(gwas_data$file)
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Should complete in under 10 seconds
  expect_lt(elapsed, 10)

  message("Validated ", gwas_data$n_snps, " SNPs in ", round(elapsed, 2), " seconds")

  unlink(gwas_data$file)
})

test_that("Large GTF file can be validated", {
  skip_on_cran()
  skip_if_not(interactive(), "Skipping long-running stress test in non-interactive mode")

  gtf_data <- generate_large_gtf(n_genes = 5000)

  start_time <- Sys.time()

  expect_silent({
    validate_gtf_file(gtf_data$file)
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  message("Validated GTF with ", gtf_data$n_genes, " genes in ", round(elapsed, 2), " seconds")

  unlink(gtf_data$file)
})

test_that("Chromosome normalization handles large datasets", {
  skip_on_cran()

  gwas_data <- generate_large_gwas(n_snps = 10000)

  start_time <- Sys.time()

  result <- normalize_chromosomes(gwas_data$file)

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  expect_equal(nrow(result), 10000)
  expect_lt(elapsed, 5)

  unlink(gwas_data$file)
})

# ==============================================================================
# Edge Case Stress Tests
# ==============================================================================

test_that("All SNPs significant is handled correctly", {
  skip_on_cran()

  # Create GWAS where every SNP is significant
  n_snps <- 1000
  df <- data.frame(
    chr = sample(1:22, n_snps, replace = TRUE),
    rs = paste0("rs", 1:n_snps),
    ps = sample(1:250000000, n_snps),
    p_wald = runif(n_snps, 1e-50, 1e-10),  # All very significant
    stringsAsFactors = FALSE
  )

  temp_file <- tempfile(fileext = ".tsv")
  write.table(df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  n_sig <- check_significant_snps(temp_file, pval = 5e-8)
  expect_equal(n_sig, n_snps)

  unlink(temp_file)
})

test_that("Single SNP GWAS file is handled", {
  df <- data.frame(
    chr = 1,
    rs = "rs12345",
    ps = 50000000,
    p_wald = 1e-10
  )

  temp_file <- tempfile(fileext = ".tsv")
  write.table(df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  expect_silent(validate_gwas_input(temp_file))

  n_sig <- check_significant_snps(temp_file, pval = 5e-8)
  expect_equal(n_sig, 1)

  unlink(temp_file)
})

test_that("Very small p-values are handled correctly", {
  # Test with extremely small p-values (underflow potential)
  df <- data.frame(
    chr = c(1, 2, 3),
    rs = c("rs1", "rs2", "rs3"),
    ps = c(1000, 2000, 3000),
    p_wald = c(1e-300, 1e-200, 1e-100)
  )

  temp_file <- tempfile(fileext = ".tsv")
  write.table(df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  n_sig <- check_significant_snps(temp_file, pval = 5e-8)
  expect_equal(n_sig, 3)

  unlink(temp_file)
})

test_that("Many overlapping genomic regions are handled", {
  skip_on_cran()

  # Create many SNPs in overlapping regions
  n_snps <- 1000
  df <- data.frame(
    chr = rep(1, n_snps),  # All on chromosome 1
    rs = paste0("rs", 1:n_snps),
    ps = seq(50000000, 51000000, length.out = n_snps),  # 1Mb region
    p_wald = runif(n_snps, 1e-15, 1e-8)
  )

  temp_file <- tempfile(fileext = ".tsv")
  write.table(df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  expect_silent(validate_gwas_input(temp_file))

  unlink(temp_file)
})

# ==============================================================================
# Memory Efficiency Tests
# ==============================================================================

test_that("Validation does not load entire file into memory when not needed", {
  skip_on_cran()
  skip_if_not(interactive(), "Skipping memory test")

  gwas_data <- generate_large_gwas(n_snps = 100000)

  # Get baseline memory
  gc()
  mem_before <- sum(gc()[, 2])

  # Validate file
  validate_gwas_input(gwas_data$file)

  # Check memory didn't increase dramatically
  gc()
  mem_after <- sum(gc()[, 2])

  # Memory increase should be minimal (less than 100MB)
  mem_increase <- mem_after - mem_before
  expect_lt(mem_increase, 100)

  unlink(gwas_data$file)
})

# ==============================================================================
# Concurrent/Batch Processing Tests
# ==============================================================================

test_that("Multiple files can be validated in sequence", {
  skip_on_cran()

  n_files <- 10
  temp_files <- character(n_files)

  # Create multiple GWAS files
  for (i in 1:n_files) {
    df <- data.frame(
      chr = sample(1:22, 100, replace = TRUE),
      rs = paste0("rs", ((i-1)*100 + 1):(i*100)),
      ps = sample(1:250000000, 100),
      p_wald = runif(100, 1e-10, 0.5)
    )
    temp_files[i] <- tempfile(fileext = ".tsv")
    write.table(df, temp_files[i], sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # Validate all files
  start_time <- Sys.time()

  for (f in temp_files) {
    expect_silent(validate_gwas_input(f))
  }

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  message("Validated ", n_files, " files in ", round(elapsed, 2), " seconds")

  # Cleanup
  sapply(temp_files, unlink)
})

# ==============================================================================
# Data Integrity Tests
# ==============================================================================

test_that("Unicode and special characters in file paths are handled", {
  skip_on_cran()
  skip_on_os("windows")  # Windows has different path handling

  # Create file with unicode in path (if supported)
  temp_dir <- file.path(tempdir(), "test_ü_path")
  dir.create(temp_dir, showWarnings = FALSE)

  df <- data.frame(
    chr = 1,
    rs = "rs12345",
    ps = 50000000,
    p_wald = 1e-10
  )

  temp_file <- file.path(temp_dir, "test_file.tsv")
  write.table(df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  if (file.exists(temp_file)) {
    expect_silent(validate_gwas_input(temp_file))
    unlink(temp_file)
  }

  unlink(temp_dir, recursive = TRUE)
})

test_that("Files with different line endings are handled", {
  df <- data.frame(
    chr = c(1, 2),
    rs = c("rs1", "rs2"),
    ps = c(1000, 2000),
    p_wald = c(1e-10, 1e-8)
  )

  # Test Unix line endings (LF)
  temp_unix <- tempfile(fileext = ".tsv")
  write.table(df, temp_unix, sep = "\t", row.names = FALSE, quote = FALSE, eol = "\n")
  expect_silent(validate_gwas_input(temp_unix))
  unlink(temp_unix)

  # Test Windows line endings (CRLF)
  temp_windows <- tempfile(fileext = ".tsv")
  write.table(df, temp_windows, sep = "\t", row.names = FALSE, quote = FALSE, eol = "\r\n")
  expect_silent(validate_gwas_input(temp_windows))
  unlink(temp_windows)
})
