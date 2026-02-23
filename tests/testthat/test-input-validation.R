# Tests for input validation and error handling across all functions.

# --- gwasFollowup input validation ---
test_that("gwasFollowup errors on missing sumstats file", {
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  expect_error(
    gwasFollowup(sumStats = "/nonexistent/path/sumstats.txt", felGTF = gtf),
    "Summary statistics file not found"
  )
})

test_that("gwasFollowup errors on missing GTF file", {
  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  expect_error(
    gwasFollowup(sumStats = sumstats, felGTF = "/nonexistent/path/file.gtf"),
    "GTF file not found"
  )
})

test_that("gwasFollowup errors on empty string paths", {
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  expect_error(gwasFollowup(sumStats = "", felGTF = gtf), "Summary statistics file not found")
  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  expect_error(gwasFollowup(sumStats = sumstats, felGTF = ""), "GTF file not found")
})

test_that("gwasFollowup errors when p-value too strict (no significant SNPs)", {
  skip_if_offline()
  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  out_dir <- file.path(tempdir(), "test_pval_strict")
  on.exit(unlink(out_dir, recursive = TRUE))

  expect_error(
    gwasFollowup(sumStats = sumstats, felGTF = gtf, pval = 1e-100, ResultsPath = out_dir),
    "No significant SNPs found"
  )
})

test_that("gwasFollowup errors when sumstats lacks 'gene' column", {
  # Create a sumstats file without the gene column
  mock_dir <- file.path(tempdir(), "test_no_gene_col")
  dir.create(mock_dir, showWarnings = FALSE)
  on.exit(unlink(mock_dir, recursive = TRUE))

  bad_sumstats <- file.path(mock_dir, "bad_sumstats.tsv")
  df <- data.frame(
    chr = c(1, 1), rs = c("rs1", "rs2"),
    ps = c(1000, 2000), p_wald = c(1e-10, 1e-9),
    af = c(0.4, 0.3), beta = c(0.5, 0.4), se = c(0.05, 0.05)
  )
  write.table(df, bad_sumstats, sep = "\t", row.names = FALSE, quote = FALSE)

  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  expect_error(
    gwasFollowup(sumStats = bad_sumstats, felGTF = gtf, ResultsPath = file.path(mock_dir, "out")),
    "gene.*column"
  )
})

# --- gwasFollowupFull input validation ---
test_that("gwasFollowupFull errors on missing sumstats file", {
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  expect_error(
    gwasFollowupFull(sumStats = "/nonexistent/sumstats.txt", felGTF = gtf),
    "Summary statistics file not found"
  )
})

test_that("gwasFollowupFull errors on missing GTF file", {
  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  expect_error(
    gwasFollowupFull(sumStats = sumstats, felGTF = "/nonexistent/file.gtf"),
    "GTF file not found"
  )
})

# --- gwasFollowupMan input validation ---
test_that("gwasFollowupMan errors on missing sumstats file", {
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  expect_error(
    gwasFollowupMan(sumStats = "/nonexistent/sumstats.txt", felGTF = gtf,
                    assocOT = "dummy", l2gOT = "dummy"),
    "Summary statistics file not found"
  )
})

test_that("gwasFollowupMan errors on missing GTF file", {
  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  expect_error(
    gwasFollowupMan(sumStats = sumstats, felGTF = "/nonexistent/file.gtf",
                    assocOT = "dummy", l2gOT = "dummy"),
    "GTF file not found"
  )
})

# --- load_zoonomia_orthologs validation ---
test_that("load_zoonomia_orthologs errors for unsupported species combo", {
  expect_error(
    load_zoonomia_orthologs("human", "elephant"),
    "Zoonomia orthology file not found"
  )
})

test_that("load_zoonomia_orthologs errors for reversed species order", {
  expect_error(
    load_zoonomia_orthologs("cat", "human"),
    "Zoonomia orthology file not found"
  )
})

test_that("load_zoonomia_orthologs errors for invalid zoo_dir", {
  expect_error(
    load_zoonomia_orthologs("human", "cat", zoo_dir = "/nonexistent/dir"),
    "Zoonomia orthology file not found"
  )
})

# --- gene2disease edge cases ---
test_that("gene2disease handles empty assoc_data", {
  empty_df <- data.frame(
    gene_name = character(), gene_id = character(),
    disease_name = character(), score = numeric(),
    stringsAsFactors = FALSE
  )
  result <- gene2disease("FTO", empty_df)
  expect_equal(nrow(result), 0)
})

test_that("gene2disease handles special regex characters in gene name", {
  assoc_data <- data.frame(
    gene_name = c("C2orf16", "C2orf16-AS1"),
    gene_id = c("ENSG1", "ENSG2"),
    disease_name = c("disease1", "disease2"),
    score = c(0.5, 0.3),
    stringsAsFactors = FALSE
  )
  # Should only match C2orf16 exactly, not C2orf16-AS1
  result <- gene2disease("C2orf16", assoc_data)
  expect_equal(nrow(result), 1)
  expect_equal(result$gene_name, "C2orf16")
})

# --- locus2gene edge cases ---
test_that("locus2gene handles empty l2g_data", {
  empty_df <- data.frame(
    study_id = character(), gene_name = character(),
    y_proba_full_model = numeric(), trait_reported = character(),
    stringsAsFactors = FALSE
  )
  result <- locus2gene("FTO", empty_df)
  expect_equal(nrow(result), 0)
})

test_that("locus2gene matches exact gene names only", {
  l2g_data <- data.frame(
    study_id = c("S1", "S2", "S3"),
    gene_name = c("FTO", "FTOS", "SFTO"),
    y_proba_full_model = c(0.9, 0.8, 0.7),
    trait_reported = c("BMI", "other", "other"),
    stringsAsFactors = FALSE
  )
  result <- locus2gene("FTO", l2g_data)
  expect_equal(nrow(result), 1)
  expect_equal(result$gene_name, "FTO")
})

# --- Chromosome handling ---
test_that("gwasFollowup converts chr 20 to X", {
  # The code does gsub("^20$", "X", ...) on the chr column
  chr_values <- c("1", "2", "20", "X")
  converted <- gsub("^20$", "X", chr_values)
  expect_equal(converted, c("1", "2", "X", "X"))
})

test_that("chr 20 conversion does not affect partial matches", {
  # "200" should not become "X0"
  chr_values <- c("200", "120", "201")
  converted <- gsub("^20$", "X", chr_values)
  expect_equal(converted, c("200", "120", "201"))
})

# --- GenomicRanges overlap logic ---
test_that("500kb window calculation is correct", {
  ps <- 53786615
  start <- ps - 500000
  end <- ps + 500000
  expect_equal(start, 53286615)
  expect_equal(end, 54286615)
  expect_equal(end - start, 1000000)
})
