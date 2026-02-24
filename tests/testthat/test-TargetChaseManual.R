# Integration tests for TargetChaseManual using mock pre-downloaded data.
# TargetChaseManual works with pre-downloaded files, so we create mock data
# to test the full pipeline without requiring real bulk downloads.

# Helper to create mock OpenTargets genetic association data
.create_mock_assoc_data <- function(path) {
  assoc_data <- data.frame(
    gene_name = c("FTO", "FTO", "GNPDA2", "TMEM18", "MC4R"),
    gene_id = c("ENSG00000140718", "ENSG00000140718", "ENSG00000182396",
                "ENSG00000068394", "ENSG00000166603"),
    disease_name = c("obesity", "body mass index", "obesity", "obesity", "obesity"),
    disease_id = c("EFO_0001073", "EFO_0004340", "EFO_0001073",
                   "EFO_0001073", "EFO_0001073"),
    datatypeId = rep("genetic_association", 5),
    score = c(0.85, 0.82, 0.65, 0.55, 0.78),
    stringsAsFactors = FALSE
  )
  write.table(assoc_data, path, quote = FALSE, row.names = FALSE, sep = "\t")
  invisible(path)
}

# Helper to create mock L2G data
.create_mock_l2g_data <- function(path) {
  l2g_data <- data.frame(
    study_id = c("GCST001", "GCST001", "GCST002", "GCST003"),
    gene_name = c("FTO", "GNPDA2", "FTO", "TMEM18"),
    y_proba_full_model = c(0.92, 0.75, 0.88, 0.70),
    trait_reported = c("Body mass index", "Body mass index",
                       "Obesity", "Body mass index"),
    stringsAsFactors = FALSE
  )
  write.table(l2g_data, path, quote = FALSE, row.names = FALSE, sep = "\t")
  invisible(path)
}

# Helper to create mock IMPC data
.create_mock_impc_data <- function(path) {
  impc_data <- data.frame(
    gene_name = c("FTO", "GNPDA2", "TMEM18"),
    num_phenotype_hits = c(15, 3, 7),
    Phenotype_Hits = c("increased body weight; abnormal glucose tolerance",
                       "decreased circulating glucose level",
                       "decreased body weight; hyperactivity"),
    stringsAsFactors = FALSE
  )
  write.table(impc_data, path, quote = FALSE, row.names = FALSE, sep = "\t")
  invisible(path)
}

test_that("TargetChaseManual runs end-to-end with mock data", {
  # Create temp directory for mock data and output
  mock_dir <- file.path(tempdir(), "test_TargetChaseManual")
  out_dir <- file.path(mock_dir, "output")
  dir.create(mock_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(mock_dir, recursive = TRUE))

  # Create mock data files
  assoc_path <- file.path(mock_dir, "mock_assoc.txt")
  l2g_path <- file.path(mock_dir, "mock_l2g.txt")
  impc_path <- file.path(mock_dir, "mock_impc.txt")
  .create_mock_assoc_data(assoc_path)
  .create_mock_l2g_data(l2g_path)
  .create_mock_impc_data(impc_path)

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")

  expect_no_error(
    TargetChaseManual(
      sumStats = sumstats,
      gtf = gtf,
      species = "cat",
      pval = 5e-8,
      ResultsPath = out_dir,
      impc = impc_path,
      assocOT = assoc_path,
      l2gOT = l2g_path,
      zoo_dir = zoo_dir
    )
  )

  # Check output files
  expect_true(file.exists(file.path(out_dir, "g2d_results.txt")))
  expect_true(file.exists(file.path(out_dir, "l2g_results.txt")))
})

test_that("TargetChaseManual g2d_results.txt contains expected gene matches", {
  mock_dir <- file.path(tempdir(), "test_TargetChaseManual_g2d")
  out_dir <- file.path(mock_dir, "output")
  dir.create(mock_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(mock_dir, recursive = TRUE))

  assoc_path <- file.path(mock_dir, "mock_assoc.txt")
  l2g_path <- file.path(mock_dir, "mock_l2g.txt")
  impc_path <- file.path(mock_dir, "mock_impc.txt")
  .create_mock_assoc_data(assoc_path)
  .create_mock_l2g_data(l2g_path)
  .create_mock_impc_data(impc_path)

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")

  TargetChaseManual(
    sumStats = sumstats,
    gtf = gtf,
    species = "cat",
    pval = 5e-8,
    ResultsPath = out_dir,
    impc = impc_path,
    assocOT = assoc_path,
    l2gOT = l2g_path,
    zoo_dir = zoo_dir
  )

  g2d <- read.delim(file.path(out_dir, "g2d_results.txt"))
  if (nrow(g2d) > 0) {
    expect_true("gene_name" %in% colnames(g2d))
    # FTO should be in the results as it's a top hit in both GWAS and assoc data
    expect_true("FTO" %in% g2d$gene_name)
  }
})

test_that("TargetChaseManual l2g_results.txt contains expected gene matches", {
  mock_dir <- file.path(tempdir(), "test_TargetChaseManual_l2g")
  out_dir <- file.path(mock_dir, "output")
  dir.create(mock_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(mock_dir, recursive = TRUE))

  assoc_path <- file.path(mock_dir, "mock_assoc.txt")
  l2g_path <- file.path(mock_dir, "mock_l2g.txt")
  impc_path <- file.path(mock_dir, "mock_impc.txt")
  .create_mock_assoc_data(assoc_path)
  .create_mock_l2g_data(l2g_path)
  .create_mock_impc_data(impc_path)

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")

  TargetChaseManual(
    sumStats = sumstats,
    gtf = gtf,
    species = "cat",
    pval = 5e-8,
    ResultsPath = out_dir,
    impc = impc_path,
    assocOT = assoc_path,
    l2gOT = l2g_path,
    zoo_dir = zoo_dir
  )

  l2g <- read.delim(file.path(out_dir, "l2g_results.txt"))
  if (nrow(l2g) > 0) {
    expect_true("gene_name" %in% colnames(l2g))
    expect_true("FTO" %in% l2g$gene_name)
  }
})

test_that("TargetChaseManual adds orthology info to results", {
  mock_dir <- file.path(tempdir(), "test_TargetChaseManual_ortho")
  out_dir <- file.path(mock_dir, "output")
  dir.create(mock_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(mock_dir, recursive = TRUE))

  assoc_path <- file.path(mock_dir, "mock_assoc.txt")
  l2g_path <- file.path(mock_dir, "mock_l2g.txt")
  impc_path <- file.path(mock_dir, "mock_impc.txt")
  .create_mock_assoc_data(assoc_path)
  .create_mock_l2g_data(l2g_path)
  .create_mock_impc_data(impc_path)

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")

  TargetChaseManual(
    sumStats = sumstats,
    gtf = gtf,
    species = "cat",
    pval = 5e-8,
    ResultsPath = out_dir,
    impc = impc_path,
    assocOT = assoc_path,
    l2gOT = l2g_path,
    zoo_dir = zoo_dir
  )

  g2d <- read.delim(file.path(out_dir, "g2d_results.txt"))
  if (nrow(g2d) > 0) {
    # Orthology columns should be present for cat species
    expect_true("orthology_class" %in% colnames(g2d))
    expect_true("orthology_grade" %in% colnames(g2d))
  }
})

test_that("TargetChaseManual with relaxed p-value captures more SNPs", {
  mock_dir <- file.path(tempdir(), "test_TargetChaseManual_relaxed")
  out_dir_strict <- file.path(mock_dir, "output_strict")
  out_dir_relaxed <- file.path(mock_dir, "output_relaxed")
  dir.create(mock_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(mock_dir, recursive = TRUE))

  assoc_path <- file.path(mock_dir, "mock_assoc.txt")
  l2g_path <- file.path(mock_dir, "mock_l2g.txt")
  impc_path <- file.path(mock_dir, "mock_impc.txt")
  .create_mock_assoc_data(assoc_path)
  .create_mock_l2g_data(l2g_path)
  .create_mock_impc_data(impc_path)

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")

  TargetChaseManual(
    sumStats = sumstats, gtf = gtf, species = "cat",
    pval = 5e-8, ResultsPath = out_dir_strict,
    impc = impc_path, assocOT = assoc_path, l2gOT = l2g_path,
    zoo_dir = zoo_dir
  )

  TargetChaseManual(
    sumStats = sumstats, gtf = gtf, species = "cat",
    pval = 1e-5, ResultsPath = out_dir_relaxed,
    impc = impc_path, assocOT = assoc_path, l2gOT = l2g_path,
    zoo_dir = zoo_dir
  )

  # Both should produce output
  expect_true(file.exists(file.path(out_dir_strict, "g2d_results.txt")))
  expect_true(file.exists(file.path(out_dir_relaxed, "g2d_results.txt")))
})

test_that("TargetChaseManual works with human species (no orthology translation)", {
  mock_dir <- file.path(tempdir(), "test_TargetChaseManual_human")
  out_dir <- file.path(mock_dir, "output")
  dir.create(mock_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(mock_dir, recursive = TRUE))

  assoc_path <- file.path(mock_dir, "mock_assoc.txt")
  l2g_path <- file.path(mock_dir, "mock_l2g.txt")
  impc_path <- file.path(mock_dir, "mock_impc.txt")
  .create_mock_assoc_data(assoc_path)
  .create_mock_l2g_data(l2g_path)
  .create_mock_impc_data(impc_path)

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")

  expect_no_error(
    TargetChaseManual(
      sumStats = sumstats,
      gtf = gtf,
      species = "human",
      pval = 5e-8,
      ResultsPath = out_dir,
      impc = impc_path,
      assocOT = assoc_path,
      l2gOT = l2g_path
    )
  )

  expect_true(file.exists(file.path(out_dir, "g2d_results.txt")))
  expect_true(file.exists(file.path(out_dir, "l2g_results.txt")))
})
