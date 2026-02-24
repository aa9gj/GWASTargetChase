# Extended tests for Zoonomia orthology functions

test_that("load_zoonomia_orthologs loads all four species combinations", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")

  combos <- list(
    c("human", "cat"), c("human", "dog"),
    c("mouse", "cat"), c("mouse", "dog")
  )
  for (combo in combos) {
    result <- load_zoonomia_orthologs(combo[1], combo[2], zoo_dir)
    expect_true(is.data.frame(result),
                info = paste("Failed for", combo[1], "->", combo[2]))
    expect_true(nrow(result) > 0,
                info = paste("Empty for", combo[1], "->", combo[2]))
    expect_true("gene_symbol" %in% colnames(result))
    expect_true("orthology_class" %in% colnames(result))
    expect_true("orthology_grade" %in% colnames(result))
  }
})

test_that("load_zoonomia_orthologs extracts gene symbols without ENST prefix", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- load_zoonomia_orthologs("human", "cat", zoo_dir)
  # No gene_symbol should start with ENST or ENSG

  expect_false(any(grepl("^ENST", result$gene_symbol)))
  expect_false(any(grepl("^ENSG", result$gene_symbol)))
})

test_that("load_zoonomia_orthologs returns expected genes for human-cat", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- load_zoonomia_orthologs("human", "cat", zoo_dir)
  expected_genes <- c("FTO", "GNPDA2", "TMEM18", "MC4R", "LEP", "BRCA1", "BRCA2", "TP53")
  for (g in expected_genes) {
    expect_true(g %in% result$gene_symbol,
                info = paste("Expected gene", g, "not found in human-cat orthologs"))
  }
})

test_that("load_zoonomia_orthologs returns mouse gene names for mouse-cat", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- load_zoonomia_orthologs("mouse", "cat", zoo_dir)
  # Mouse gene names use lowercase convention (e.g., Fto, not FTO)
  expect_true("Fto" %in% result$gene_symbol)
  expect_true("Gnpda2" %in% result$gene_symbol)
  expect_true("Tmem18" %in% result$gene_symbol)
})

test_that("load_zoonomia_orthologs preserves orthology_class values", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- load_zoonomia_orthologs("human", "cat", zoo_dir)
  expect_true("one2one" %in% result$orthology_class)
  expect_true("many2many" %in% result$orthology_class)
})

test_that("load_zoonomia_orthologs preserves orthology_grade values", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- load_zoonomia_orthologs("human", "cat", zoo_dir)
  valid_grades <- c("A", "B", "C")
  expect_true(all(result$orthology_grade %in% valid_grades))
})

test_that("translate_genes returns all expected columns", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- translate_genes(c("FTO"), "human", "cat", zoo_dir)
  expected_cols <- c("original_gene", "target_gene", "orthology_class",
                     "V1", "V3", "orthology_grade")
  for (col in expected_cols) {
    expect_true(col %in% colnames(result),
                info = paste("Missing column:", col))
  }
})

test_that("translate_genes maps cat genes to correct human orthologs", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- translate_genes(c("FTO", "GNPDA2", "TMEM18"), "human", "cat", zoo_dir)
  expect_true("FTO" %in% result$target_gene)
  expect_true("GNPDA2" %in% result$target_gene)
  expect_true("TMEM18" %in% result$target_gene)
  expect_equal(nrow(result), 3)
})

test_that("translate_genes maps cat genes to correct mouse orthologs", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- translate_genes(c("FTO", "GNPDA2", "TMEM18"), "mouse", "cat", zoo_dir)
  expect_true("Fto" %in% result$target_gene)
  expect_true("Gnpda2" %in% result$target_gene)
  expect_true("Tmem18" %in% result$target_gene)
})

test_that("translate_genes handles mixed case input", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result_lower <- translate_genes(c("fto", "gnpda2"), "human", "cat", zoo_dir)
  result_upper <- translate_genes(c("FTO", "GNPDA2"), "human", "cat", zoo_dir)
  expect_equal(nrow(result_lower), nrow(result_upper))
})

test_that("translate_genes returns empty df for genes not in zoonomia", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- translate_genes(c("FAKEGENE1", "NOTREAL"), "human", "cat", zoo_dir)
  expect_equal(nrow(result), 0)
  expect_true("original_gene" %in% colnames(result))
  expect_true("target_gene" %in% colnames(result))
})

test_that("translate_genes handles mix of known and unknown genes", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- translate_genes(c("FTO", "NOTREAL", "GNPDA2"), "human", "cat", zoo_dir)
  expect_equal(nrow(result), 2)
  expect_true("FTO" %in% result$target_gene)
  expect_true("GNPDA2" %in% result$target_gene)
})

test_that("translate_genes deduplicates results", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  # FTO appears twice in the zoonomia file with different transcripts
  result <- translate_genes(c("FTO"), "human", "cat", zoo_dir)
  # Should deduplicate to 1 unique result
  expect_equal(nrow(result), 1)
})

test_that("add_orthology_info merges correctly with matching genes", {
  results_df <- data.frame(
    gene_name = c("FTO", "GNPDA2", "TMEM18"),
    score = c(0.9, 0.8, 0.7),
    stringsAsFactors = FALSE
  )
  ortho_df <- data.frame(
    target_gene = c("FTO", "GNPDA2", "TMEM18"),
    orthology_class = c("one2one", "one2one", "one2one"),
    V1 = c("GENE", "GENE", "GENE"),
    V3 = c("I", "I", "I"),
    orthology_grade = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  result <- add_orthology_info(results_df, ortho_df)
  expect_equal(nrow(result), 3)
  expect_true("orthology_class" %in% colnames(result))
  expect_true("orthology_grade" %in% colnames(result))
  # FTO should have grade A
  fto_row <- result[result$gene_name == "FTO", ]
  expect_equal(fto_row$orthology_grade, "A")
  # TMEM18 should have grade B
  tmem_row <- result[result$gene_name == "TMEM18", ]
  expect_equal(tmem_row$orthology_grade, "B")
})

test_that("add_orthology_info sets NA for unmatched genes", {
  results_df <- data.frame(
    gene_name = c("FTO", "UNKNOWN_GENE"),
    score = c(0.9, 0.5),
    stringsAsFactors = FALSE
  )
  ortho_df <- data.frame(
    target_gene = c("FTO"),
    orthology_class = c("one2one"),
    V1 = c("GENE"),
    V3 = c("I"),
    orthology_grade = c("A"),
    stringsAsFactors = FALSE
  )
  result <- add_orthology_info(results_df, ortho_df)
  expect_equal(nrow(result), 2)
  unknown <- result[result$gene_name == "UNKNOWN_GENE", ]
  expect_true(is.na(unknown$orthology_class))
  expect_true(is.na(unknown$orthology_grade))
})

test_that("add_orthology_info returns unchanged df when ortho_df is empty", {
  results_df <- data.frame(
    gene_name = c("FTO"),
    score = c(0.9),
    stringsAsFactors = FALSE
  )
  empty_ortho <- data.frame(
    target_gene = character(),
    orthology_class = character(),
    V1 = character(),
    V3 = character(),
    orthology_grade = character(),
    stringsAsFactors = FALSE
  )
  result <- add_orthology_info(results_df, empty_ortho)
  expect_equal(result, results_df)
})
