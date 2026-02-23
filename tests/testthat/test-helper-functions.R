test_that("gene2disease returns correct results for known genes", {
  assoc_data <- data.frame(
    gene_name = c("FTO", "FTO", "MC4R", "BRCA1"),
    gene_id = c("ENSG00000140718", "ENSG00000140718", "ENSG00000166603", "ENSG00000012048"),
    disease_name = c("obesity", "BMI", "obesity", "breast cancer"),
    score = c(0.85, 0.82, 0.78, 0.95),
    stringsAsFactors = FALSE
  )

  result <- gene2disease("FTO", assoc_data)
  expect_equal(nrow(result), 2)
  expect_true(all(result$gene_name == "FTO"))

  result_mc4r <- gene2disease("MC4R", assoc_data)
  expect_equal(nrow(result_mc4r), 1)
  expect_equal(result_mc4r$gene_name, "MC4R")
})

test_that("gene2disease returns empty data.frame for unknown genes", {
  assoc_data <- data.frame(
    gene_name = c("FTO", "MC4R"),
    gene_id = c("ENSG00000140718", "ENSG00000166603"),
    disease_name = c("obesity", "obesity"),
    score = c(0.85, 0.78),
    stringsAsFactors = FALSE
  )

  result <- gene2disease("NONEXISTENT", assoc_data)
  expect_equal(nrow(result), 0)
})

test_that("gene2disease uses word boundaries correctly", {
  assoc_data <- data.frame(
    gene_name = c("FTO", "FTOS", "SFTO"),
    gene_id = c("ENSG1", "ENSG2", "ENSG3"),
    disease_name = c("obesity", "other", "other"),
    score = c(0.85, 0.5, 0.5),
    stringsAsFactors = FALSE
  )

  result <- gene2disease("FTO", assoc_data)
  expect_equal(nrow(result), 1)
  expect_equal(result$gene_name, "FTO")
})

test_that("gene2diseaseMan returns correct results", {
  assoc_data <- data.frame(
    gene_name = c("BRCA1", "BRCA2", "TP53"),
    gene_id = c("ENSG1", "ENSG2", "ENSG3"),
    disease_name = c("breast cancer", "breast cancer", "cancer"),
    score = c(0.95, 0.94, 0.88),
    stringsAsFactors = FALSE
  )

  result <- gene2diseaseMan("BRCA1", assoc_data)
  expect_equal(nrow(result), 1)
  expect_equal(result$gene_name, "BRCA1")
})

test_that("locus2gene returns correct results for known genes", {
  l2g_data <- data.frame(
    study_id = c("GCST1", "GCST1", "GCST2"),
    gene_name = c("FTO", "FTO", "MC4R"),
    y_proba_full_model = c(0.92, 0.88, 0.85),
    trait_reported = c("BMI", "BMI", "BMI"),
    stringsAsFactors = FALSE
  )

  result <- locus2gene("FTO", l2g_data)
  expect_equal(nrow(result), 2)
  expect_true(all(result$gene_name == "FTO"))
})

test_that("locus2gene returns empty data.frame for unknown genes", {
  l2g_data <- data.frame(
    study_id = c("GCST1"),
    gene_name = c("FTO"),
    y_proba_full_model = c(0.92),
    trait_reported = c("BMI"),
    stringsAsFactors = FALSE
  )

  result <- locus2gene("NONEXISTENT", l2g_data)
  expect_equal(nrow(result), 0)
})

test_that("locus2geneMan returns correct results", {
  l2g_data <- data.frame(
    study_id = c("GCST1", "GCST2"),
    gene_name = c("TCF7L2", "PCSK9"),
    y_proba_full_model = c(0.91, 0.94),
    trait_reported = c("T2D", "LDL"),
    stringsAsFactors = FALSE
  )

  result <- locus2geneMan("TCF7L2", l2g_data)
  expect_equal(nrow(result), 1)
  expect_equal(result$gene_name, "TCF7L2")
})

test_that("load_zoonomia_orthologs extracts gene symbols correctly", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- load_zoonomia_orthologs("human", "cat", zoo_dir)
  expect_true(is.data.frame(result))
  expect_true("gene_symbol" %in% colnames(result))
  expect_true("orthology_class" %in% colnames(result))
  expect_true("V1" %in% colnames(result))
  expect_true("V3" %in% colnames(result))
  expect_true("hills_grade" %in% colnames(result))
  expect_true(nrow(result) > 0)
  # Check that gene symbols are extracted (no ENST prefix)
  expect_false(any(grepl("^ENST", result$gene_symbol)))
  # Check known gene is present
  expect_true("FTO" %in% result$gene_symbol)
})

test_that("load_zoonomia_orthologs fails for missing files", {
  expect_error(load_zoonomia_orthologs("human", "elephant"),
               "Zoonomia orthology file not found")
})

test_that("translate_genes returns correct human orthologs for cat genes", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- translate_genes(c("FTO", "GNPDA2", "TMEM18"), "human", "cat", zoo_dir)
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  expect_true("original_gene" %in% colnames(result))
  expect_true("target_gene" %in% colnames(result))
  expect_true("orthology_class" %in% colnames(result))
  expect_true("hills_grade" %in% colnames(result))
  expect_true("FTO" %in% result$target_gene)
})

test_that("translate_genes returns correct mouse orthologs for cat genes", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- translate_genes(c("FTO", "GNPDA2"), "mouse", "cat", zoo_dir)
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  # Mouse gene names may differ in case (e.g., Fto vs FTO)
  expect_true("Fto" %in% result$target_gene || "FTO" %in% toupper(result$target_gene))
})

test_that("translate_genes handles case-insensitive matching", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- translate_genes(c("fto", "gnpda2"), "human", "cat", zoo_dir)
  expect_true(nrow(result) > 0)
})

test_that("translate_genes returns empty data.frame for unknown genes", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  result <- translate_genes(c("NONEXISTENT_GENE"), "human", "cat", zoo_dir)
  expect_equal(nrow(result), 0)
  expect_true("original_gene" %in% colnames(result))
  expect_true("target_gene" %in% colnames(result))
})

test_that("add_orthology_info adds metadata columns to results", {
  results_df <- data.frame(
    gene_name = c("FTO", "MC4R", "UNKNOWN"),
    score = c(0.9, 0.8, 0.7),
    stringsAsFactors = FALSE
  )
  ortho_df <- data.frame(
    target_gene = c("FTO", "MC4R"),
    orthology_class = c("one2one", "one2one"),
    V1 = c("GENE", "GENE"),
    V3 = c("I", "I"),
    hills_grade = c("A", "A"),
    stringsAsFactors = FALSE
  )
  result <- add_orthology_info(results_df, ortho_df)
  expect_true("orthology_class" %in% colnames(result))
  expect_true("hills_grade" %in% colnames(result))
  expect_equal(nrow(result), 3)
  # UNKNOWN gene should have NA for orthology
  unknown_row <- result[result$gene_name == "UNKNOWN", ]
  expect_true(is.na(unknown_row$orthology_class))
})

test_that("add_orthology_info handles NULL ortho_df gracefully", {
  results_df <- data.frame(
    gene_name = c("FTO"),
    score = c(0.9),
    stringsAsFactors = FALSE
  )
  result <- add_orthology_info(results_df, NULL)
  expect_equal(result, results_df)
})
