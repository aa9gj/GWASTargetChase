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
