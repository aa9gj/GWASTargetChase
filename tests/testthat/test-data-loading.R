test_that("package data files exist", {
  expect_true(file.exists(system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")))
  expect_true(file.exists(system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")))
})

test_that("example GWAS summary statistics has expected structure", {
  sumstats_path <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gwas <- read.delim(sumstats_path)
  expect_true("chr" %in% colnames(gwas))
  expect_true("rs" %in% colnames(gwas))
  expect_true("ps" %in% colnames(gwas))
  expect_true("p_wald" %in% colnames(gwas))
  expect_true("gene" %in% colnames(gwas))
  expect_true(nrow(gwas) > 0)
})

test_that("zoonomia orthology files exist", {
  expect_true(file.exists(system.file("extdata", "zoo_human_cat.txt", package = "GWASTargetChase")))
  expect_true(file.exists(system.file("extdata", "zoo_human_dog.txt", package = "GWASTargetChase")))
  expect_true(file.exists(system.file("extdata", "zoo_mouse_cat.txt", package = "GWASTargetChase")))
  expect_true(file.exists(system.file("extdata", "zoo_mouse_dog.txt", package = "GWASTargetChase")))
})

test_that("zoonomia files have expected structure", {
  zoo_path <- system.file("extdata", "zoo_human_cat.txt", package = "GWASTargetChase")
  zoo <- read.delim(zoo_path, sep = "\t")
  expect_true("t_gene" %in% colnames(zoo))
  expect_true("t_transcript" %in% colnames(zoo))
  expect_true("q_gene" %in% colnames(zoo))
  expect_true("q_transcript" %in% colnames(zoo))
  expect_true("orthology_class" %in% colnames(zoo))
  expect_true("V1" %in% colnames(zoo))
  expect_true("V3" %in% colnames(zoo))
  expect_true("orthology_grade" %in% colnames(zoo))
  expect_true(nrow(zoo) > 0)
})
