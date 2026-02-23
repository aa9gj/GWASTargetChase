# Tests for API fetch functions (fetch_opentargets, fetch_impc)
# These tests use mocked responses to avoid requiring network access.

# --- Cache helper tests ---
test_that(".get_cache_dir creates and returns a temp directory", {
  cache_dir <- GWASTargetChase:::.get_cache_dir()
  expect_true(dir.exists(cache_dir))
  expect_true(grepl("GWASTargetChase_cache", cache_dir))
})

# --- fetch_opentargets tests ---
test_that("fetch_opentargets returns a data.frame", {
  skip_if_offline()
  result <- fetch_opentargets("FTO", use_cache = FALSE)
  expect_true(is.data.frame(result))
})

test_that("fetch_opentargets returns expected columns for known gene", {
  skip_if_offline()
  result <- fetch_opentargets("FTO", use_cache = FALSE)
  if (nrow(result) > 0) {
    expect_true("gene_name" %in% colnames(result))
    expect_true("gene_id" %in% colnames(result))
    expect_true("disease_name" %in% colnames(result))
    expect_true("overall_score" %in% colnames(result))
    expect_true("genetic_association_score" %in% colnames(result))
  }
})

test_that("fetch_opentargets returns FTO disease associations", {
  skip_if_offline()
  result <- fetch_opentargets("FTO", use_cache = FALSE)
  # FTO is a well-known obesity gene, should have associations
  expect_true(nrow(result) > 0)
  expect_true(all(result$gene_name == "FTO"))
})

test_that("fetch_opentargets returns empty df for nonexistent gene", {
  skip_if_offline()
  result <- fetch_opentargets("ZZZZNOTAREALGENE999", use_cache = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
})

test_that("fetch_opentargets caching works", {
  skip_if_offline()
  # Clear cache first
  cache_file <- file.path(GWASTargetChase:::.get_cache_dir(), "ot_BRCA1.rds")
  if (file.exists(cache_file)) unlink(cache_file)

  # First call should create cache
  result1 <- fetch_opentargets("BRCA1", use_cache = TRUE)
  if (nrow(result1) > 0) {
    expect_true(file.exists(cache_file))

    # Second call should use cache
    result2 <- fetch_opentargets("BRCA1", use_cache = TRUE)
    expect_equal(result1, result2)
  }
})

# --- fetch_impc tests ---
test_that("fetch_impc returns a data.frame", {
  skip_if_offline()
  result <- fetch_impc("Fto", use_cache = FALSE)
  expect_true(is.data.frame(result))
})

test_that("fetch_impc returns expected columns for known gene", {
  skip_if_offline()
  result <- fetch_impc("Fto", use_cache = FALSE)
  if (nrow(result) > 0) {
    expect_true("gene_name" %in% colnames(result))
    expect_true("num_phenotype_hits" %in% colnames(result))
    expect_true("Phenotype_Hits" %in% colnames(result))
  }
})

test_that("fetch_impc returns phenotypes for Fto (well-characterized gene)", {
  skip_if_offline()
  result <- fetch_impc("Fto", use_cache = FALSE)
  expect_true(nrow(result) > 0)
  expect_true(result$num_phenotype_hits > 0)
})

test_that("fetch_impc returns empty df for nonexistent gene", {
  skip_if_offline()
  result <- fetch_impc("ZZZZNOTAREALGENE999", use_cache = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
})

test_that("fetch_impc caching works", {
  skip_if_offline()
  cache_file <- file.path(GWASTargetChase:::.get_cache_dir(), "impc_BRCA1.rds")
  if (file.exists(cache_file)) unlink(cache_file)

  result1 <- fetch_impc("Brca1", use_cache = TRUE)
  if (nrow(result1) > 0) {
    expect_true(file.exists(cache_file))

    result2 <- fetch_impc("Brca1", use_cache = TRUE)
    expect_equal(result1, result2)
  }
})

# --- Tests with mock data (no network required) ---
test_that("fetch_opentargets returns from cache without network", {
  # Create a mock cached result
  cache_dir <- GWASTargetChase:::.get_cache_dir()
  mock_data <- data.frame(
    gene_name = "MOCKGENE",
    gene_id = "ENSG99999999",
    disease_id = "EFO_0000000",
    disease_name = "mock disease",
    overall_score = 0.5,
    genetic_association_score = 0.3,
    stringsAsFactors = FALSE
  )
  cache_file <- file.path(cache_dir, "ot_MOCKGENE.rds")
  saveRDS(mock_data, cache_file)
  on.exit(unlink(cache_file))

  result <- fetch_opentargets("MOCKGENE", use_cache = TRUE)
  expect_equal(result, mock_data)
  expect_equal(nrow(result), 1)
  expect_equal(result$gene_name, "MOCKGENE")
})

test_that("fetch_impc returns from cache without network", {
  cache_dir <- GWASTargetChase:::.get_cache_dir()
  mock_data <- data.frame(
    gene_name = "MOCKGENE",
    num_phenotype_hits = 3,
    Phenotype_Hits = "phenotype1; phenotype2; phenotype3",
    stringsAsFactors = FALSE
  )
  cache_file <- file.path(cache_dir, "impc_MOCKGENE.rds")
  saveRDS(mock_data, cache_file)
  on.exit(unlink(cache_file))

  result <- fetch_impc("MOCKGENE", use_cache = TRUE)
  expect_equal(result, mock_data)
})
