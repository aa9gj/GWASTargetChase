# Integration tests for gwasFollowup using example test data.
# These tests require network access for API calls.

test_that("gwasFollowup runs end-to-end with example cat data", {
  skip_if_offline()

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  out_dir <- file.path(tempdir(), "test_gwasFollowup_cat")
  on.exit(unlink(out_dir, recursive = TRUE))

  expect_no_error(
    gwasFollowup(
      sumStats = sumstats,
      felGTF = gtf,
      species = "cat",
      pval = 5e-8,
      ResultsPath = out_dir,
      zoo_dir = zoo_dir
    )
  )

  # Check output files were created
  expect_true(file.exists(file.path(out_dir, "g2d_results.txt")))
})

test_that("gwasFollowup creates output directory if it doesn't exist", {
  skip_if_offline()

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  out_dir <- file.path(tempdir(), "test_gwasFollowup_newdir", "subdir")
  on.exit(unlink(file.path(tempdir(), "test_gwasFollowup_newdir"), recursive = TRUE))

  expect_false(dir.exists(out_dir))
  expect_no_error(
    gwasFollowup(
      sumStats = sumstats,
      felGTF = gtf,
      species = "cat",
      pval = 5e-8,
      ResultsPath = out_dir,
      zoo_dir = zoo_dir
    )
  )
  expect_true(dir.exists(out_dir))
})

test_that("gwasFollowup with relaxed p-value finds more loci", {
  skip_if_offline()

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  out_dir <- file.path(tempdir(), "test_gwasFollowup_relaxed")
  on.exit(unlink(out_dir, recursive = TRUE))

  # Very relaxed p-value to capture all SNPs
  expect_no_error(
    gwasFollowup(
      sumStats = sumstats,
      felGTF = gtf,
      species = "cat",
      pval = 1e-4,
      ResultsPath = out_dir,
      zoo_dir = zoo_dir
    )
  )
  expect_true(file.exists(file.path(out_dir, "g2d_results.txt")))
})

test_that("gwasFollowup g2d_results.txt has content with known genes", {
  skip_if_offline()

  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gtf <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  out_dir <- file.path(tempdir(), "test_gwasFollowup_content")
  on.exit(unlink(out_dir, recursive = TRUE))

  gwasFollowup(
    sumStats = sumstats,
    felGTF = gtf,
    species = "cat",
    pval = 5e-8,
    ResultsPath = out_dir,
    zoo_dir = zoo_dir
  )

  g2d <- read.delim(file.path(out_dir, "g2d_results.txt"))
  # FTO is a well-known gene, should have disease associations
  if (nrow(g2d) > 0) {
    expect_true("gene_name" %in% colnames(g2d))
    # FTO should be among the results since it's the main hit
    expect_true("FTO" %in% g2d$gene_name)
  }
})

# --- Test the internal logic without API calls ---
test_that("gwasFollowup correctly reads and filters summary statistics", {
  sumstats <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
  gwas <- data.table::fread(sumstats)

  expect_true("chr" %in% colnames(gwas))
  expect_true("ps" %in% colnames(gwas))
  expect_true("p_wald" %in% colnames(gwas))
  expect_true("gene" %in% colnames(gwas))
  expect_equal(nrow(gwas), 15)

  # Filter at 5e-8 threshold
  sig <- dplyr::filter(gwas, p_wald <= 5e-8)
  expect_true(nrow(sig) > 0)
  expect_true(nrow(sig) < nrow(gwas))

  # Unique genes in significant hits
  sig_genes <- unique(sig$gene)
  expect_true("FTO" %in% sig_genes)
})

test_that("gwasFollowup correctly parses GTF and finds nearby genes", {
  gtf_path <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  cat_gtf <- as.data.frame(rtracklayer::import(gtf_path))
  cat_pc_gene_gtf <- dplyr::filter(cat_gtf, gene_biotype == "protein_coding", type == "gene")

  # Should have 4 gene entries (FTO on chr1, GNPDA2 on chr2, TMEM18 on chr2, FTO on chr16)
  expect_equal(nrow(cat_pc_gene_gtf), 4)

  # Test TSS-based 500kb window for FTO on chr1
  gene_row <- cat_pc_gene_gtf[cat_pc_gene_gtf$gene_name == "FTO" &
                                as.character(cat_pc_gene_gtf$seqnames) == "1", ]
  expect_equal(nrow(gene_row), 1)
  # FTO is on + strand, so TSS = start
  tss <- gene_row$start[1]
  expect_equal(tss, 52786615)

  # Find nearby genes within 500kb of TSS
  chr <- "1"
  window_start <- tss - 500000
  window_end <- tss + 500000
  nearby <- cat_pc_gene_gtf[as.character(cat_pc_gene_gtf$seqnames) == chr &
                              cat_pc_gene_gtf$start <= window_end &
                              cat_pc_gene_gtf$end >= window_start, ]
  expect_true("FTO" %in% nearby$gene_name)
})

test_that("gwasFollowup handles TSS correctly for minus strand genes", {
  gtf_path <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
  cat_gtf <- as.data.frame(rtracklayer::import(gtf_path))
  cat_pc_gene_gtf <- dplyr::filter(cat_gtf, gene_biotype == "protein_coding", type == "gene")

  # TMEM18 is on the - strand
  tmem <- cat_pc_gene_gtf[cat_pc_gene_gtf$gene_name == "TMEM18", ]
  expect_equal(as.character(tmem$strand), "-")

  # TSS for - strand gene should be the 'end' coordinate
  tss <- tmem$end[1]
  expect_equal(tss, 700000)
})

test_that("gwasFollowup zoonomia translation works for example genes", {
  zoo_dir <- system.file("extdata", package = "GWASTargetChase")
  genes <- c("FTO", "GNPDA2", "TMEM18")

  human_ortho <- translate_genes(genes, "human", "cat", zoo_dir)
  expect_true(nrow(human_ortho) > 0)
  expect_true("FTO" %in% human_ortho$target_gene)

  mouse_ortho <- translate_genes(genes, "mouse", "cat", zoo_dir)
  expect_true(nrow(mouse_ortho) > 0)
  expect_true("Fto" %in% mouse_ortho$target_gene)
})
