#' gwasFollowupMan
#'
#' main package extension function that will take GWAS sumstats, a gtf file, and a pvalue of interest to returns information on the genes of interest but also the genetic
#' association file and l2g file data prepped with l2gPrep and geneticAssocPrep functions in this package. Only use this function when you have updated genetic and l2g data
#' from opentargets
#'
#' @param sumStats GWAS summary statistics file. It assumed a ps, and chr columns.
#' @param felGTF GTF file
#' @param species Species of the GWAS data: "human", "cat", or "dog". For non-human species, Zoonomia orthology is used to translate genes. Default is "cat".
#' @param pval default is 5*10-8 and you could change but make sure to use standard form
#' @param ResultsPath default is the working directory but you can provide a path of your own ensuring it ends with a /
#' @param impc impc data file path
#' @param assocOT results from geneticAssocPrep function (file path)
#' @param l2gOT results from l2gPrep function (file path)
#' @param zoo_dir Directory containing Zoonomia orthology files. If NULL, uses package extdata.
#' @return Writes g2d_results.txt and l2g_results.txt to ResultsPath
#' @importFrom dplyr filter left_join
#' @importFrom data.table fread
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @importFrom S4Vectors queryHits subjectHits
#' @export gwasFollowupMan

gwasFollowupMan <- function(sumStats, felGTF, species = "cat", pval = 0.00000005, ResultsPath = ".", impc = NULL, assocOT, l2gOT, zoo_dir = NULL) {
  # Create output directory if it doesn't exist
  if (!dir.exists(ResultsPath)) {
    dir.create(ResultsPath, recursive = TRUE)
  }
  print("Reading OpenTargets Genetic Association file")
  association_data <- fread(assocOT)
  print("Reading OpenTargets locus to gene file")
  l2g_data <- fread(l2gOT, sep = "\t", quote = "")
  print("Reading IMPC phenotype data")
  impc_data <- fread(impc)
  print("Reading the gtf file and filtering for protein coding genes, this might take a while")
  cat_gtf <- as.data.frame(rtracklayer::import(felGTF))
  cat_pc_gene_gtf <- filter(cat_gtf, gene_biotype == "protein_coding", type == "gene")
  print("reading the summary statistics file, filtering based on the input p_value")
  gwas <- fread(sumStats)
  gwas$chr <- gsub("^20$", "X", gwas$chr)
  sig_gwas <- filter(gwas, p_wald <= pval)
  print("Overlap genes within significant loci (500Kb upstream and downstream)")
  sig_gwas$start <- sig_gwas$ps - 500000
  sig_gwas$end <- sig_gwas$ps + 500000
  tryCatch(sig_gwas_ranges <- makeGRangesFromDataFrame(sig_gwas, keep.extra.columns = T, seqnames.field = "chr"),
           error = function(e)
             stop("nothing is being reported using this p-value, consider increasing it"))
  cat_pc_gene_ranges <- makeGRangesFromDataFrame(cat_pc_gene_gtf, keep.extra.columns = T)
  hits <- findOverlaps(cat_pc_gene_ranges, sig_gwas_ranges)
  olap <- pintersect(cat_pc_gene_ranges[queryHits(hits)],sig_gwas_ranges[subjectHits(hits)])
  overlap_genes <- as.vector(na.exclude(olap$gene_name))
  genes <- overlap_genes

  # Zoonomia orthology translation for non-human species
  human_ortho <- NULL
  mouse_ortho <- NULL
  if (species != "human") {
    print("Translating genes to human orthologs using Zoonomia data for OpenTargets")
    human_ortho <- translate_genes(genes, "human", species, zoo_dir)
    if (nrow(human_ortho) > 0) {
      ot_genes <- unique(human_ortho$target_gene)
    } else {
      warning("No human orthologs found in Zoonomia data. Using original gene names.")
      ot_genes <- genes
    }
    print("Translating genes to mouse orthologs using Zoonomia data for IMPC")
    mouse_ortho <- translate_genes(genes, "mouse", species, zoo_dir)
    if (nrow(mouse_ortho) > 0) {
      impc_genes <- unique(mouse_ortho$original_gene)
    } else {
      warning("No mouse orthologs found in Zoonomia data.")
      impc_genes <- character(0)
    }
  } else {
    ot_genes <- genes
    impc_genes <- genes
  }

  print("Genes to Disease GWAS follow-up")
  results_g2d <- list()
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(ot_genes), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for (i in seq_along(ot_genes)) {
    results_g2d[[i]] <- gene2diseaseMan(ot_genes[i], association_data)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_g2d_df <- as.data.frame(do.call("rbind", results_g2d))
  if (nrow(results_g2d_df) > 0) {
    # Only join IMPC for genes with confirmed mouse orthologs
    if (species != "human" && length(impc_genes) > 0) {
      impc_filtered <- impc_data[toupper(impc_data$gene_name) %in% toupper(impc_genes), , drop = FALSE]
      results_g2d_df <- left_join(results_g2d_df, impc_filtered, by = "gene_name")
    } else {
      results_g2d_df <- left_join(results_g2d_df, impc_data, by = "gene_name")
    }
    # Add Zoonomia orthology metadata
    results_g2d_df <- add_orthology_info(results_g2d_df, human_ortho)
  }
  print("writing results")
  write.table(results_g2d_df, file.path(ResultsPath, "g2d_results.txt"), quote = F, row.names = F, sep = "\t")
  print("Locus to Genes GWAS follow-up")
  results_l2g <- list()
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(ot_genes), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for (i in seq_along(ot_genes)) {
    results_l2g[[i]] <- locus2geneMan(ot_genes[i], l2g_data)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_l2g_df <- as.data.frame(do.call("rbind", results_l2g))
  if (nrow(results_l2g_df) > 0) {
    # Only join IMPC for genes with confirmed mouse orthologs
    if (species != "human" && length(impc_genes) > 0) {
      impc_filtered <- impc_data[toupper(impc_data$gene_name) %in% toupper(impc_genes), , drop = FALSE]
      results_l2g_df <- left_join(results_l2g_df, impc_filtered, by = "gene_name")
    } else {
      results_l2g_df <- left_join(results_l2g_df, impc_data, by = "gene_name")
    }
    # Add Zoonomia orthology metadata
    results_l2g_df <- add_orthology_info(results_l2g_df, human_ortho)
  }
  print("writing results")
  write.table(results_l2g_df, file.path(ResultsPath, "l2g_results.txt"), quote = F, row.names = F, sep = "\t")
  print("All done, happy exploring")
}
