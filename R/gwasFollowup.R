#' gwasFollowup
#'
#' main package function that will take GWAS sumstats, a gtf file, and a pvalue of interest to returns information on the genes of interest
#'
#' @param sumStats GWAS summary statistics file. It assumed a ps, and chr columns.
#' @param felGTF GTF file
#' @param species Species of the GWAS data: "human", "cat", or "dog". For non-human species, Zoonomia orthology is used to translate genes. Default is "cat".
#' @param pval default is 5*10-8 and you could change but make sure to use standard form
#' @param ResultsPath default is the working directory but you can provide a path of your own ensuring it ends with a /
#' @param zoo_dir Directory containing Zoonomia orthology files. If NULL, uses package extdata.
#' @return Writes g2d_results.txt, l2g_results.txt, l2g_filtered_results.txt, and plots.pdf to ResultsPath
#' @importFrom dplyr filter left_join inner_join %>%
#' @importFrom data.table fread
#' @import GenomicRanges
#' @import ggplot2
#' @import patchwork
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices pdf dev.off
#' @importFrom rtracklayer import
#' @importFrom S4Vectors queryHits subjectHits
#' @export gwasFollowup

gwasFollowup <- function(sumStats, felGTF, species = "cat", pval = 0.00000005, ResultsPath = ".", zoo_dir = NULL) {
  # Validate input files
  if (!nzchar(sumStats) || !file.exists(sumStats)) {
    stop("Summary statistics file not found: '", sumStats, "'. ",
         "Check the file path or verify the package is installed correctly.")
  }
  if (!nzchar(felGTF) || !file.exists(felGTF)) {
    stop("GTF file not found: '", felGTF, "'. ",
         "Check the file path or verify the package is installed correctly.")
  }
  # Create output directory if it doesn't exist
  if (!dir.exists(ResultsPath)) {
    dir.create(ResultsPath, recursive = TRUE)
  }
  # Load bundled package data
  print("Loading OpenTargets Genetic Association file from this package")
  suppressWarnings(data("disease_target_genetic_association", envir = environment()))
  if (!exists("disease_target_genetic_association", envir = environment())) {
    stop("Required data 'disease_target_genetic_association' not found. ",
         "The package may not be installed correctly. Try reinstalling:\n",
         "  remove.packages('GWASTargetChase')\n",
         "  devtools::install_github('aa9gj/GWASTargetChase')")
  }
  print("Loading OpenTargets locus to gene file from this package")
  suppressWarnings(data("l2g_annotated_full", envir = environment()))
  if (!exists("l2g_annotated_full", envir = environment())) {
    stop("Required data 'l2g_annotated_full' not found. ",
         "The package may not be installed correctly. Try reinstalling:\n",
         "  remove.packages('GWASTargetChase')\n",
         "  devtools::install_github('aa9gj/GWASTargetChase')")
  }
  print("loading IMPC phenotype data from this package")
  suppressWarnings(data("impc", envir = environment()))
  if (!exists("impc", envir = environment())) {
    stop("Required data 'impc' not found. ",
         "The package may not be installed correctly. Try reinstalling:\n",
         "  remove.packages('GWASTargetChase')\n",
         "  devtools::install_github('aa9gj/GWASTargetChase')")
  }
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
  sig_gwas$locus <- paste0("locus_", seq_len(nrow(sig_gwas)))
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

  print("Getting SNP information")
  hits2 <- findOverlaps(sig_gwas_ranges, cat_pc_gene_ranges)
  olap2 <- pintersect(sig_gwas_ranges[queryHits(hits2)],cat_pc_gene_ranges[subjectHits(hits2)])
  olap2_df <- as.data.frame(olap2)
  olap2_df$join <- paste0(olap2_df$seqnames, "_", olap2_df$start, "_", olap2_df$end)
  hits3 <- findOverlaps(cat_pc_gene_ranges, sig_gwas_ranges)
  olap3 <- pintersect(cat_pc_gene_ranges[queryHits(hits3)],sig_gwas_ranges[subjectHits(hits3)])
  olap3_df <- as.data.frame(olap3)
  # Keep essential columns: seqnames, start, end, width, strand, and gene_name
  keep_cols <- c("seqnames", "start", "end", "width", "strand", "gene_name")
  keep_cols <- keep_cols[keep_cols %in% colnames(olap3_df)]
  olap3_df <- olap3_df[, keep_cols, drop = FALSE]
  olap3_df$join <- paste0(olap3_df$seqnames, "_", olap3_df$start, "_", olap3_df$end)
  sig_SNP_info <- inner_join(olap2_df, olap3_df, by = "join")
  # Keep relevant columns by name instead of hardcoded indices
  snp_keep <- intersect(c("rs", "ps", "p_wald", "af", "beta", "se", "locus",
                           "gene_name", "seqnames.y", "start.y", "end.y"),
                        colnames(sig_SNP_info))
  sig_SNP_info <- sig_SNP_info[, snp_keep, drop = FALSE]
  print("Genes to Disease GWAS follow-up")
  results_g2d <- list()
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(ot_genes), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for (i in seq_along(ot_genes)) {
    results_g2d[[i]] <- gene2disease(ot_genes[i], disease_target_genetic_association)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_g2d_df <- as.data.frame(do.call("rbind", results_g2d))
  if (nrow(results_g2d_df) > 0) {
    # Only join IMPC for genes with confirmed mouse orthologs
    if (species != "human" && length(impc_genes) > 0) {
      impc_filtered <- impc[toupper(impc$gene_name) %in% toupper(impc_genes), , drop = FALSE]
      results_g2d_df <- left_join(results_g2d_df, impc_filtered, by = "gene_name")
    } else {
      results_g2d_df <- left_join(results_g2d_df, impc, by = "gene_name")
    }
    if ("gene_name" %in% colnames(sig_SNP_info)) {
      results_g2d_df <- inner_join(results_g2d_df, sig_SNP_info, by = "gene_name")
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
    results_l2g[[i]] <- locus2gene(ot_genes[i], l2g_annotated_full)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_l2g_df <- as.data.frame(do.call("rbind", results_l2g))
  if (nrow(results_l2g_df) > 0) {
    # Only join IMPC for genes with confirmed mouse orthologs
    if (species != "human" && length(impc_genes) > 0) {
      impc_filtered <- impc[toupper(impc$gene_name) %in% toupper(impc_genes), , drop = FALSE]
      results_l2g_df <- left_join(results_l2g_df, impc_filtered, by = "gene_name")
    } else {
      results_l2g_df <- left_join(results_l2g_df, impc, by = "gene_name")
    }
    if ("gene_name" %in% colnames(sig_SNP_info)) {
      results_l2g_df <- inner_join(results_l2g_df, sig_SNP_info, by = "gene_name")
    }
    # Add Zoonomia orthology metadata
    results_l2g_df <- add_orthology_info(results_l2g_df, human_ortho)
    # Create filtered L2G results with key columns
    available_cols <- colnames(results_l2g_df)
    filter_cols <- intersect(c("gene_name", "study_id", "variant_id", "gene_id",
                                "y_proba_full_model", "y_proba_logi_distance",
                                "y_proba_logi_interaction", "y_proba_logi_molecularQTL",
                                "y_proba_logi_pathogenicity", "trait_reported",
                                "trait_category", "Phenotype_Hits",
                                "rs", "ps", "p_wald", "locus",
                                "orthology_class", "V1", "V3", "orthology_grade"), available_cols)
    results_l2g_filter <- results_l2g_df[, filter_cols, drop = FALSE]
    results_l2g_filter$gene_cards <- paste0("https://www.genecards.org/Search/Keyword?queryString=", results_l2g_filter$gene_name)
    if ("y_proba_full_model" %in% colnames(results_l2g_filter)) {
      names(results_l2g_filter)[names(results_l2g_filter) == 'y_proba_full_model'] <- 'full_l2g_score'
    }
    if ("Phenotype_Hits" %in% colnames(results_l2g_filter)) {
      names(results_l2g_filter)[names(results_l2g_filter) == 'Phenotype_Hits'] <- 'IMPC_results'
    }
  } else {
    results_l2g_filter <- results_l2g_df
  }
  print("writing results")
  write.table(results_l2g_df, file.path(ResultsPath, "l2g_results.txt"), quote = F, row.names = F, sep = "\t")
  write.table(results_l2g_filter, file.path(ResultsPath, "l2g_filtered_results.txt"), quote = F, row.names = F, sep = "\t")
  print("Creating plots per significant locus, this might take a while")
  gwas_sum_ranges <- makeGRangesFromDataFrame(gwas, seqnames.field = "chr", start.field = "ps", end.field = "ps", keep.extra.columns = T)
  hits4 <- findOverlaps(sig_gwas_ranges, gwas_sum_ranges)
  olap4 <- pintersect(sig_gwas_ranges[queryHits(hits4)],gwas_sum_ranges[subjectHits(hits4)])
  olap4_df <- as.data.frame(olap4)
  hits5 <- findOverlaps(gwas_sum_ranges, sig_gwas_ranges)
  olap5 <- pintersect(gwas_sum_ranges[queryHits(hits5)],sig_gwas_ranges[subjectHits(hits5)])
  olap5_df <- as.data.frame(olap5)
  test_input_for_plot <- inner_join(olap4_df, olap5_df, by = "start")
  data.list <- list()
  snp.list <- list()
  plot.list <- list()
  unique_loci <- unique(test_input_for_plot$locus)
  for (i in seq_along(unique_loci)) {
    data.list[[i]] <- filter(test_input_for_plot, locus == unique_loci[i])
    snp.list[[i]] <- filter(sig_SNP_info, locus == unique_loci[i])
    p3 <- ggplot(data = data.list[[i]]) + geom_point(aes(start,-log10(p_wald.y))) +
      geom_text(data=data.list[[i]] %>% filter(-log10(p_wald.y) >= 5),
                aes(start, -log10(p_wald.y),label=rs.y)) + theme(text = element_text(size = 6)) + ylab("-Log10 Corrected P-value")
    p4 <- ggplot(data = snp.list[[i]]) +
      geom_linerange(aes(x = gene_name, ymin = start.y, ymax = end.y)) +
      coord_flip() +
      xlab("Gene Name") + theme(text = element_text(size = 6))
    p3b <- p3 + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
    plot.list[[i]] <- p3b + p4 + plot_layout(ncol = 1, heights = c(6, 6))
  }
  pdf(file.path(ResultsPath, "plots.pdf"))
  for (i in seq_along(plot.list)) {
    print(plot.list[[i]])
  }
  dev.off()
  print("All done, happy exploring")
}
