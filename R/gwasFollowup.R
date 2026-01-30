#' gwasFollowup
#'
#' main package function that will take GWAS sumstats, a gtf file, and a pvalue of interest to returns information on the genes of interest
#'
#' @param sumStats GWAS summary statistics file. It assumed a ps, and chr columns.
#' @param felGTF GTF file
#' @param pval default is 5*10-8 and you could change but make sure to use standard form
#' @param ResultsPath default is the working directory but you can provide a path of your own ensuring it ends with a /
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

gwasFollowup <- function(sumStats, felGTF, pval = 0.00000005, ResultsPath = ".") {
  # Read in GTF file for cats and keep only pc genes
  print("Loading OpenTargets Genetic Association file from this package")
  tryCatch({
    data("disease_target_genetic_association", envir = environment())
  }, error = function(e) {
    stop("Required data 'disease_target_genetic_association' not found. ",
         "Please use gwasFollowupMan() with your own prepared data files, ",
         "or prepare the data using geneticAssocPrep(). See ?geneticAssocPrep for details.")
  })
  print("Loading OpenTargets locus to gene file from this package")
  tryCatch({
    data("l2g_annotated_full", envir = environment())
  }, error = function(e) {
    stop("Required data 'l2g_annotated_full' not found. ",
         "Please use gwasFollowupMan() with your own prepared data files, ",
         "or prepare the data using l2gPrep(). See ?l2gPrep for details.")
  })
  print("loading IMPC phenotype data from this package")
  tryCatch({
    data("impc", envir = environment())
  }, error = function(e) {
    stop("Required data 'impc' not found. ",
         "Please use gwasFollowupMan() with your own prepared data files, ",
         "or prepare the data using IMPCprep(). See ?IMPCprep for details.")
  })
  print("Reading the feline gtf file and filtering for protein coding genes, this might take a while")
  cat_gtf <- as.data.frame(rtracklayer::import(felGTF))
  cat_pc_gene_gtf <- filter(cat_gtf, gene_biotype == "protein_coding", type == "gene")
  print("reading the summary statistics file, filtering based on the input p_value")
  gwas <- fread(sumStats)
  gwas$chr <- gsub("^20$", "X", gwas$chr)
  sig_gwas <- filter(gwas, p_wald <= pval)
  print("Overlap genes within significant loci (1Mb)")
  sig_gwas$start <- sig_gwas$ps - 1000000
  sig_gwas$end <- sig_gwas$ps + 1000000
  sig_gwas$locus <- paste0("locus_", seq_len(nrow(sig_gwas)))
  tryCatch(sig_gwas_ranges <- makeGRangesFromDataFrame(sig_gwas, keep.extra.columns = T, seqnames.field = "chr"),
  error = function(e)
    stop("nothing is being reported using this p-value, consider increasing it"))
  cat_pc_gene_ranges <- makeGRangesFromDataFrame(cat_pc_gene_gtf, keep.extra.columns = T)
  hits <- findOverlaps(cat_pc_gene_ranges, sig_gwas_ranges)
  olap <- pintersect(cat_pc_gene_ranges[queryHits(hits)],sig_gwas_ranges[subjectHits(hits)])
  overlap_genes <- as.vector(na.exclude(olap$gene_name))
  genes <- overlap_genes
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
                       max = length(genes), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for (i in seq_along(genes)) {
    results_g2d[[i]] <- gene2disease(genes[i], disease_target_genetic_association)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_g2d_df <- as.data.frame(do.call("rbind", results_g2d))
  if (nrow(results_g2d_df) > 0) {
    results_g2d_df <- left_join(results_g2d_df, impc, by = "gene_name")
    if ("gene_name" %in% colnames(sig_SNP_info)) {
      results_g2d_df <- inner_join(results_g2d_df, sig_SNP_info, by = "gene_name")
    }
  }
  print("writing results")
  write.table(results_g2d_df, file.path(ResultsPath, "g2d_results.txt"), quote = F, row.names = F, sep = "\t")
  print("Locus to Genes GWAS follow-up")
  results_l2g <- list()
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(genes), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for (i in seq_along(genes)) {
    results_l2g[[i]] <- locus2gene(genes[i], l2g_annotated_full)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_l2g_df <- as.data.frame(do.call("rbind", results_l2g))
  if (nrow(results_l2g_df) > 0) {
    results_l2g_df <- left_join(results_l2g_df, impc, by = "gene_name")
    if ("gene_name" %in% colnames(sig_SNP_info)) {
      results_l2g_df <- inner_join(results_l2g_df, sig_SNP_info, by = "gene_name")
    }
    # Create filtered L2G results with key columns
    available_cols <- colnames(results_l2g_df)
    filter_cols <- intersect(c("gene_name", "study_id", "variant_id", "gene_id",
                                "y_proba_full_model", "y_proba_logi_distance",
                                "y_proba_logi_interaction", "y_proba_logi_molecularQTL",
                                "y_proba_logi_pathogenicity", "trait_reported",
                                "trait_category", "Phenotype_Hits",
                                "rs", "ps", "p_wald", "locus"), available_cols)
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
