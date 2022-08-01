#' gwasFollowup
#'
#' main package function that will take GWAS sumstats, a gtf file, and a pvalue of interest to returns information on the genes of interest
#'
#' @param sumStats GWAS summary statistics file. It assumed a ps, and chr columns.
#' @param felGTF GTF file
#' @param pval default is 5*10-8 and you could change but make sure to use standard form
#' @param ResultsPath default is the working directory but you can provide a path of your own ensuring it ends with a /
#' @importFrom dplyr filter left_join inner_join
#' @importFrom data.table fread
#' @import GenomicRanges
#' @import ggplot2
#' @import patchwork
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices pdf
#' @importFrom rtracklayer import
#' @importFrom S4Vectors queryHits subjectHits
#' @export gwasFollowup

gwasFollowup <- function(sumStats, felGTF, pval = 0.00000005, ResultsPath = ".") {
  # Read in GTF file for cats and keep only pc genes
  print("Loading OpenTargets Genetic Association file from this package")
  data("disease_target_genetic_association")
  print("Loading OpenTargets locus to gene file from this package")
  data("l2g_annotated_full")
  print("loading IMPC phenotype data from this package")
  data("impc")
  print("Reading the feline gtf file and filtering for protein coding genes, this might take a while")
  cat_gtf <- as.data.frame(rtracklayer::import(felGTF))
  cat_pc_gene_gtf <- filter(cat_gtf, gene_biotype == "protein_coding", type == "gene")
  print("reading the summary statistics file, filtering based on the input p_value")
  gwas <- fread(sumStats)
  gwas$chr <- gsub("20", "X", gwas$chr)
  sig_gwas <- filter(gwas, p_wald <= pval)
  print("Overlap genes within significant loci (1Mb)")
  sig_gwas$start <- sig_gwas$ps - 1000000
  sig_gwas$end <- sig_gwas$ps + 1000000
  sig_gwas$locus <- rep("x", nrow(sig_gwas))
  for (i in seq_along(sig_gwas$locus)) {
    sig_gwas$locus[i] <- paste0("locus_", i)
  }
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
  olap3_df <- olap3_df[, c(1:5,24)]
  olap3_df$join <- paste0(olap3_df$seqnames, "_", olap3_df$start, "_", olap3_df$end)
  sig_SNP_info <- inner_join(olap2_df, olap3_df, by = "join")
  sig_SNP_info <- sig_SNP_info[,c(6:17,20:22,25)]
  print("Genes to Disease GWAS follow-up")
  results_g2d <- list()
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(genes), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for (i in seq_along(genes)) {
    results_g2d[[i]] <- gene2disease(genes[i])
    setTxtProgressBar(pb, i)
  }
  results_g2d_df <- as.data.frame(do.call("rbind", results_g2d))
  results_g2d_df <- left_join(results_g2d_df, impc, by = "gene_name")
  results_g2d_df <- inner_join(results_g2d_df, sig_SNP_info, by = "gene_name")
  print("writing results")
  write.table(results_g2d_df, paste0(ResultsPath, "g2d_results.txt"), quote = F, row.names = F, sep = "\t")
  print("Locus to Genes GWAS follow-up")
  results_l2g <- list()
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(genes), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for (i in seq_along(genes)) {
    results_l2g[[i]] <- locus2gene(genes[i])
    setTxtProgressBar(pb, i)
  }
  results_l2g_df <- as.data.frame(do.call("rbind", results_l2g))
  results_l2g_df <- left_join(results_l2g_df, impc, by = "gene_name")
  results_l2g_df <- inner_join(results_l2g_df, sig_SNP_info, by = "gene_name")
  results_l2g_filter <- results_l2g_df[,c(1,6,12:16,21,27,36,39:51)]
  results_l2g_filter$gene_cards <- paste0("https://www.genecards.org/Search/Keyword?queryString=", results_l2g_filter$gene_name)
  names(results_l2g_filter)[names(results_l2g_filter) == 'y_proba_full_model'] <- 'full_l2g_score'
  names(results_l2g_filter)[names(results_l2g_filter) == 'Phenotype_Hits'] <- 'IMPC_results'
  print("writing results")
  write.table(results_l2g_df, paste0(ResultsPath, "l2g_results.txt"), quote = F, row.names = F, sep = "\t")
  write.table(results_l2g_filter, paste0(ResultsPath, "l2g_filtered_results.txt"), quote = F, row.names = F, sep = "\t")
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
  for (i in 1:length(unique(test_input_for_plot$locus))) {
    data.list[[i]] <- filter(test_input_for_plot, locus == paste0("locus_", i))
    snp.list[[i]] <- filter(sig_SNP_info, locus == paste0("locus_", i))
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
  setwd(ResultsPath)
  pdf("plots.pdf")
  for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
  }
  dev.off()
  print("All done, happy exploring")
}
