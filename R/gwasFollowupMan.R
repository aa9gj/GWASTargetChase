#' gwasFollowupMan
#'
#' main package extension function that will take GWAS sumstats, a gtf file, and a pvalue of interest to returns information on the genes of interest but also the genetic
#' association file and l2g file data prepped with l2gPrep and geneticAssocPrep functions in this package. Only use this function when you have updated genetic and l2g data
#' from opentargets
#'
#' @param sumStats GWAS summary statistics file. It assumed a ps, and chr columns.
#' @param felGTF GTF file
#' @param pval default is 5*10-8 and you could change but make sure to use standard form
#' @param ResultsPath default is the working directory but you can provide a path of your own ensuring it ends with a /
#' @param impc impc data file path
#' @param assocOT results from geneticAssocPrep function (file path)
#' @param l2gOT results from l2gPrep function (file path)
#' @return Writes g2d_results.txt and l2g_results.txt to ResultsPath
#' @importFrom dplyr filter left_join
#' @importFrom data.table fread
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @importFrom S4Vectors queryHits subjectHits
#' @export gwasFollowupMan

gwasFollowupMan <- function(sumStats, felGTF, pval = 0.00000005, ResultsPath = ".", impc = NULL, assocOT, l2gOT) {
  print("Reading OpenTargets Genetic Association file")
  association_data <- fread(assocOT)
  print("Reading OpenTargets locus to gene file")
  l2g_data <- fread(l2gOT, sep = "\t", quote = "")
  print("Reading IMPC phenotype data")
  impc_data <- fread(impc)
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
  tryCatch(sig_gwas_ranges <- makeGRangesFromDataFrame(sig_gwas, keep.extra.columns = T, seqnames.field = "chr"),
           error = function(e)
             stop("nothing is being reported using this p-value, consider increasing it"))
  cat_pc_gene_ranges <- makeGRangesFromDataFrame(cat_pc_gene_gtf, keep.extra.columns = T)
  hits <- findOverlaps(cat_pc_gene_ranges, sig_gwas_ranges)
  olap <- pintersect(cat_pc_gene_ranges[queryHits(hits)],sig_gwas_ranges[subjectHits(hits)])
  overlap_genes <- as.vector(na.exclude(olap$gene_name))
  genes <- overlap_genes
  print("Genes to Disease GWAS follow-up")
  results_g2d <- list()
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(genes), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for (i in seq_along(genes)) {
    results_g2d[[i]] <- gene2diseaseMan(genes[i], association_data)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_g2d_df <- as.data.frame(do.call("rbind", results_g2d))
  if (nrow(results_g2d_df) > 0) {
    results_g2d_df <- left_join(results_g2d_df, impc_data, by = "gene_name")
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
    results_l2g[[i]] <- locus2geneMan(genes[i], l2g_data)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_l2g_df <- as.data.frame(do.call("rbind", results_l2g))
  if (nrow(results_l2g_df) > 0) {
    results_l2g_df <- left_join(results_l2g_df, impc_data, by = "gene_name")
  }
  print("writing results")
  write.table(results_l2g_df, file.path(ResultsPath, "l2g_results.txt"), quote = F, row.names = F, sep = "\t")
  print("All done, happy exploring")
}
