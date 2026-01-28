#' gwasFollowuptest
#'
#' main package function that will take GWAS sumstats, a gtf file, and a pvalue of interest to returns information on the genes of interest
#'
#' @param sumStats GWAS summary statistics file. It assumed a ps, and chr columns.
#' @param felGTF GTF file
#' @param pval default is 5*10-8 and you could change but make sure to use standard form
#' @param ResultsPath default is the working directory but you can provide a path of your own ensuring it ends with a /
#' @param phenomePath path to phenomeXcan data if available default NULL
#' @param twasPath path to twas data if available default NULL
#' @importFrom dplyr filter left_join bind_rows select_if if_all
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @import  GenomicRanges
#' @importFrom rtracklayer import
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom tibble rownames_to_column
#' @export gwasFollowuptest

gwasFollowuptest <- function(sumStats, felGTF, pval = 0.00000005, ResultsPath = ".",
                             phenomePath = NULL, twasPath = NULL) {
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
  print("Loading IMPC phenotype data")
  tryCatch({
    data("impc", envir = environment())
  }, error = function(e) {
    stop("Required data 'impc' not found. ",
         "Please use gwasFollowupMan() with your own prepared data files, ",
         "or prepare the data using IMPCprep(). See ?IMPCprep for details.")
  })
  if (!is.null(phenomePath)) {
    print("Reading phenomeXcan data")
    phenomeXcan <- fread(phenomePath)
  } else {
    phenomeXcan <- NULL
    print("No phenomeXcan data has been provided")}
  if (!is.null(twasPath)) {
    print("Reading TWAS data")
    twas <- fread(twasPath)
    twas$`ENSG ID` <- gsub("\\..*", "", twas$`ENSG ID`)
    twas <- twas[, c(2:15)]
  } else {
    twas <- NULL
    print("No TWAS data has been provided")}
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
    results_g2d[[i]] <- gene2disease(genes[i])
    setTxtProgressBar(pb, i)
  }
  results_g2d_df <- as.data.frame(do.call("rbind", results_g2d))
  results_g2d_df <- left_join(results_g2d_df, impc, by = "gene_name")
  if (!is.null(phenomeXcan)) {
    print("Wrangling phenomeXcan data, this might take a while")
    pheList <- list()
    for (i in seq_along(genes)){
      pheList[[i]] <- phenomeXcan[grep(genes[i], phenomeXcan$gene_name),]
      pheList[[i]] <- pheList[[i]] %>%
        select_if(~ !any(is.na(.)))
    }
    columns <- list()
    for (g in seq_along(genes)) {
      columns[[g]] <- pheList[[g]]$gene_name
    }
    for (j in seq_along(genes)) {
      pheList[[j]] <- subset(pheList[[j]], select = -gene_name)
      pheList[[j]] <- as.data.frame(t(pheList[[j]]))
      colnames(pheList[[j]]) <- columns[[j]]
      pheList[[j]] <- pheList[[j]] %>% filter(if_all(.fns = ~. >= 0.101))
      pheList[[j]] <- as.data.frame(t(pheList[[j]]))
      pheList[[j]] <- tibble::rownames_to_column(pheList[[j]], "gene_name")
    }
    phenomeXcan <- do.call("bind_rows", pheList)
  }
  if (!is.null(phenomeXcan)) {
    results_g2d_df <- left_join(results_g2d_df, phenomeXcan, by = "gene_name")} else {results_g2d_df <- results_g2d_df}
  if (!is.null(twas)){
    results_g2d_df <- left_join(results_g2d_df, twas, by = c("gene_name" = "Gene"))} else {results_g2d_df<- results_g2d_df}
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
  if (!is.null(phenomeXcan)) {
    results_l2g_df <- left_join(results_l2g_df, phenomeXcan, by = "gene_name")} else {results_l2g_df <- results_l2g_df}
  if (!is.null(twas)){
    results_l2g_df <- left_join(results_l2g_df, twas, by = c("gene_name" = "Gene"))} else {results_l2g_df<- results_l2g_df}
  print("writing results")
  write.table(results_l2g_df, paste0(ResultsPath, "l2g_results.txt"), quote = F, row.names = F, sep = "\t")
  print("All done, happy exploring")
}
