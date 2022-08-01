#' locus2geneMan
#'
#' A function that takes the prepared l2g data from opentargets match them with gene names from feline GWAS
#'
#' @param gene gene name
locus2geneMan <- function(gene) {
  df <- l2g_file[grep(paste0("\\b",gene,"\\b"), l2g_file$gene_name),]
  return(df)
}
