#' locus2gene
#'
#' A function that takes the prepared l2g data from opentargets match them with gene names from feline GWAS
#'
#' @param gene gene name
#' @param l2g_data data.frame of locus-to-gene data
#' @return data.frame of matching L2G entries for the given gene
locus2gene <- function(gene, l2g_data) {
  df <- l2g_data[grep(paste0("\\b",gene,"\\b"), l2g_data$gene_name), , drop = FALSE]
  return(df)
}
