#' locus2gene
#'
#' A function that takes the prepared l2g data from opentargets match them with gene names from feline GWAS
#'
#' @param gene gene name
locus2gene <- function(gene) {
  df <- l2g_annotated_full[grep(paste0("\\b",gene,"\\b"), l2g_annotated_full$gene_name),]
  return(df)
}
