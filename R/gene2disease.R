#' gene2disease
#'
#' A function that takes the prepared disease target genetic association and match them with gene names from feline GWAS
#'
#' @param gene gene name
gene2disease <- function(gene) {
  df <- disease_target_genetic_association[grep(paste0("\\b",gene,"\\b"), disease_target_genetic_association$gene_name)]
  return(df)
}
