#' gene2diseaseMan
#'
#' A function that takes the prepared disease target genetic association and match them with gene names from feline GWAS
#'
#' @param gene gene name
gene2diseaseMan <- function(gene) {
  df <- association_file[grep(paste0("\\b",gene,"\\b"), association_file$gene_name)]
  return(df)
}
