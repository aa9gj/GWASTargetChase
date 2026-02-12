#' gene2diseaseMan
#'
#' A function that takes the prepared disease target genetic association and match them with gene names from feline GWAS
#'
#' @param gene gene name
#' @param assoc_data data.frame of disease-target genetic associations
#' @return data.frame of matching associations for the given gene
gene2diseaseMan <- function(gene, assoc_data) {
  df <- assoc_data[grep(paste0("\\b",gene,"\\b"), assoc_data$gene_name), , drop = FALSE]
  return(df)
}
