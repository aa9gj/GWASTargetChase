#' IMPCprep
#'
#' Function to prepare IMPC data downloaded through the the command line
#'
#' @importFrom data.table fread
#' @export IMPCprep

IMPCprep <- function(impcData, ResultsPath = ".") {
  impc_phenotype <- fread(impcData)
  impc_phenotype$`Gene Symbol` <- toupper(impc_phenotype$`Gene Symbol`)
  colnames(impc_phenotype) <- c("gene_name", "MGI_Gene_id", "#phenotype_hits", "Phenotype_Hits")
  write.table(impc_phenotype, paste0(ResultsPath, "impcData"), quote = F, row.names = F, col.names = T, sep = "\t")
}


