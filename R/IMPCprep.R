#' IMPCprep
#'
#' Function to prepare IMPC data downloaded through the the command line
#'
#' @param impcData path to IMPC phenotype hits per gene file
#' @param ResultsPath default is the working directory but you can provide a path of your own
#' @return Writes impcData file to ResultsPath
#' @importFrom data.table fread
#' @export IMPCprep

IMPCprep <- function(impcData, ResultsPath = ".") {
  impc_phenotype <- fread(impcData)
  impc_phenotype$`Gene Symbol` <- toupper(impc_phenotype$`Gene Symbol`)
  colnames(impc_phenotype) <- c("gene_name", "MGI_Gene_id", "#phenotype_hits", "Phenotype_Hits")
  write.table(impc_phenotype, file.path(ResultsPath, "impcData"), quote = F, row.names = F, col.names = T, sep = "\t")
}


