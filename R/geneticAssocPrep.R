#' geneticAssocPrep
#'
#' Preparation of OpenTarget data evidence by data source file for gwasFollowupMan function
#'
#' @param assocOT evidence by data source file full path (see external documentation for data download from OpenTargets).
#' @param diseaseOT Disease data full path(see external documentation for data download from OpenTargets)
#' @param humanGTF human GTF file from GENCODE
#' @param ResultsPath default is the working directory but you can provide a path of your own ensuring it ends with a /
#' @importFrom dplyr filter inner_join
#' @importFrom data.table fread
#' @importFrom rtracklayer import
#' @import jsonlite
#' @export geneticAssocPrep

geneticAssocPrep <- function(assocOT, diseaseOT, humanGTF, ResultsPath = ".") {
  print("Reading in OpenTarget evidence data, this might take a while")
  assocOT <- stream_in(file(assocOT))
  print("Reading in OpenTarget disease annotation data")
  diseaseOT <- stream_in(file(diseaseOT))
  diseaseOT <- diseaseOT[,c(1,5)]
  print("Preparing data")
  part1_ot <- inner_join(diseaseOT, assocOT, by = c("id" = "diseaseId"))
  gencode <- as.data.frame(rtracklayer::import(humanGTF))
  gencode_pc <- filter(gencode, gene_type == "protein_coding", type == "gene")
  gencode_pc$gene_id <- gsub("\\..*", "", gencode_pc$gene_id)
  gencode_pc <- gencode_pc[,-c(8,9,13:26)]
  part2_ot <- inner_join(part1_ot, gencode_pc, by = c("targetId" = "gene_id"))
  disease_target_associations <- filter(part2_ot, datatypeId == "genetic_association")
  print("writing results")
  write.table(disease_target_associations, file.path(ResultsPath, "disease_target_genetic_associations"), quote = F, row.names = F, sep = "\t", col.names = T)
  print("All done!!")
}
