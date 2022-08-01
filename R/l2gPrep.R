#' l2gPrep
#'
#' Preparation of OpenTarget l2g gwasFollowupMan function
#'
#' @param l2gOT l2g Data full path (see external documentation for data download from OpenTargets).
#' @param studyOT study index data full path(see external documentation for data download from OpenTargets)
#' @param humanGTF human GTF file from GENCODE
#' @param ResultsPath default is the working directory but you can provide a path of your own ensuring it ends with a /
#' @import arrow
#' @importFrom dplyr filter inner_join left_join
#' @importFrom data.table fread
#' @importFrom rtracklayer import
#' @import jsonlite
#' @export l2gPrep

l2gPrep <- function(l2gOT, studyOT, humanGTF, ResultsPath = ".") {
  DS <- arrow::open_dataset(sources = l2gOT)
  SO <- Scanner$create(DS)
  AT <- SO$ToTable()
  l2g_data <- as.data.frame(AT)
  study_index <- stream_in(file(studyOT))
  study_index <- study_index[,c(1,4,7:11,14)]
  print("Wrangling data")
  l2g_data_annotated_studies <- left_join(l2g_data, study_index, by = "study_id")
  l2g_data_annotated_studies <- do.call(data.frame, l2g_data_annotated_studies)
  print("Reading human gtf file")
  gencode <- as.data.frame(rtracklayer::import(humanGTF))
  print("wrangling data")
  gencode_pc <- filter(gencode, gene_type == "protein_coding", type == "gene")
  gencode_pc$gene_id <- gsub("\\..*", "", gencode_pc$gene_id)
  gencode_pc <- gencode_pc[,-c(8,9,13:26)]
  l2g_annotated <- inner_join(l2g_data_annotated_studies, gencode_pc, by = "gene_id")
  print("Writing restuls, this might take a while")
  write.table(l2g_annotated, paste0(ResultsPath, "l2g_annotated_full"), quote = F, row.names = F, sep = "\t")
  print("that's a lot of data, all done!")
}
