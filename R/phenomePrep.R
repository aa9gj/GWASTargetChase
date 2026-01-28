#' phenomePrep
#'
#' Prepare PhenomeXcan data if present
#'
#' @param phenome path to bulk phenomeXcan data
#' @param gtf path to human gtf file (use recommended)
#' @param ResultsPath default is the working directory but you can provide a path of your own ensuring it ends with a /
#' @importFrom dplyr filter left_join inner_join %>%
#' @importFrom data.table fread
#' @importFrom rtracklayer import
#' @export phenomePrep

phenomePrep <- function(phenome, gtf, ResultsPath = ".") {
  phenome <- as.data.frame(fread(phenome))
  gtf <- as.data.frame(rtracklayer::import(gtf))
  gtf <- filter(gtf, type == "gene", gene_type == "protein_coding")
  gtf <- gtf[, c("gene_id", "gene_name")]
  gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)
  phenome <- inner_join(phenome, gtf, by = "gene_id")
  phenome <- phenome[, -1]
  write.table(phenome, paste0(ResultsPath, "phenomeXcan"), quote = F, row.names = F)
}



