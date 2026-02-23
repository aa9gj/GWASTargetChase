#' Load Zoonomia Orthology Data
#'
#' Loads a Zoonomia orthology file and extracts gene symbols from the
#' t_transcript column (everything after the first dot).
#'
#' @param target_species Target species: "human" or "mouse"
#' @param query_species Query species: "cat" or "dog"
#' @param zoo_dir Directory containing zoonomia files. If NULL, uses package extdata.
#' @return data.frame with columns: gene_symbol, t_gene, orthology_class, V1, V3, hills_grade
#' @importFrom data.table fread
#' @export
load_zoonomia_orthologs <- function(target_species, query_species, zoo_dir = NULL) {
  zoo_filename <- paste0("zoo_", target_species, "_", query_species, ".txt")

  if (is.null(zoo_dir)) {
    zoo_file <- system.file("extdata", zoo_filename, package = "GWASTargetChase")
  } else {
    zoo_file <- file.path(zoo_dir, zoo_filename)
  }

  if (!file.exists(zoo_file) || zoo_file == "") {
    stop("Zoonomia orthology file not found: ", zoo_filename, call. = FALSE)
  }

  zoo_data <- data.table::fread(zoo_file)

  # Extract gene symbol from t_transcript (everything after the first dot)
  zoo_data$gene_symbol <- sub("^[^.]+\\.", "", zoo_data$t_transcript)

  # Keep relevant columns and deduplicate per gene symbol
  result <- unique(as.data.frame(
    zoo_data[, c("gene_symbol", "t_gene", "orthology_class", "V1", "V3", "hills_grade")]
  ))

  return(result)
}

#' Translate Gene Names Using Zoonomia Orthology
#'
#' Translates gene names between species using Zoonomia orthology data.
#' For non-human GWAS, use target_species="human" to get human orthologs
#' for OpenTargets, or target_species="mouse" for IMPC lookups.
#' Gene symbols are extracted from the t_transcript column of the zoonomia
#' file (everything after the first dot).
#'
#' @param genes Character vector of gene names to translate
#' @param target_species Target species: "human" or "mouse"
#' @param query_species Query species: "cat" or "dog"
#' @param zoo_dir Directory containing zoonomia files. If NULL, uses package extdata.
#' @return data.frame with columns: original_gene, target_gene, orthology_class, V1, V3, hills_grade
#' @export
translate_genes <- function(genes, target_species, query_species, zoo_dir = NULL) {
  zoo <- load_zoonomia_orthologs(target_species, query_species, zoo_dir)

  # Case-insensitive matching
  zoo$gene_symbol_upper <- toupper(zoo$gene_symbol)
  genes_upper <- toupper(genes)

  matched <- zoo[zoo$gene_symbol_upper %in% genes_upper, ]

  if (nrow(matched) > 0) {
    result <- data.frame(
      original_gene = matched$gene_symbol_upper,
      target_gene = matched$gene_symbol,
      orthology_class = matched$orthology_class,
      V1 = matched$V1,
      V3 = matched$V3,
      hills_grade = matched$hills_grade,
      stringsAsFactors = FALSE
    )
    result <- unique(result)
    return(result)
  } else {
    return(data.frame(
      original_gene = character(),
      target_gene = character(),
      orthology_class = character(),
      V1 = character(),
      V3 = character(),
      hills_grade = character(),
      stringsAsFactors = FALSE
    ))
  }
}

#' Get Zoonomia Orthology Info for Results
#'
#' Adds orthology metadata columns (orthology_class, V1, V3, hills_grade)
#' to a results data frame by matching on gene_name.
#'
#' @param results_df data.frame with a gene_name column
#' @param ortho_df data.frame from translate_genes() with orthology metadata
#' @return results_df with added orthology columns
#' @keywords internal
add_orthology_info <- function(results_df, ortho_df) {
  if (is.null(ortho_df) || nrow(ortho_df) == 0 || nrow(results_df) == 0) {
    return(results_df)
  }

  ortho_info <- unique(ortho_df[, c("target_gene", "orthology_class", "V1", "V3", "hills_grade")])
  names(ortho_info)[1] <- "gene_name"

  merged <- merge(results_df, ortho_info, by = "gene_name", all.x = TRUE)
  return(merged)
}
