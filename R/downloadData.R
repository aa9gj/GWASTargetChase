#' downloadData
#'
#' Downloads and prepares real data from IMPC, OpenTargets Platform, and
#' OpenTargets Genetics for use with GWASTargetChase. After running this
#' function, use the output files with \code{gwasFollowupMan()}.
#'
#' This function downloads:
#' \itemize{
#'   \item IMPC phenotype hits per gene
#'   \item OpenTargets Platform genetic association evidence and disease annotations
#'   \item OpenTargets Genetics locus-to-gene (L2G) scores and study index
#'   \item GENCODE human GTF (needed to map gene IDs to gene names)
#' }
#'
#' The downloads are large (several GB total). Ensure you have sufficient disk
#' space and a stable internet connection.
#'
#' @param destdir Directory to save downloaded and processed files. Will be
#'   created if it does not exist. Default is "GWASTargetChase_data" in the
#'   current working directory.
#' @param OT_platform_version OpenTargets Platform release version (e.g. "24.06").
#'   Check \url{https://platform.opentargets.org/downloads} for available versions.
#' @param gencode_version GENCODE human release version (e.g. "44").
#'   Check \url{https://www.gencodegenes.org/human/} for available versions.
#' @param skip_existing If TRUE (default), skip files that have already been downloaded.
#' @param process If TRUE (default), process the raw downloads into analysis-ready
#'   files after downloading. Set to FALSE to only download without processing.
#' @return Invisibly returns a list of output file paths.
#'
#' @details
#' After running this function, use the output files with \code{gwasFollowupMan}:
#' \preformatted{
#'   paths <- downloadData(destdir = "my_data")
#'   gwasFollowupMan(
#'     sumStats = "my_gwas_sumstats.txt",
#'     felGTF   = "Felis_catus.Felis_catus_9.0.106.gtf",
#'     impc     = paths$impc,
#'     assocOT  = paths$assoc,
#'     l2gOT    = paths$l2g
#'   )
#' }
#'
#' @importFrom data.table fread
#' @importFrom dplyr filter inner_join left_join
#' @import jsonlite
#' @import arrow
#' @importFrom rtracklayer import
#' @export downloadData

downloadData <- function(destdir = "GWASTargetChase_data",
                         OT_platform_version = "24.06",
                         gencode_version = "44",
                         skip_existing = TRUE,
                         process = TRUE) {

  # Create directories
  raw_dir <- file.path(destdir, "raw")
  processed_dir <- file.path(destdir, "processed")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

  output_paths <- list()

  # ============================================================================
  # 1. IMPC Data
  # ============================================================================
  cat("==================================================\n")
  cat("1/4: Downloading IMPC phenotype data...\n")
  cat("==================================================\n")

  impc_url <- "https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/phenotypeHitsPerGene.csv.gz"
  impc_gz <- file.path(raw_dir, "phenotypeHitsPerGene.csv.gz")
  impc_csv <- file.path(raw_dir, "phenotypeHitsPerGene.csv")

  if (skip_existing && file.exists(impc_csv)) {
    cat("  Already exists, skipping download.\n")
  } else {
    cat("  Downloading from IMPC FTP...\n")
    tryCatch({
      download.file(impc_url, impc_gz, mode = "wb", quiet = FALSE)
      R.utils::gunzip(impc_gz, remove = TRUE, overwrite = TRUE)
      cat("  IMPC data downloaded successfully.\n")
    }, error = function(e) {
      # Try without R.utils
      tryCatch({
        system2("gunzip", args = c("-f", impc_gz))
        cat("  IMPC data downloaded and extracted.\n")
      }, error = function(e2) {
        cat("  WARNING: Could not extract .gz file. Install R.utils package or gunzip manually:\n")
        cat("    R.utils::gunzip('", impc_gz, "')\n")
        cat("  Or from command line: gunzip ", impc_gz, "\n")
      })
    })
  }
  cat("\n")

  # ============================================================================
  # 2. OpenTargets Platform - Association Evidence + Diseases
  # ============================================================================
  cat("==================================================\n")
  cat("2/4: Downloading OpenTargets Platform data...\n")
  cat("     Version: ", OT_platform_version, "\n")
  cat("==================================================\n")

  # --- 2a. Association by datasource ---
  assoc_dir <- file.path(raw_dir, "associationByDatasourceDirect")
  assoc_json <- file.path(raw_dir, "associationByDatasourceDirect.json")

  if (skip_existing && file.exists(assoc_json)) {
    cat("  Association data already exists, skipping.\n")
  } else {
    assoc_url <- paste0("http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/",
                        OT_platform_version,
                        "/output/etl/json/associationByDatasourceDirect/")
    cat("  Downloading association evidence data...\n")
    cat("  Source: ", assoc_url, "\n")
    cat("  This is a large download with multiple part files.\n")
    dir.create(assoc_dir, showWarnings = FALSE)
    .download_ot_parts(assoc_url, assoc_dir, assoc_json)
  }

  # --- 2b. Disease annotations ---
  disease_dir <- file.path(raw_dir, "diseases")
  disease_json <- file.path(raw_dir, "diseases.json")

  if (skip_existing && file.exists(disease_json)) {
    cat("  Disease data already exists, skipping.\n")
  } else {
    disease_url <- paste0("http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/",
                          OT_platform_version,
                          "/output/etl/json/diseases/")
    cat("  Downloading disease annotation data...\n")
    cat("  Source: ", disease_url, "\n")
    dir.create(disease_dir, showWarnings = FALSE)
    .download_ot_parts(disease_url, disease_dir, disease_json)
  }
  cat("\n")

  # ============================================================================
  # 3. OpenTargets Genetics - L2G + Study Index
  # ============================================================================
  cat("==================================================\n")
  cat("3/4: Downloading OpenTargets Genetics data...\n")
  cat("==================================================\n")

  # --- 3a. L2G data (Parquet format) ---
  l2g_dir <- file.path(raw_dir, "l2g")
  if (skip_existing && dir.exists(l2g_dir) && length(list.files(l2g_dir, pattern = "\\.parquet$")) > 0) {
    cat("  L2G data already exists, skipping.\n")
  } else {
    l2g_url <- "http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/l2g/"
    cat("  Downloading L2G parquet data...\n")
    cat("  Source: ", l2g_url, "\n")
    dir.create(l2g_dir, showWarnings = FALSE)
    .download_ot_parts(l2g_url, l2g_dir, output_combined = NULL)
  }

  # --- 3b. Study index ---
  study_dir <- file.path(raw_dir, "study-index")
  study_json <- file.path(raw_dir, "study_index.json")

  if (skip_existing && file.exists(study_json)) {
    cat("  Study index already exists, skipping.\n")
  } else {
    study_url <- "http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/lut/study-index/"
    cat("  Downloading study index data...\n")
    cat("  Source: ", study_url, "\n")
    dir.create(study_dir, showWarnings = FALSE)
    .download_ot_parts(study_url, study_dir, study_json)
  }
  cat("\n")

  # ============================================================================
  # 4. GENCODE Human GTF
  # ============================================================================
  cat("==================================================\n")
  cat("4/4: Downloading GENCODE human GTF...\n")
  cat("     Version: ", gencode_version, "\n")
  cat("==================================================\n")

  gtf_name <- paste0("gencode.v", gencode_version, ".annotation.gtf")
  gtf_gz <- file.path(raw_dir, paste0(gtf_name, ".gz"))
  gtf_file <- file.path(raw_dir, gtf_name)

  if (skip_existing && file.exists(gtf_file)) {
    cat("  GENCODE GTF already exists, skipping.\n")
  } else {
    gtf_url <- paste0("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
                       gencode_version, "/", gtf_name, ".gz")
    cat("  Downloading from: ", gtf_url, "\n")
    tryCatch({
      download.file(gtf_url, gtf_gz, mode = "wb", quiet = FALSE)
      tryCatch(
        R.utils::gunzip(gtf_gz, remove = TRUE, overwrite = TRUE),
        error = function(e) system2("gunzip", args = c("-f", gtf_gz))
      )
      cat("  GENCODE GTF downloaded successfully.\n")
    }, error = function(e) {
      cat("  ERROR downloading GTF: ", conditionMessage(e), "\n")
    })
  }
  cat("\n")

  # ============================================================================
  # Process raw data into analysis-ready files
  # ============================================================================
  if (!process) {
    cat("Skipping processing (process = FALSE).\n")
    cat("Raw files are in: ", raw_dir, "\n")
    cat("Run downloadData() again with process = TRUE to process.\n")
    return(invisible(list(raw_dir = raw_dir)))
  }

  cat("==================================================\n")
  cat("Processing downloaded data...\n")
  cat("==================================================\n")

  # --- Process IMPC ---
  impc_output <- file.path(processed_dir, "impcData")
  if (file.exists(impc_csv)) {
    cat("  Processing IMPC data...\n")
    impc_phenotype <- fread(impc_csv)
    impc_phenotype[[1]] <- toupper(impc_phenotype[[1]])
    colnames(impc_phenotype) <- c("gene_name", "MGI_Gene_id", "num_phenotype_hits", "Phenotype_Hits")
    write.table(impc_phenotype, impc_output, quote = FALSE, row.names = FALSE,
                col.names = TRUE, sep = "\t")
    output_paths$impc <- impc_output
    cat("  Saved: ", impc_output, "\n")
  } else {
    cat("  WARNING: IMPC CSV not found, skipping.\n")
  }

  # --- Process OpenTargets Genetic Associations ---
  assoc_output <- file.path(processed_dir, "disease_target_genetic_associations")
  if (file.exists(assoc_json) && file.exists(disease_json) && file.exists(gtf_file)) {
    cat("  Processing OpenTargets association data (this may take several minutes)...\n")
    cat("    Reading association evidence...\n")
    assoc_data <- stream_in(file(assoc_json), verbose = FALSE)
    cat("    Reading disease annotations...\n")
    disease_data <- stream_in(file(disease_json), verbose = FALSE)
    disease_data <- disease_data[, c(1, 5)]
    cat("    Joining evidence with diseases...\n")
    part1_ot <- inner_join(disease_data, assoc_data, by = c("id" = "diseaseId"))
    cat("    Reading GENCODE GTF...\n")
    gencode <- as.data.frame(rtracklayer::import(gtf_file))
    gencode_pc <- filter(gencode, gene_type == "protein_coding", type == "gene")
    gencode_pc$gene_id <- gsub("\\..*", "", gencode_pc$gene_id)
    gencode_pc <- gencode_pc[, -c(8, 9, 13:26)]
    cat("    Intersecting with protein-coding genes...\n")
    part2_ot <- inner_join(part1_ot, gencode_pc, by = c("targetId" = "gene_id"))
    disease_target_associations <- filter(part2_ot, datatypeId == "genetic_association")
    write.table(disease_target_associations, assoc_output, quote = FALSE,
                row.names = FALSE, sep = "\t", col.names = TRUE)
    output_paths$assoc <- assoc_output
    cat("  Saved: ", assoc_output, " (", nrow(disease_target_associations), " rows)\n")
  } else {
    missing <- c()
    if (!file.exists(assoc_json)) missing <- c(missing, "association JSON")
    if (!file.exists(disease_json)) missing <- c(missing, "disease JSON")
    if (!file.exists(gtf_file)) missing <- c(missing, "GENCODE GTF")
    cat("  WARNING: Skipping association processing. Missing: ", paste(missing, collapse = ", "), "\n")
  }

  # --- Process OpenTargets L2G ---
  l2g_output <- file.path(processed_dir, "l2g_annotated_full")
  if (dir.exists(l2g_dir) && length(list.files(l2g_dir, pattern = "\\.parquet$")) > 0 &&
      file.exists(study_json) && file.exists(gtf_file)) {
    cat("  Processing OpenTargets L2G data (this may take several minutes)...\n")
    cat("    Reading L2G parquet files...\n")
    DS <- arrow::open_dataset(sources = l2g_dir)
    SO <- Scanner$create(DS)
    AT <- SO$ToTable()
    l2g_data <- as.data.frame(AT)
    cat("    Reading study index...\n")
    study_index <- stream_in(file(study_json), verbose = FALSE)
    study_index <- study_index[, c(1, 4, 7:11, 14)]
    cat("    Merging L2G with study annotations...\n")
    l2g_annotated <- left_join(l2g_data, study_index, by = "study_id")
    l2g_annotated <- do.call(data.frame, l2g_annotated)
    cat("    Reading GENCODE GTF...\n")
    if (!exists("gencode_pc")) {
      gencode <- as.data.frame(rtracklayer::import(gtf_file))
      gencode_pc <- filter(gencode, gene_type == "protein_coding", type == "gene")
      gencode_pc$gene_id <- gsub("\\..*", "", gencode_pc$gene_id)
      gencode_pc <- gencode_pc[, -c(8, 9, 13:26)]
    }
    cat("    Intersecting with protein-coding genes...\n")
    l2g_final <- inner_join(l2g_annotated, gencode_pc, by = "gene_id")
    write.table(l2g_final, l2g_output, quote = FALSE, row.names = FALSE, sep = "\t")
    output_paths$l2g <- l2g_output
    cat("  Saved: ", l2g_output, " (", nrow(l2g_final), " rows)\n")
  } else {
    missing <- c()
    if (!dir.exists(l2g_dir) || length(list.files(l2g_dir, pattern = "\\.parquet$")) == 0)
      missing <- c(missing, "L2G parquet files")
    if (!file.exists(study_json)) missing <- c(missing, "study index JSON")
    if (!file.exists(gtf_file)) missing <- c(missing, "GENCODE GTF")
    cat("  WARNING: Skipping L2G processing. Missing: ", paste(missing, collapse = ", "), "\n")
  }

  # ============================================================================
  # Summary
  # ============================================================================
  cat("\n")
  cat("==================================================\n")
  cat("  Download and processing complete!\n")
  cat("==================================================\n")
  cat("\n")
  cat("Processed files in: ", processed_dir, "\n")
  if (length(output_paths) > 0) {
    for (nm in names(output_paths)) {
      fsize <- file.size(output_paths[[nm]])
      cat("  ", nm, ": ", output_paths[[nm]],
          " (", format(fsize, big.mark = ","), " bytes)\n")
    }
  }
  cat("\n")
  cat("Usage with gwasFollowupMan():\n")
  cat("  gwasFollowupMan(\n")
  cat("    sumStats = 'your_gwas_sumstats.txt',\n")
  cat("    felGTF   = 'your_species.gtf',\n")
  if (!is.null(output_paths$impc))
    cat("    impc     = '", output_paths$impc, "',\n", sep = "")
  if (!is.null(output_paths$assoc))
    cat("    assocOT  = '", output_paths$assoc, "',\n", sep = "")
  if (!is.null(output_paths$l2g))
    cat("    l2gOT    = '", output_paths$l2g, "'\n", sep = "")
  cat("  )\n")

  invisible(output_paths)
}


# Internal helper: download part files from an OpenTargets FTP directory
# and optionally concatenate JSON parts into one file.
.download_ot_parts <- function(base_url, dest_dir, output_combined = NULL) {
  # Try wget first (most reliable for recursive FTP)
  wget_available <- Sys.which("wget") != ""
  curl_available <- Sys.which("curl") != ""

  if (wget_available) {
    cat("  Using wget for recursive download...\n")
    cmd <- paste0("wget -r -np -nH --cut-dirs=8 -P '", dest_dir, "' '", base_url, "' 2>&1")
    result <- system(cmd, intern = TRUE)
    if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
      cat("  wget failed. Trying alternative method...\n")
    }
  } else if (curl_available) {
    cat("  wget not found. Trying curl...\n")
    # curl can't do recursive FTP easily; try to list and download
    cat("  NOTE: curl cannot do recursive FTP downloads.\n")
    cat("  Please install wget or download manually from:\n")
    cat("    ", base_url, "\n")
    cat("  Save files to: ", dest_dir, "\n")
    return(invisible(NULL))
  } else {
    cat("  Neither wget nor curl found.\n")
    cat("  Please download manually from:\n")
    cat("    ", base_url, "\n")
    cat("  Save files to: ", dest_dir, "\n")
    return(invisible(NULL))
  }

  # Concatenate JSON parts if requested
  if (!is.null(output_combined)) {
    json_files <- list.files(dest_dir, pattern = "\\.(json|parquet)$", full.names = TRUE)
    if (length(json_files) > 0 && grepl("\\.json$", output_combined)) {
      cat("  Concatenating ", length(json_files), " part files...\n")
      # Use system cat for efficiency with large files
      cmd <- paste0("cat '", paste(json_files, collapse = "' '"), "' > '", output_combined, "'")
      system(cmd)
      cat("  Combined file: ", output_combined, "\n")
    } else if (length(json_files) == 0) {
      # Check subdirectories (wget sometimes creates nested dirs)
      all_files <- list.files(dest_dir, pattern = "\\.(json|parquet)$",
                              full.names = TRUE, recursive = TRUE)
      if (length(all_files) > 0 && grepl("\\.json$", output_combined)) {
        cat("  Found ", length(all_files), " files in subdirectories, concatenating...\n")
        cmd <- paste0("cat '", paste(all_files, collapse = "' '"), "' > '", output_combined, "'")
        system(cmd)
        cat("  Combined file: ", output_combined, "\n")
      } else {
        cat("  WARNING: No part files found in ", dest_dir, "\n")
        cat("  You may need to download manually from: ", base_url, "\n")
      }
    }
  }

  invisible(NULL)
}
