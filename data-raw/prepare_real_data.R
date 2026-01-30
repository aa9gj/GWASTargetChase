# prepare_real_data.R
# ===================
# Processes downloaded IMPC, OpenTargets Platform, and OpenTargets Genetics
# data into .rda files for the GWASTargetChase package.
#
# Prerequisites:
#   1. Run download_data.sh first (or download files manually)
#   2. Required R packages: data.table, dplyr, jsonlite, arrow, rtracklayer
#
# Usage:
#   setwd("path/to/GWASTargetChase")
#   source("data-raw/prepare_real_data.R")
#
# This script will:
#   - Read and process the raw downloaded data
#   - Use the package's prep functions (or equivalent logic)
#   - Save the final .rda files to data/

library(data.table)
library(dplyr)
library(jsonlite)
library(arrow)
library(rtracklayer)

DOWNLOAD_DIR <- "data-raw/downloaded_data"

# Check that download directory exists
if (!dir.exists(DOWNLOAD_DIR)) {
  stop("Download directory not found: ", DOWNLOAD_DIR, "\n",
       "Run data-raw/download_data.sh first, or set DOWNLOAD_DIR to your download path.")
}

# Create data/ directory
if (!dir.exists("data")) {
  dir.create("data")
}

# ==============================================================================
# 1. IMPC Data
# ==============================================================================
cat("=== Processing IMPC data ===\n")

impc_path <- file.path(DOWNLOAD_DIR, "phenotypeHitsPerGene.csv")
if (!file.exists(impc_path)) {
  stop("IMPC file not found: ", impc_path, "\n",
       "Download it from: https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/phenotypeHitsPerGene.csv.gz")
}

impc_raw <- fread(impc_path)
cat("  Raw IMPC data:", nrow(impc_raw), "genes\n")

# The IMPC file has columns: Gene Symbol, MGI Gene Id, # Phenotype Hits, Phenotype Hits
# Standardize column names and convert gene symbols to uppercase
impc_raw[[1]] <- toupper(impc_raw[[1]])
colnames(impc_raw) <- c("gene_name", "MGI_Gene_id", "num_phenotype_hits", "Phenotype_Hits")

# For the package internal data, we only need gene_name and Phenotype_Hits
impc <- as.data.frame(impc_raw[, c("gene_name", "Phenotype_Hits")])
cat("  Processed IMPC data:", nrow(impc), "genes\n")

save(impc, file = "data/impc.rda", compress = "xz")
cat("  Saved to data/impc.rda\n\n")

# ==============================================================================
# 2. OpenTargets Platform - Genetic Associations
# ==============================================================================
cat("=== Processing OpenTargets genetic association data ===\n")

assoc_path <- file.path(DOWNLOAD_DIR, "associationByDatasourceDirect.json")
disease_path <- file.path(DOWNLOAD_DIR, "diseases.json")
gencode_gtf_path <- Sys.glob(file.path(DOWNLOAD_DIR, "gencode.v*.annotation.gtf"))

if (!file.exists(assoc_path)) {
  stop("Association file not found: ", assoc_path, "\n",
       "Run download_data.sh or download manually from OpenTargets Platform.")
}
if (!file.exists(disease_path)) {
  stop("Disease file not found: ", disease_path, "\n",
       "Run download_data.sh or download manually from OpenTargets Platform.")
}
if (length(gencode_gtf_path) == 0) {
  stop("GENCODE GTF file not found in: ", DOWNLOAD_DIR, "\n",
       "Download from: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/")
}
gencode_gtf_path <- gencode_gtf_path[1]

cat("  Reading association evidence data (this may take several minutes)...\n")
assoc_data <- stream_in(file(assoc_path))

cat("  Reading disease annotation data...\n")
disease_data <- stream_in(file(disease_path))
# Keep disease id and name columns
disease_data <- disease_data[, c("id", "name")]

cat("  Joining evidence with disease annotations...\n")
merged <- inner_join(disease_data, assoc_data, by = c("id" = "diseaseId"))

cat("  Reading GENCODE GTF (this may take a few minutes)...\n")
gencode <- as.data.frame(rtracklayer::import(gencode_gtf_path))
gencode_pc <- filter(gencode, gene_type == "protein_coding", type == "gene")
gencode_pc$gene_id <- gsub("\\..*", "", gencode_pc$gene_id)
# Keep key columns: gene_id, gene_name, and genomic coordinates
gencode_keep <- gencode_pc[, c("seqnames", "start", "end", "strand",
                                "gene_id", "gene_name")]
colnames(gencode_keep)[1] <- "chr"

cat("  Filtering for genetic associations and intersecting with protein-coding genes...\n")
assoc_filtered <- filter(merged, datatypeId == "genetic_association")
disease_target_genetic_association <- inner_join(
  assoc_filtered, gencode_keep,
  by = c("targetId" = "gene_id")
)

# Keep the most useful columns for the package
keep_cols <- intersect(
  c("name", "id", "targetId", "gene_name", "score",
    "datatypeId", "datasourceId", "chr", "start", "end", "strand"),
  colnames(disease_target_genetic_association)
)
disease_target_genetic_association <- disease_target_genetic_association[, keep_cols]

# Rename for clarity
if ("name" %in% colnames(disease_target_genetic_association)) {
  colnames(disease_target_genetic_association)[colnames(disease_target_genetic_association) == "name"] <- "disease_name"
}
if ("id" %in% colnames(disease_target_genetic_association)) {
  colnames(disease_target_genetic_association)[colnames(disease_target_genetic_association) == "id"] <- "disease_id"
}
if ("targetId" %in% colnames(disease_target_genetic_association)) {
  colnames(disease_target_genetic_association)[colnames(disease_target_genetic_association) == "targetId"] <- "gene_id"
}

disease_target_genetic_association <- as.data.frame(disease_target_genetic_association)
cat("  Processed association data:", nrow(disease_target_genetic_association), "rows,",
    length(unique(disease_target_genetic_association$gene_name)), "unique genes\n")

save(disease_target_genetic_association,
     file = "data/disease_target_genetic_association.rda", compress = "xz")
cat("  Saved to data/disease_target_genetic_association.rda\n\n")

# ==============================================================================
# 3. OpenTargets Genetics - Locus-to-Gene (L2G)
# ==============================================================================
cat("=== Processing OpenTargets Genetics L2G data ===\n")

l2g_dir <- file.path(DOWNLOAD_DIR, "l2g")
study_index_path <- file.path(DOWNLOAD_DIR, "study_index.json")

if (!dir.exists(l2g_dir)) {
  stop("L2G directory not found: ", l2g_dir, "\n",
       "Run download_data.sh or download the L2G parquet files manually.")
}
if (!file.exists(study_index_path)) {
  stop("Study index file not found: ", study_index_path, "\n",
       "Run download_data.sh or download manually from OpenTargets Genetics.")
}

cat("  Reading L2G parquet data (this may take a few minutes)...\n")
DS <- arrow::open_dataset(sources = l2g_dir)
SO <- Scanner$create(DS)
AT <- SO$ToTable()
l2g_data <- as.data.frame(AT)
cat("  L2G raw data:", nrow(l2g_data), "rows\n")

cat("  Reading study index...\n")
study_index <- stream_in(file(study_index_path))
# Keep key study columns
study_cols <- intersect(
  c("study_id", "trait_reported", "trait_efos", "trait_category",
    "pub_id", "pub_author", "ancestry_initial", "n_total"),
  colnames(study_index)
)
study_index <- study_index[, study_cols]

cat("  Merging L2G with study annotations...\n")
l2g_annotated <- left_join(l2g_data, study_index, by = "study_id")
l2g_annotated <- do.call(data.frame, l2g_annotated)

cat("  Intersecting with protein-coding genes from GENCODE...\n")
# Reuse gencode_keep from section 2
l2g_annotated_full <- inner_join(l2g_annotated, gencode_keep, by = "gene_id")

l2g_annotated_full <- as.data.frame(l2g_annotated_full)
cat("  Processed L2G data:", nrow(l2g_annotated_full), "rows,",
    length(unique(l2g_annotated_full$gene_name)), "unique genes\n")

save(l2g_annotated_full, file = "data/l2g_annotated_full.rda", compress = "xz")
cat("  Saved to data/l2g_annotated_full.rda\n\n")

# ==============================================================================
# Summary
# ==============================================================================
cat("============================================\n")
cat("  Data preparation complete!\n")
cat("============================================\n")
cat("\n")
cat("Saved files:\n")
cat("  data/impc.rda                                -", format(file.size("data/impc.rda"), big.mark = ","), "bytes\n")
cat("  data/disease_target_genetic_association.rda   -", format(file.size("data/disease_target_genetic_association.rda"), big.mark = ","), "bytes\n")
cat("  data/l2g_annotated_full.rda                   -", format(file.size("data/l2g_annotated_full.rda"), big.mark = ","), "bytes\n")
cat("\n")
cat("Dataset sizes:\n")
cat("  impc:                                  ", nrow(impc), "rows\n")
cat("  disease_target_genetic_association:     ", nrow(disease_target_genetic_association), "rows\n")
cat("  l2g_annotated_full:                    ", nrow(l2g_annotated_full), "rows\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Rebuild the package: R CMD build .\n")
cat("  2. Or use devtools::install() to install from source\n")
cat("  3. The new data will be available via data() or LazyData loading\n")
cat("\n")
cat("NOTE: The .rda files may be large. If the package exceeds CRAN's 5MB limit,\n")
cat("      consider hosting the data externally and downloading on first use.\n")
