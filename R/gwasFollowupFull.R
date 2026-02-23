#' gwasFollowupFull
#'
#' Extended version of gwasFollowup that also supports PhenomeXcan and TWAS data.
#' Takes GWAS summary statistics (with a closest gene column), a GTF file,
#' and a p-value threshold. For each significant locus it:
#' 1. Finds the closest gene (from the 'gene' column in sumstats)
#' 2. Extracts all protein-coding genes within 500kb of that gene's TSS from the GTF
#' 3. Translates genes via Zoonomia orthology (non-human species)
#' 4. Queries OpenTargets API for disease associations
#' 5. Queries IMPC API for mouse phenotypes
#' 6. Optionally joins PhenomeXcan and TWAS data
#'
#' @param sumStats GWAS summary statistics file. Must have columns: chr, ps, p_wald, gene
#' @param felGTF GTF file for the species
#' @param species Species of the GWAS data: "human", "cat", or "dog". Default is "cat".
#' @param pval P-value threshold. Default is 5e-8.
#' @param ResultsPath Directory to write output files. Created if it doesn't exist.
#' @param phenomePath path to phenomeXcan data if available default NULL
#' @param twasPath path to twas data if available default NULL
#' @param zoo_dir Directory containing Zoonomia orthology files. If NULL, uses package extdata.
#' @return Writes g2d_results.txt and impc_results.txt to ResultsPath
#' @importFrom dplyr filter left_join bind_rows select_if if_all
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @import  GenomicRanges
#' @importFrom rtracklayer import
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom tibble rownames_to_column
#' @export gwasFollowupFull

gwasFollowupFull <- function(sumStats, felGTF, species = "cat", pval = 0.00000005, ResultsPath = ".",
                             phenomePath = NULL, twasPath = NULL, zoo_dir = NULL) {
  # Validate input files
  if (!nzchar(sumStats) || !file.exists(sumStats)) {
    stop("Summary statistics file not found: '", sumStats, "'. ",
         "Check the file path or verify the package is installed correctly.")
  }
  if (!nzchar(felGTF) || !file.exists(felGTF)) {
    stop("GTF file not found: '", felGTF, "'. ",
         "Check the file path or verify the package is installed correctly.")
  }
  # Create output directory if it doesn't exist
  if (!dir.exists(ResultsPath)) {
    dir.create(ResultsPath, recursive = TRUE)
  }

  # Optional extra data
  if (!is.null(phenomePath)) {
    print("Reading phenomeXcan data")
    phenomeXcan <- fread(phenomePath)
  } else {
    phenomeXcan <- NULL
    print("No phenomeXcan data has been provided")
  }
  if (!is.null(twasPath)) {
    print("Reading TWAS data")
    twas <- fread(twasPath)
    twas$`ENSG ID` <- gsub("\\..*", "", twas$`ENSG ID`)
    twas <- twas[, c(2:15)]
  } else {
    twas <- NULL
    print("No TWAS data has been provided")
  }

  # --- Step 1: Read sumstats and filter by p-value ---
  print("Reading the summary statistics file, filtering based on the input p_value")
  gwas <- fread(sumStats)
  if (!"gene" %in% colnames(gwas)) {
    stop("Summary statistics file must contain a 'gene' column with the closest gene per locus.")
  }
  gwas$chr <- gsub("^20$", "X", gwas$chr)
  sig_gwas <- filter(gwas, p_wald <= pval)
  if (nrow(sig_gwas) == 0) {
    stop("No significant SNPs found at p-value threshold ", pval, ". Consider increasing it.")
  }

  # Get unique closest genes from significant loci
  closest_genes <- unique(sig_gwas$gene)
  print(paste0("Found ", length(closest_genes), " unique closest gene(s): ",
               paste(closest_genes, collapse = ", ")))

  # --- Step 2: Read GTF and find genes within 500kb of each closest gene's TSS ---
  print("Reading the GTF file and filtering for protein coding genes")
  cat_gtf <- as.data.frame(rtracklayer::import(felGTF))
  cat_pc_gene_gtf <- filter(cat_gtf, gene_biotype == "protein_coding", type == "gene")

  print("Finding genes within 500kb of each closest gene's TSS")
  all_nearby_genes <- character()
  for (cg in closest_genes) {
    gene_row <- cat_pc_gene_gtf[toupper(cat_pc_gene_gtf$gene_name) == toupper(cg), , drop = FALSE]
    if (nrow(gene_row) == 0) {
      warning("Closest gene '", cg, "' not found in GTF. Skipping.")
      next
    }
    tss <- ifelse(gene_row$strand[1] == "-", gene_row$end[1], gene_row$start[1])
    chr <- as.character(gene_row$seqnames[1])
    window_start <- tss - 500000
    window_end <- tss + 500000
    nearby <- cat_pc_gene_gtf[as.character(cat_pc_gene_gtf$seqnames) == chr &
                                cat_pc_gene_gtf$start <= window_end &
                                cat_pc_gene_gtf$end >= window_start, ]
    all_nearby_genes <- c(all_nearby_genes, nearby$gene_name)
  }
  genes <- unique(all_nearby_genes)

  if (length(genes) == 0) {
    stop("No protein-coding genes found within 500kb of any closest gene TSS.")
  }
  print(paste0("Found ", length(genes), " unique gene(s) within 500kb windows: ",
               paste(genes, collapse = ", ")))

  # --- Step 3: Zoonomia orthology translation ---
  human_ortho <- NULL
  mouse_ortho <- NULL
  if (species != "human") {
    print("Translating genes to human orthologs using Zoonomia data for OpenTargets")
    human_ortho <- translate_genes(genes, "human", species, zoo_dir)
    if (nrow(human_ortho) > 0) {
      ot_genes <- unique(human_ortho$target_gene)
    } else {
      warning("No human orthologs found in Zoonomia data. Using original gene names.")
      ot_genes <- genes
    }
    print("Translating genes to mouse orthologs using Zoonomia data for IMPC")
    mouse_ortho <- translate_genes(genes, "mouse", species, zoo_dir)
    if (nrow(mouse_ortho) > 0) {
      impc_genes <- unique(mouse_ortho$target_gene)
    } else {
      warning("No mouse orthologs found in Zoonomia data.")
      impc_genes <- character(0)
    }
  } else {
    ot_genes <- genes
    impc_genes <- genes
  }

  # --- Step 4: Query OpenTargets API ---
  print("Querying OpenTargets API for gene-disease associations")
  results_g2d <- list()
  pb <- txtProgressBar(min = 0, max = length(ot_genes), style = 3, width = 50, char = "=")
  for (i in seq_along(ot_genes)) {
    results_g2d[[i]] <- fetch_opentargets(ot_genes[i])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_g2d_df <- do.call("rbind", results_g2d)
  if (is.null(results_g2d_df)) results_g2d_df <- data.frame()

  # --- Step 5: Query IMPC API ---
  print("Querying IMPC API for mouse phenotype data")
  impc_results <- list()
  if (length(impc_genes) > 0) {
    pb <- txtProgressBar(min = 0, max = length(impc_genes), style = 3, width = 50, char = "=")
    for (i in seq_along(impc_genes)) {
      impc_results[[i]] <- fetch_impc(impc_genes[i])
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  impc_df <- do.call("rbind", impc_results)
  if (is.null(impc_df)) impc_df <- data.frame()

  # --- Step 6: Join results ---
  if (nrow(results_g2d_df) > 0 && nrow(impc_df) > 0) {
    results_g2d_df <- left_join(results_g2d_df, impc_df, by = "gene_name")
  }

  # PhenomeXcan integration
  if (!is.null(phenomeXcan) && nrow(results_g2d_df) > 0) {
    print("Wrangling phenomeXcan data")
    pheList <- list()
    for (i in seq_along(ot_genes)) {
      pheList[[i]] <- phenomeXcan[grep(ot_genes[i], phenomeXcan$gene_name), ]
      pheList[[i]] <- pheList[[i]] %>% select_if(~ !any(is.na(.)))
    }
    columns <- list()
    for (g in seq_along(ot_genes)) {
      columns[[g]] <- pheList[[g]]$gene_name
    }
    for (j in seq_along(ot_genes)) {
      pheList[[j]] <- subset(pheList[[j]], select = -gene_name)
      pheList[[j]] <- as.data.frame(t(pheList[[j]]))
      colnames(pheList[[j]]) <- columns[[j]]
      pheList[[j]] <- pheList[[j]] %>% filter(if_all(.fns = ~. >= 0.101))
      pheList[[j]] <- as.data.frame(t(pheList[[j]]))
      pheList[[j]] <- tibble::rownames_to_column(pheList[[j]], "gene_name")
    }
    phenomeXcan <- do.call("bind_rows", pheList)
    results_g2d_df <- left_join(results_g2d_df, phenomeXcan, by = "gene_name")
  }

  # TWAS integration
  if (!is.null(twas) && nrow(results_g2d_df) > 0) {
    results_g2d_df <- left_join(results_g2d_df, twas, by = c("gene_name" = "Gene"))
  }

  # Add Zoonomia orthology metadata and clean up
  if (nrow(results_g2d_df) > 0) {
    results_g2d_df <- results_g2d_df[, colSums(is.na(results_g2d_df)) < nrow(results_g2d_df), drop = FALSE]
    results_g2d_df <- add_orthology_info(results_g2d_df, human_ortho)
  }

  # --- Step 7: Write results ---
  print("Writing results")
  write.table(results_g2d_df, file.path(ResultsPath, "g2d_results.txt"),
              quote = FALSE, row.names = FALSE, sep = "\t")
  if (nrow(impc_df) > 0) {
    write.table(impc_df, file.path(ResultsPath, "impc_results.txt"),
                quote = FALSE, row.names = FALSE, sep = "\t")
  }

  print("All done, happy exploring")
}
