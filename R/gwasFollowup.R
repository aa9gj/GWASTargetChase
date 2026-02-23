#' gwasFollowup
#'
#' Main package function. Takes GWAS summary statistics (with a closest gene
#' column), a GTF file, and a p-value threshold. For each significant locus it:
#' 1. Finds the closest gene (from the 'gene' column in sumstats)
#' 2. Extracts all protein-coding genes within 500kb of that gene's TSS from the GTF
#' 3. Translates genes via Zoonomia orthology (non-human species)
#' 4. Queries OpenTargets API for disease associations
#' 5. Queries IMPC API for mouse phenotypes
#'
#' @param sumStats GWAS summary statistics file. Must have columns: chr, ps, p_wald, gene
#' @param felGTF GTF file for the species
#' @param species Species of the GWAS data: "human", "cat", or "dog". Default is "cat".
#' @param pval P-value threshold. Default is 5e-8.
#' @param ResultsPath Directory to write output files. Created if it doesn't exist.
#' @param zoo_dir Directory containing Zoonomia orthology files. If NULL, uses package extdata.
#' @return Writes g2d_results.txt, l2g_results.txt, and plots.pdf to ResultsPath
#' @importFrom dplyr filter left_join inner_join %>%
#' @importFrom data.table fread
#' @import GenomicRanges
#' @import ggplot2
#' @import patchwork
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices pdf dev.off
#' @importFrom rtracklayer import
#' @importFrom S4Vectors queryHits subjectHits
#' @export gwasFollowup

gwasFollowup <- function(sumStats, felGTF, species = "cat", pval = 0.00000005, ResultsPath = ".", zoo_dir = NULL) {
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
  sig_gwas$locus <- paste0("locus_", seq_len(nrow(sig_gwas)))

  # Get unique closest genes from significant loci
  closest_genes <- unique(sig_gwas$gene)
  print(paste0("Found ", length(closest_genes), " unique closest gene(s) across significant loci: ",
               paste(closest_genes, collapse = ", ")))

  # --- Step 2: Read GTF and find genes within 500kb of each closest gene's TSS ---
  print("Reading the GTF file and filtering for protein coding genes")
  cat_gtf <- as.data.frame(rtracklayer::import(felGTF))
  cat_pc_gene_gtf <- filter(cat_gtf, gene_biotype == "protein_coding", type == "gene")

  print("Finding genes within 500kb of each closest gene's TSS")
  all_nearby_genes <- character()
  for (cg in closest_genes) {
    # Find this gene in the GTF
    gene_row <- cat_pc_gene_gtf[toupper(cat_pc_gene_gtf$gene_name) == toupper(cg), , drop = FALSE]
    if (nrow(gene_row) == 0) {
      warning("Closest gene '", cg, "' not found in GTF. Skipping.")
      next
    }
    # Use TSS (start for + strand, end for - strand)
    tss <- ifelse(gene_row$strand[1] == "-", gene_row$end[1], gene_row$start[1])
    chr <- as.character(gene_row$seqnames[1])

    # Find all protein-coding genes within 500kb of TSS
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

  # --- Step 3: Zoonomia orthology translation for non-human species ---
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

  # --- Step 4: Query OpenTargets API for disease associations ---
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

  # --- Step 5: Query IMPC API for mouse phenotypes ---
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

  # --- Step 6: Join IMPC phenotypes with disease results ---
  if (nrow(results_g2d_df) > 0 && nrow(impc_df) > 0) {
    results_g2d_df <- left_join(results_g2d_df, impc_df, by = "gene_name")
  }

  # Add Zoonomia orthology metadata
  if (nrow(results_g2d_df) > 0) {
    results_g2d_df <- add_orthology_info(results_g2d_df, human_ortho)
  }

  # --- Step 7: Write results ---
  print("Writing results")
  write.table(results_g2d_df, file.path(ResultsPath, "g2d_results.txt"),
              quote = FALSE, row.names = FALSE, sep = "\t")

  # Write a summary of IMPC results
  if (nrow(impc_df) > 0) {
    write.table(impc_df, file.path(ResultsPath, "impc_results.txt"),
                quote = FALSE, row.names = FALSE, sep = "\t")
  }

  # --- Step 8: Create plots ---
  print("Creating plots per significant locus")
  sig_gwas$start <- sig_gwas$ps - 500000
  sig_gwas$end <- sig_gwas$ps + 500000
  sig_gwas_ranges <- makeGRangesFromDataFrame(sig_gwas, keep.extra.columns = TRUE, seqnames.field = "chr")
  cat_pc_gene_ranges <- makeGRangesFromDataFrame(cat_pc_gene_gtf, keep.extra.columns = TRUE)

  # SNP to gene mapping for plots
  hits <- findOverlaps(sig_gwas_ranges, cat_pc_gene_ranges)
  if (length(hits) > 0) {
    olap_snp <- as.data.frame(sig_gwas_ranges[queryHits(hits)])
    olap_gene <- as.data.frame(cat_pc_gene_ranges[subjectHits(hits)])
    sig_SNP_info <- data.frame(
      rs = olap_snp$rs, ps = olap_snp$ps, p_wald = olap_snp$p_wald,
      locus = olap_snp$locus, gene_name = olap_gene$gene_name,
      gene_start = olap_gene$start, gene_end = olap_gene$end,
      stringsAsFactors = FALSE
    )

    # Full GWAS data for manhattan-style plots
    gwas_sum_ranges <- makeGRangesFromDataFrame(gwas, seqnames.field = "chr",
                                                 start.field = "ps", end.field = "ps",
                                                 keep.extra.columns = TRUE)
    hits4 <- findOverlaps(sig_gwas_ranges, gwas_sum_ranges)

    plot.list <- list()
    unique_loci <- unique(sig_SNP_info$locus)
    for (i in seq_along(unique_loci)) {
      locus_snps <- filter(sig_SNP_info, locus == unique_loci[i])
      # Get all GWAS SNPs in this locus window
      locus_idx <- which(sig_gwas$locus == unique_loci[i])
      if (length(locus_idx) == 0) next
      locus_range <- sig_gwas_ranges[locus_idx[1]]
      locus_hits <- findOverlaps(locus_range, gwas_sum_ranges)
      if (length(locus_hits) == 0) next
      locus_gwas <- as.data.frame(gwas_sum_ranges[subjectHits(locus_hits)])

      p3 <- ggplot(data = locus_gwas) +
        geom_point(aes(start, -log10(p_wald))) +
        theme(text = element_text(size = 6)) +
        ylab("-Log10 P-value") +
        xlab("") +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

      p4 <- ggplot(data = locus_snps) +
        geom_linerange(aes(x = gene_name, ymin = gene_start, ymax = gene_end)) +
        coord_flip() +
        xlab("Gene Name") + theme(text = element_text(size = 6))

      plot.list[[i]] <- p3 + p4 + plot_layout(ncol = 1, heights = c(6, 6))
    }

    if (length(plot.list) > 0) {
      pdf(file.path(ResultsPath, "plots.pdf"))
      for (i in seq_along(plot.list)) {
        if (!is.null(plot.list[[i]])) print(plot.list[[i]])
      }
      dev.off()
    }
  }

  print("All done, happy exploring")
}
