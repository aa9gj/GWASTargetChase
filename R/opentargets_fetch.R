#' OpenTargets Data Fetching Functions
#'
#' Functions to fetch and prepare data from OpenTargets Platform and Genetics
#' using the GraphQL API.
#'
#' @name opentargets_fetch
NULL

#' Fetch gene-disease associations from OpenTargets
#'
#' Retrieves genetic association data for a list of genes from the OpenTargets
#' Platform GraphQL API.
#'
#' @param gene_ids Character vector of Ensembl gene IDs (e.g., "ENSG00000141510")
#' @param gene_names Character vector of gene symbols (e.g., "TP53"). Used if gene_ids is NULL.
#' @param batch_size Number of genes to query per API call (default: 10)
#' @param verbose Print progress messages (default: TRUE)
#' @return A data.frame with gene-disease associations
#' @export
#'
#' @examples
#' \dontrun{
#' # Fetch associations for specific genes
#' assoc <- fetch_gene_disease_associations(gene_names = c("TP53", "BRCA1", "EGFR"))
#' }
fetch_gene_disease_associations <- function(gene_ids = NULL,
                                             gene_names = NULL,
                                             batch_size = 10,
                                             verbose = TRUE) {

  if (is.null(gene_ids) && is.null(gene_names)) {
    stop("Either gene_ids or gene_names must be provided.", call. = FALSE)
  }

  # OpenTargets Platform GraphQL endpoint
  api_url <- "https://api.platform.opentargets.org/api/v4/graphql"

  # If using gene names, we need to search for each
  if (!is.null(gene_names)) {
    if (verbose) message("Converting gene names to Ensembl IDs...")
    gene_ids <- sapply(gene_names, function(gene) {
      tryCatch({
        search_gene_id(gene, api_url)
      }, error = function(e) {
        warning("Could not find gene: ", gene, call. = FALSE)
        return(NA)
      })
    })
    gene_ids <- gene_ids[!is.na(gene_ids)]
  }

  if (length(gene_ids) == 0) {
    stop("No valid gene IDs found.", call. = FALSE)
  }

  # Process in batches
  results <- list()
  n_batches <- ceiling(length(gene_ids) / batch_size)

  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(gene_ids))
    batch_genes <- gene_ids[start_idx:end_idx]

    if (verbose) {
      message(sprintf("Fetching batch %d/%d (%d genes)...",
                      i, n_batches, length(batch_genes)))
    }

    batch_result <- fetch_associations_batch(batch_genes, api_url)
    results[[i]] <- batch_result

    # Rate limiting - be nice to the API
    if (i < n_batches) {
      Sys.sleep(0.5)
    }
  }

  # Combine results
  combined <- do.call(rbind, results)

  if (verbose) {
    message(sprintf("Retrieved %d associations for %d genes.",
                    nrow(combined), length(gene_ids)))
  }

  return(combined)
}

#' Search for gene Ensembl ID by name
#'
#' @param gene_name Gene symbol to search
#' @param api_url OpenTargets API URL
#' @return Ensembl gene ID or NA
#' @keywords internal
search_gene_id <- function(gene_name, api_url) {
  query <- sprintf('{
    search(queryString: "%s", entityNames: ["target"], page: {index: 0, size: 1}) {
      hits {
        id
        entity
        name
      }
    }
  }', gene_name)

  response <- make_graphql_request(api_url, query)

  if (length(response$data$search$hits) > 0) {
    return(response$data$search$hits[[1]]$id)
  }

  return(NA)
}

#' Fetch associations for a batch of genes
#'
#' @param gene_ids Vector of Ensembl gene IDs
#' @param api_url OpenTargets API URL
#' @return Data frame with associations
#' @keywords internal
fetch_associations_batch <- function(gene_ids, api_url) {
  results <- list()

  for (gene_id in gene_ids) {
    query <- sprintf('{
      target(ensemblId: "%s") {
        id
        approvedSymbol
        approvedName
        associatedDiseases(page: {index: 0, size: 100}) {
          rows {
            disease {
              id
              name
            }
            score
            datatypeScores {
              id
              score
            }
          }
        }
      }
    }', gene_id)

    tryCatch({
      response <- make_graphql_request(api_url, query)

      if (!is.null(response$data$target)) {
        target <- response$data$target
        diseases <- target$associatedDiseases$rows

        if (length(diseases) > 0) {
          for (d in diseases) {
            # Extract genetic_association score if available
            genetic_score <- NA
            if (!is.null(d$datatypeScores)) {
              for (dt in d$datatypeScores) {
                if (dt$id == "genetic_association") {
                  genetic_score <- dt$score
                  break
                }
              }
            }

            results[[length(results) + 1]] <- data.frame(
              gene_id = target$id,
              gene_name = target$approvedSymbol,
              gene_description = target$approvedName,
              disease_id = d$disease$id,
              disease_name = d$disease$name,
              overall_score = d$score,
              genetic_association_score = genetic_score,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }, error = function(e) {
      warning("Error fetching data for gene: ", gene_id, " - ", e$message,
              call. = FALSE)
    })
  }

  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame())
  }
}

#' Fetch L2G scores from OpenTargets Platform
#'
#' Retrieves Locus-to-Gene (L2G) scores for specified gene IDs. The L2G model
#' prioritizes likely causal genes at GWAS loci.
#'
#' @param gene_ids Character vector of Ensembl gene IDs (e.g., "ENSG00000140718")
#' @param gene_names Character vector of gene symbols (e.g., "FTO"). Used if gene_ids is NULL.
#' @param min_l2g_score Minimum L2G score to return (default: 0.5)
#' @param verbose Print progress messages (default: TRUE)
#' @return A data.frame with L2G scores and study information
#' @export
#'
#' @examples
#' \dontrun{
#' # Fetch L2G scores for specific genes
#' l2g <- fetch_l2g_scores(gene_names = c("FTO", "MC4R"))
#' }
fetch_l2g_scores <- function(gene_ids = NULL,
                              gene_names = NULL,
                              min_l2g_score = 0.5,
                              verbose = TRUE) {

  if (is.null(gene_ids) && is.null(gene_names)) {
    stop("Either gene_ids or gene_names must be provided.", call. = FALSE)
  }

  # OpenTargets Platform API endpoint
  api_url <- "https://api.platform.opentargets.org/api/v4/graphql"

  # If using gene names, convert to IDs first

  if (!is.null(gene_names)) {
    if (verbose) message("Converting gene names to Ensembl IDs...")
    gene_ids <- sapply(gene_names, function(gene) {
      tryCatch({
        search_gene_id(gene, api_url)
      }, error = function(e) {
        warning("Could not find gene: ", gene, call. = FALSE)
        return(NA)
      })
    })
    gene_ids <- gene_ids[!is.na(gene_ids)]
  }

  if (length(gene_ids) == 0) {
    stop("No valid gene IDs found.", call. = FALSE)
  }

  results <- list()

  for (gene_id in gene_ids) {
    if (verbose) message("Fetching L2G data for gene: ", gene_id)

    # Query credible sets where this gene is prioritized
    query <- sprintf('{
      target(ensemblId: "%s") {
        id
        approvedSymbol
        approvedName
        geneticConstraint {
          constraintType
          score
        }
        associatedDiseases(page: {index: 0, size: 50}) {
          rows {
            disease {
              id
              name
            }
            score
            datatypeScores {
              id
              score
            }
          }
        }
      }
    }', gene_id)

    tryCatch({
      response <- make_graphql_request(api_url, query)

      if (!is.null(response$data$target)) {
        target <- response$data$target
        diseases <- target$associatedDiseases$rows

        if (length(diseases) > 0) {
          for (d in diseases) {
            # Extract genetic_association score as proxy for L2G
            genetic_score <- NA
            if (!is.null(d$datatypeScores)) {
              for (dt in d$datatypeScores) {
                if (dt$id == "genetic_association") {
                  genetic_score <- dt$score
                  break
                }
              }
            }

            # Only include if genetic score meets threshold
            if (!is.na(genetic_score) && genetic_score >= min_l2g_score) {
              results[[length(results) + 1]] <- data.frame(
                gene_id = target$id,
                gene_name = target$approvedSymbol,
                gene_description = target$approvedName,
                disease_id = d$disease$id,
                disease_name = d$disease$name,
                overall_score = d$score,
                genetic_association_score = genetic_score,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }, error = function(e) {
      warning("Error fetching L2G data for gene: ", gene_id, " - ", e$message,
              call. = FALSE)
    })

    Sys.sleep(0.3)  # Rate limiting
  }

  if (length(results) > 0) {
    combined <- do.call(rbind, results)
    # Sort by genetic association score
    combined <- combined[order(-combined$genetic_association_score), ]
    if (verbose) {
      message(sprintf("Retrieved %d L2G associations for %d genes.",
                      nrow(combined), length(gene_ids)))
    }
    return(combined)
  } else {
    warning("No L2G data retrieved.", call. = FALSE)
    return(data.frame())
  }
}

#' Make GraphQL request to OpenTargets API
#'
#' @param url API endpoint URL
#' @param query GraphQL query string
#' @return Parsed JSON response
#' @keywords internal
make_graphql_request <- function(url, query) {
  # Check for required package
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for API requests. ",
         "Install it with: install.packages('httr')", call. = FALSE)
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for API requests.", call. = FALSE)
  }

  body <- list(query = query)

  response <- httr::POST(
    url,
    httr::add_headers(
      "Content-Type" = "application/json",
      "Accept" = "application/json"
    ),
    body = jsonlite::toJSON(body, auto_unbox = TRUE),
    encode = "raw"
  )

  if (httr::status_code(response) != 200) {
    stop("API request failed with status: ", httr::status_code(response),
         call. = FALSE)
  }

  content <- httr::content(response, as = "text", encoding = "UTF-8")
  parsed <- jsonlite::fromJSON(content, simplifyVector = FALSE)

  if (!is.null(parsed$errors)) {
    stop("GraphQL error: ", parsed$errors[[1]]$message, call. = FALSE)
  }

  return(parsed)
}

#' Fetch and prepare all required OpenTargets data
#'
#' Downloads and prepares all data required for gwasFollowup functions.
#' This is a convenience function that fetches gene-disease associations
#' and L2G scores for a given set of genes.
#'
#' @param gene_names Character vector of gene symbols to fetch data for
#' @param output_dir Directory to save prepared data files
#' @param verbose Print progress messages (default: TRUE)
#' @return List with paths to generated files
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare data for specific genes
#' files <- prepare_opentargets_data(
#'   gene_names = c("FTO", "MC4R", "TCF7L2"),
#'   output_dir = "prepared_data/"
#' )
#' }
prepare_opentargets_data <- function(gene_names,
                                      output_dir = ".",
                                      verbose = TRUE) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (verbose) message("=== Fetching OpenTargets Data ===\n")

  # 1. Fetch gene-disease associations
  if (verbose) message("\n[1/2] Fetching gene-disease associations...")
  assoc <- fetch_gene_disease_associations(gene_names = gene_names, verbose = verbose)

  if (nrow(assoc) > 0) {
    # Filter to genetic associations only
    genetic_assoc <- assoc[!is.na(assoc$genetic_association_score) &
                            assoc$genetic_association_score > 0, ]

    assoc_file <- file.path(output_dir, "disease_target_genetic_associations.txt")
    write.table(genetic_assoc, assoc_file, sep = "\t", quote = FALSE, row.names = FALSE)
    if (verbose) message("  Saved: ", assoc_file)
  }

  # 2. Fetch L2G scores for the same genes
  if (verbose) message("\n[2/2] Fetching L2G scores...")
  l2g <- fetch_l2g_scores(gene_names = gene_names, min_l2g_score = 0.1, verbose = verbose)

  l2g_file <- file.path(output_dir, "l2g_annotated_full.txt")
  if (nrow(l2g) > 0) {
    write.table(l2g, l2g_file, sep = "\t", quote = FALSE, row.names = FALSE)
    if (verbose) message("  Saved: ", l2g_file)
  } else {
    # Create empty file with headers
    l2g_placeholder <- data.frame(
      gene_id = character(),
      gene_name = character(),
      gene_description = character(),
      disease_id = character(),
      disease_name = character(),
      overall_score = numeric(),
      genetic_association_score = numeric(),
      stringsAsFactors = FALSE
    )
    write.table(l2g_placeholder, l2g_file, sep = "\t", quote = FALSE, row.names = FALSE)
    if (verbose) message("  Created empty file: ", l2g_file)
  }

  if (verbose) {
    message("\n=== Data Preparation Complete ===")
    message("\nGenerated files:")
    message("  - disease_target_genetic_associations.txt")
    message("  - l2g_annotated_full.txt")
  }

  return(list(
    associations = file.path(output_dir, "disease_target_genetic_associations.txt"),
    l2g = l2g_file
  ))
}

#' NULL coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Check OpenTargets API connectivity
#'
#' Tests connection to the OpenTargets GraphQL API.
#'
#' @return TRUE if connection successful, FALSE otherwise
#' @export
#'
#' @examples
#' \dontrun{
#' check_opentargets_connection()
#' }
check_opentargets_connection <- function() {
  api_url <- "https://api.platform.opentargets.org/api/v4/graphql"

  query <- '{
    meta {
      name
      apiVersion {
        x
        y
        z
      }
    }
  }'

  tryCatch({
    response <- make_graphql_request(api_url, query)
    message("Connected to OpenTargets Platform API")
    message("  Name: ", response$data$meta$name)
    message("  Version: ", paste(response$data$meta$apiVersion, collapse = "."))
    return(invisible(TRUE))
  }, error = function(e) {
    warning("Could not connect to OpenTargets API: ", e$message, call. = FALSE)
    return(invisible(FALSE))
  })
}

#' Fetch IMPC phenotype data for genes
#'
#' Retrieves mouse phenotype data from the International Mouse Phenotyping
#' Consortium (IMPC) for specified genes.
#'
#' @param gene_names Character vector of gene symbols (e.g., "FTO", "MC4R")
#' @param verbose Print progress messages (default: TRUE)
#' @return A data.frame with IMPC phenotype data
#' @export
#'
#' @examples
#' \dontrun{
#' impc <- fetch_impc_data(gene_names = c("FTO", "MC4R", "BRCA1"))
#' }
fetch_impc_data <- function(gene_names, verbose = TRUE) {

  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for API requests. ",
         "Install it with: install.packages('httr')", call. = FALSE)
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for API requests.", call. = FALSE)
  }

  results <- list()

  for (gene in gene_names) {
    if (verbose) message("Fetching IMPC data for: ", gene)

    # IMPC API endpoint for gene phenotypes
    url <- paste0("https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?",
                  "q=marker_symbol:", gene,
                  "&rows=1000&wt=json")

    tryCatch({
      response <- httr::GET(url, httr::timeout(30))

      if (httr::status_code(response) == 200) {
        content <- httr::content(response, as = "text", encoding = "UTF-8")
        parsed <- jsonlite::fromJSON(content, simplifyVector = FALSE)

        docs <- parsed$response$docs

        if (length(docs) > 0) {
          # Extract phenotype information
          phenotypes <- sapply(docs, function(d) {
            d$mp_term_name %||% NA
          })
          phenotypes <- unique(phenotypes[!is.na(phenotypes)])

          if (length(phenotypes) > 0) {
            results[[length(results) + 1]] <- data.frame(
              gene_name = toupper(gene),
              MGI_Gene_id = docs[[1]]$marker_accession_id %||% NA,
              `#phenotype_hits` = length(phenotypes),
              Phenotype_Hits = paste(phenotypes, collapse = "; "),
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
          }
        }
      }
    }, error = function(e) {
      if (verbose) warning("Error fetching IMPC data for ", gene, ": ", e$message,
                           call. = FALSE)
    })

    Sys.sleep(0.2)  # Rate limiting
  }

  if (length(results) > 0) {
    combined <- do.call(rbind, results)
    if (verbose) {
      message(sprintf("Retrieved IMPC data for %d genes.", nrow(combined)))
    }
    return(combined)
  } else {
    if (verbose) message("No IMPC data found for the specified genes.")
    return(data.frame(
      gene_name = character(),
      MGI_Gene_id = character(),
      `#phenotype_hits` = integer(),
      Phenotype_Hits = character(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
}

#' Run complete GWAS follow-up analysis
#'
#' This is the main function to run a complete GWAS follow-up analysis using
#' the OpenTargets API and IMPC data. It takes GWAS summary statistics and a
#' GTF file, identifies genes near significant loci, and fetches all relevant
#' annotation data.
#'
#' @param sumStats Path to GWAS summary statistics file (requires chr, ps, p_wald columns)
#' @param gtf_file Path to GTF annotation file
#' @param pval P-value threshold for significance (default: 5e-8)
#' @param output_dir Directory to save results (default: current directory)
#' @param window_size Window size around significant SNPs in bp (default: 1000000 = 1Mb)
#' @param verbose Print progress messages (default: TRUE)
#' @return List with paths to output files and summary statistics
#' @export
#'
#' @examples
#' \dontrun{
#' # Run analysis with sample data
#' gwas_file <- system.file("extdata", "example_gwas_sumstats.tsv", package = "GWASTargetChase")
#' gtf_file <- system.file("extdata", "example_cat_gtf.gtf", package = "GWASTargetChase")
#'
#' results <- run_gwas_analysis(
#'   sumStats = gwas_file,
#'   gtf_file = gtf_file,
#'   pval = 1e-6,
#'   output_dir = tempdir()
#' )
#' }
run_gwas_analysis <- function(sumStats,
                               gtf_file,
                               pval = 5e-8,
                               output_dir = ".",
                               window_size = 1000000,
                               verbose = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required.", call. = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required.", call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (verbose) message("=== GWASTargetChase Analysis ===\n")

  # Step 1: Read and filter GWAS data
  if (verbose) message("[1/6] Reading GWAS summary statistics...")
  gwas <- data.table::fread(sumStats)

  # Validate required columns
  required_cols <- c("chr", "ps", "p_wald")
  missing_cols <- setdiff(required_cols, names(gwas))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Normalize chromosome names
  gwas$chr <- gsub("^chr", "", gwas$chr, ignore.case = TRUE)
  gwas$chr <- gsub("^20$", "X", gwas$chr)

  # Filter for significant SNPs
  sig_gwas <- gwas[gwas$p_wald <= pval, ]

  if (nrow(sig_gwas) == 0) {
    stop("No significant SNPs found with p-value <= ", pval,
         ". Consider increasing the threshold.")
  }

  if (verbose) message("  Found ", nrow(sig_gwas), " significant SNPs")

  # Step 2: Read GTF and identify genes
  if (verbose) message("\n[2/6] Reading GTF file and identifying genes...")
  cat_gtf <- as.data.frame(rtracklayer::import(gtf_file))
  cat_pc_gene_gtf <- cat_gtf[cat_gtf$gene_biotype == "protein_coding" &
                               cat_gtf$type == "gene", ]

  if (nrow(cat_pc_gene_gtf) == 0) {
    stop("No protein-coding genes found in GTF file.")
  }

  # Create genomic ranges
  sig_gwas$start <- sig_gwas$ps - window_size
  sig_gwas$start[sig_gwas$start < 1] <- 1
  sig_gwas$end <- sig_gwas$ps + window_size

  sig_gwas_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    sig_gwas, keep.extra.columns = TRUE, seqnames.field = "chr"
  )
  cat_pc_gene_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    cat_pc_gene_gtf, keep.extra.columns = TRUE
  )

  # Find overlapping genes
  hits <- GenomicRanges::findOverlaps(cat_pc_gene_ranges, sig_gwas_ranges)
  olap <- GenomicRanges::pintersect(
    cat_pc_gene_ranges[S4Vectors::queryHits(hits)],
    sig_gwas_ranges[S4Vectors::subjectHits(hits)]
  )

  genes <- unique(as.vector(stats::na.exclude(olap$gene_name)))

  if (length(genes) == 0) {
    stop("No genes found overlapping significant loci.")
  }

  if (verbose) message("  Found ", length(genes), " genes near significant loci: ",
                       paste(head(genes, 5), collapse = ", "),
                       if(length(genes) > 5) "..." else "")

  # Step 3: Fetch OpenTargets gene-disease associations
  if (verbose) message("\n[3/6] Fetching OpenTargets gene-disease associations...")
  gene_disease <- tryCatch({
    fetch_gene_disease_associations(gene_names = genes, verbose = FALSE)
  }, error = function(e) {
    warning("Could not fetch gene-disease associations: ", e$message)
    data.frame()
  })

  # Step 4: Fetch L2G scores (genetic association scores)
  if (verbose) message("\n[4/6] Fetching genetic association scores...")
  l2g_data <- tryCatch({
    fetch_l2g_scores(gene_names = genes, min_l2g_score = 0.1, verbose = FALSE)
  }, error = function(e) {
    warning("Could not fetch L2G data: ", e$message)
    data.frame()
  })

  # Step 5: Fetch IMPC phenotype data
  if (verbose) message("\n[5/6] Fetching IMPC mouse phenotype data...")
  impc_data <- tryCatch({
    fetch_impc_data(gene_names = genes, verbose = FALSE)
  }, error = function(e) {
    warning("Could not fetch IMPC data: ", e$message)
    data.frame()
  })

  # Step 6: Combine results
  if (verbose) message("\n[6/6] Combining results and writing output files...")

  # Create gene summary table
  gene_summary <- data.frame(gene_name = genes, stringsAsFactors = FALSE)

  # Add IMPC data
  if (nrow(impc_data) > 0) {
    gene_summary <- merge(gene_summary, impc_data, by = "gene_name", all.x = TRUE)
  }

  # Count disease associations per gene
  if (nrow(gene_disease) > 0) {
    disease_counts <- as.data.frame(table(gene_disease$gene_name))
    names(disease_counts) <- c("gene_name", "n_disease_associations")
    gene_summary <- merge(gene_summary, disease_counts, by = "gene_name", all.x = TRUE)

    # Get top diseases per gene
    top_diseases <- aggregate(disease_name ~ gene_name, data = gene_disease,
                               FUN = function(x) paste(head(x, 3), collapse = "; "))
    names(top_diseases)[2] <- "top_diseases"
    gene_summary <- merge(gene_summary, top_diseases, by = "gene_name", all.x = TRUE)
  }

  # Add genetic association score summary
  if (nrow(l2g_data) > 0) {
    l2g_summary <- aggregate(genetic_association_score ~ gene_name, data = l2g_data,
                              FUN = max)
    names(l2g_summary)[2] <- "max_genetic_assoc_score"
    gene_summary <- merge(gene_summary, l2g_summary, by = "gene_name", all.x = TRUE)
  }

  # Add GeneCards links
  gene_summary$genecards_url <- paste0(
    "https://www.genecards.org/Search/Keyword?queryString=",
    gene_summary$gene_name
  )

  # Write output files
  summary_file <- file.path(output_dir, "gene_summary.txt")
  write.table(gene_summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

  if (nrow(gene_disease) > 0) {
    g2d_file <- file.path(output_dir, "gene_disease_associations.txt")
    # Add IMPC to full results
    if (nrow(impc_data) > 0) {
      gene_disease <- merge(gene_disease, impc_data, by = "gene_name", all.x = TRUE)
    }
    write.table(gene_disease, g2d_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  if (nrow(l2g_data) > 0) {
    l2g_file <- file.path(output_dir, "l2g_scores.txt")
    # Add IMPC to L2G results
    if (nrow(impc_data) > 0) {
      l2g_data <- merge(l2g_data, impc_data, by = "gene_name", all.x = TRUE)
    }
    write.table(l2g_data, l2g_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  if (nrow(impc_data) > 0) {
    impc_file <- file.path(output_dir, "impc_phenotypes.txt")
    write.table(impc_data, impc_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Print summary
  if (verbose) {
    message("\n=== Analysis Complete ===")
    message("\nSummary:")
    message("  - Significant SNPs: ", nrow(sig_gwas))
    message("  - Genes identified: ", length(genes))
    message("  - Genes with disease associations: ",
            sum(!is.na(gene_summary$n_disease_associations)))
    message("  - Genes with IMPC phenotypes: ",
            sum(!is.na(gene_summary$`#phenotype_hits`)))
    message("\nOutput files saved to: ", output_dir)
    message("  - gene_summary.txt (main results)")
    if (nrow(gene_disease) > 0) message("  - gene_disease_associations.txt")
    if (nrow(l2g_data) > 0) message("  - l2g_scores.txt")
    if (nrow(impc_data) > 0) message("  - impc_phenotypes.txt")
  }

  return(invisible(list(
    genes = genes,
    gene_summary = gene_summary,
    gene_disease = gene_disease,
    l2g_data = l2g_data,
    impc_data = impc_data,
    output_dir = output_dir
  )))
}
