#' API Fetch Functions for OpenTargets and IMPC
#'
#' Functions to query OpenTargets Platform GraphQL API and IMPC REST API
#' with local file caching.
#'
#' @name api_fetch
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom data.table fread
NULL

# Cache directory helper
.get_cache_dir <- function() {
  cache_dir <- file.path(tempdir(), "GWASTargetChase_cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  cache_dir
}

#' Fetch gene-disease associations from OpenTargets Platform API
#'
#' Queries the OpenTargets Platform GraphQL API for genetic association
#' evidence linked to a gene. Results are cached locally per session.
#'
#' @param gene_name Gene symbol (e.g. "FTO")
#' @param use_cache Use cached results if available (default TRUE)
#' @return data.frame with disease associations or empty data.frame
#' @export
fetch_opentargets <- function(gene_name, use_cache = TRUE) {
  cache_file <- file.path(.get_cache_dir(), paste0("ot_", gene_name, ".rds"))
  if (use_cache && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  api_url <- "https://api.platform.opentargets.org/api/v4/graphql"

  # First resolve gene name to Ensembl ID
  search_query <- paste0('{ search(queryString: "', gene_name,
                         '", entityNames: ["target"], page: { size: 1, index: 0 }) {',
                         ' hits { id name } } }')

  result <- tryCatch({
    response <- .post_graphql(api_url, search_query)
    hits <- response$data$search$hits
    if (is.null(hits) || length(hits) == 0) return(data.frame())
    ensembl_id <- hits$id[1]

    # Query for genetic associations
    assoc_query <- paste0(
      '{ target(ensemblId: "', ensembl_id, '") {',
      '  approvedSymbol',
      '  associatedDiseases(page: { size: 100, index: 0 }) {',
      '    rows {',
      '      disease { id name }',
      '      score',
      '      datasourceScores {',
      '        componentId: id',
      '        score',
      '      }',
      '    }',
      '  }',
      '} }'
    )

    response2 <- .post_graphql(api_url, assoc_query)
    target <- response2$data$target
    if (is.null(target)) return(data.frame())

    rows <- target$associatedDiseases$rows
    if (is.null(rows) || length(rows) == 0) return(data.frame())

    # Build result data.frame
    df <- data.frame(
      gene_name = target$approvedSymbol,
      gene_id = ensembl_id,
      disease_id = sapply(rows, function(r) r$disease$id),
      disease_name = sapply(rows, function(r) r$disease$name),
      overall_score = sapply(rows, function(r) r$score),
      stringsAsFactors = FALSE
    )

    # Extract genetic_association datasource score if available
    df$genetic_association_score <- sapply(rows, function(r) {
      ds <- r$datasourceScores
      if (is.null(ds) || length(ds) == 0) return(NA_real_)
      ga <- ds[sapply(ds, function(d) d$componentId) == "ot_genetics_portal", ]
      if (length(ga) == 0 || is.null(ga[[1]])) return(NA_real_)
      ga[[1]]$score
    })

    # Keep only rows with genetic association evidence
    df <- df[!is.na(df$genetic_association_score) & df$genetic_association_score > 0, ]

    if (use_cache && nrow(df) > 0) saveRDS(df, cache_file)
    df
  }, error = function(e) {
    warning("OpenTargets API query failed for ", gene_name, ": ", e$message)
    data.frame()
  })

  result
}

#' Fetch phenotype data from IMPC API
#'
#' Queries the IMPC REST API for mouse phenotype data linked to a gene.
#' Results are cached locally per session.
#'
#' @param gene_name Gene symbol (e.g. "FTO", "Fto")
#' @param use_cache Use cached results if available (default TRUE)
#' @return data.frame with phenotype hits or empty data.frame
#' @export
fetch_impc <- function(gene_name, use_cache = TRUE) {
  cache_file <- file.path(.get_cache_dir(), paste0("impc_", toupper(gene_name), ".rds"))
  if (use_cache && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  # IMPC gene search endpoint
  base_url <- "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select"

  result <- tryCatch({
    query_url <- paste0(base_url, "?q=marker_symbol:", gene_name,
                        "&rows=500&wt=json")
    response <- jsonlite::fromJSON(query_url, simplifyVector = TRUE)
    docs <- response$response$docs

    if (is.null(docs) || length(docs) == 0 || nrow(docs) == 0) {
      # Try uppercase
      query_url2 <- paste0(base_url, "?q=marker_symbol:", toupper(gene_name),
                           "&rows=500&wt=json")
      response2 <- jsonlite::fromJSON(query_url2, simplifyVector = TRUE)
      docs <- response2$response$docs
    }

    if (is.null(docs) || length(docs) == 0 || nrow(docs) == 0) {
      return(data.frame())
    }

    # Extract key columns
    keep_cols <- intersect(c("marker_symbol", "marker_accession_id",
                              "mp_term_name", "mp_term_id",
                              "top_level_mp_term_name",
                              "p_value", "effect_size",
                              "procedure_name", "parameter_name"),
                           colnames(docs))
    df <- docs[, keep_cols, drop = FALSE]
    colnames(df)[colnames(df) == "marker_symbol"] <- "gene_name"

    # Summarize: unique phenotypes per gene
    phenotypes <- unique(df$mp_term_name)
    summary_df <- data.frame(
      gene_name = toupper(gene_name),
      num_phenotype_hits = length(phenotypes),
      Phenotype_Hits = paste(phenotypes, collapse = "; "),
      stringsAsFactors = FALSE
    )

    if (use_cache) saveRDS(summary_df, cache_file)
    summary_df
  }, error = function(e) {
    warning("IMPC API query failed for ", gene_name, ": ", e$message)
    data.frame()
  })

  result
}

# Internal: POST a GraphQL query
.post_graphql <- function(url, query) {
  con <- url(url, method = "libcurl")
  on.exit(close(con))

  body <- jsonlite::toJSON(list(query = query), auto_unbox = TRUE)

  response <- tryCatch({
    tmp <- tempfile(fileext = ".json")
    on.exit(unlink(tmp), add = TRUE)
    writeBin(charToRaw(as.character(body)), tmp)

    # Use base R URL connection for POST
    h <- curl::new_handle()
    curl::handle_setopt(h,
      post = TRUE,
      postfields = as.character(body),
      httpheader = "Content-Type: application/json"
    )
    req <- curl::curl_fetch_memory(url, handle = h)
    jsonlite::fromJSON(rawToChar(req$content), simplifyVector = TRUE)
  }, error = function(e) {
    # Fallback: use system curl
    tmp_body <- tempfile(fileext = ".json")
    tmp_resp <- tempfile(fileext = ".json")
    on.exit(unlink(c(tmp_body, tmp_resp)), add = TRUE)
    writeLines(as.character(body), tmp_body)
    system2("curl", args = c("-s", "-X", "POST",
                              "-H", "'Content-Type: application/json'",
                              "-d", paste0("@", tmp_body),
                              url, "-o", tmp_resp),
            stdout = FALSE, stderr = FALSE)
    jsonlite::fromJSON(tmp_resp, simplifyVector = TRUE)
  })

  response
}
