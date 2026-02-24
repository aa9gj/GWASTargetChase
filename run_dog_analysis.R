#!/usr/bin/env Rscript
# =============================================================================
# Dog Body Weight GWAS Analysis using GWASTargetChase
# =============================================================================
#
# Based on published dog body weight GWAS literature:
#   - Hayward et al. 2016, Nature Communications (doi:10.1038/ncomms10460)
#   - Rimbault et al. 2013, Genome Research
#
# Loci: FTO (CFA2), MC4R (CFA1), LEP (CFA20)
# Species: dog (uses Zoonomia orthology for gene translation)
#
# Usage: source("run_dog_analysis.R") or Rscript run_dog_analysis.R
# =============================================================================

library(GWASTargetChase)

# --- Input files (bundled with the package) ---
sumstats <- system.file("extdata", "example_dog_gwas_sumstats.tsv",
                        package = "GWASTargetChase")
gtf      <- system.file("extdata", "example_dog_gtf.gtf",
                        package = "GWASTargetChase")

stopifnot(file.exists(sumstats), file.exists(gtf))

cat("=== Dog Body Weight GWAS — TargetChase Analysis ===\n\n")
cat("Summary stats:", sumstats, "\n")
cat("GTF file:     ", gtf, "\n\n")

# --- Preview input data ---
cat("--- Summary Statistics ---\n")
print(data.table::fread(sumstats))
cat("\n")

# --- Run TargetChase ---
results_dir <- file.path(tempdir(), "dog_analysis_results")
cat("Results will be written to:", results_dir, "\n\n")

results <- TargetChase(
  sumStats    = sumstats,
  felGTF      = gtf,
  species     = "dog",
  pval        = 5e-8,
  ResultsPath = results_dir
)

# --- Show results ---
cat("\n=== Output Files ===\n")
output_files <- list.files(results_dir, recursive = TRUE, full.names = TRUE)
for (f in output_files) {
  cat(" ", f, "\n")
}

cat("\n=== Gene-to-Disease Results ===\n")
g2d <- read.delim(file.path(results_dir, "g2d_results.txt"))
if (nrow(g2d) > 0) {
  print(head(g2d, 20))
} else {
  cat("No gene-disease associations found.\n")
}

impc_file <- file.path(results_dir, "impc_results.txt")
if (file.exists(impc_file)) {
  cat("\n=== IMPC Mouse Phenotype Results ===\n")
  impc <- read.delim(impc_file)
  if (nrow(impc) > 0) {
    print(head(impc, 20))
  } else {
    cat("No IMPC phenotype data found.\n")
  }
}

cat("\n=== Done! ===\n")
cat("Full results are in:", results_dir, "\n")
