# Data preparation script for HillsGWASfollowup package
#
# This script documents how to prepare the required data files for this package.
# The data files are required for the main functions (gwasFollowup, gwasFollowupFull,
# gwasFollowuptest) to work properly.
#
# If you don't have these data files, use gwasFollowupMan() instead with your own
# prepared data files.

# ==============================================================================
# 1. IMPC (International Mouse Phenotyping Consortium) Data
# ==============================================================================
#
# Download from: ftp://ftp.ebi.ac.uk/pub/databases/impc//all-data-releases/release-16.0/results/phenotypeHitsPerGene.csv.gz
#
# Preparation code:
# impc_raw <- data.table::fread("phenotypeHitsPerGene.csv")
# impc <- impc_raw[, .(gene_name = marker_symbol, Phenotype_Hits = phenotype_hits)]
# usethis::use_data(impc, overwrite = TRUE)

# ==============================================================================
# 2. OpenTargets Genetic Association Data (disease_target_genetic_association)
# ==============================================================================
#
# Download from:
# wget --recursive --no-parent --no-host-directories --cut-dirs 8 \
#   ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/json/associationByDatasourceDirect
# wget --recursive --no-parent --no-host-directories --cut-dirs 8 \
#   ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/json/diseases
#
# Then prepare using geneticAssocPrep() function from this package

# ==============================================================================
# 3. OpenTargets Locus2Gene Data (l2g_annotated_full)
# ==============================================================================
#
# Download from:
# wget ftp://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/l2g/*
# wget ftp://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/lut/study-index/*
#
# Then prepare using l2gPrep() function from this package

# ==============================================================================
# Creating placeholder data (minimal example for package installation)
# ==============================================================================
# Note: These are minimal placeholders. For actual analysis, you need to prepare
# the real data files using the instructions above.

# Create placeholder IMPC data
impc <- data.frame(
  gene_name = character(0),
  Phenotype_Hits = character(0),
  stringsAsFactors = FALSE
)

# Create placeholder disease_target_genetic_association data
disease_target_genetic_association <- data.table::data.table(
  gene_name = character(0),
  disease_id = character(0),
  disease_name = character(0),
  score = numeric(0)
)

# Create placeholder l2g_annotated_full data
l2g_annotated_full <- data.frame(
  study_id = character(0),
  gene_id = character(0),
  gene_name = character(0),
  y_proba_full_model = numeric(0),
  y_proba_logi_distance = numeric(0),
  y_proba_logi_interaction = numeric(0),
  y_proba_logi_molecularQTL = numeric(0),
  y_proba_logi_pathogenicity = numeric(0),
  trait_reported = character(0),
  trait_category = character(0),
  stringsAsFactors = FALSE
)

# Save placeholder data
# usethis::use_data(impc, overwrite = TRUE)
# usethis::use_data(disease_target_genetic_association, overwrite = TRUE)
# usethis::use_data(l2g_annotated_full, overwrite = TRUE)
