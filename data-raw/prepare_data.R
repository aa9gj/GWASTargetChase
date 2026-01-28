# Data preparation script for GWASTargetChase package
#
# This script creates sample datasets for testing and demonstrates how to
# prepare the required data files for this package.
#
# Run this script to generate the sample data files:
# source("data-raw/prepare_data.R")

library(data.table)

# ==============================================================================
# SAMPLE DATA FOR TESTING
# ==============================================================================
# These are example datasets with realistic structure that allow users to
# test the package functionality.

# ------------------------------------------------------------------------------
# 1. Sample IMPC Data
# ------------------------------------------------------------------------------
# IMPC (International Mouse Phenotyping Consortium) phenotype data
# Real data from: ftp://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/

impc <- data.frame(

  gene_name = c("BRCA1", "BRCA2", "TP53", "EGFR", "KRAS", "MYC", "PTEN", "RB1",
                "APC", "CDKN2A", "FTO", "MC4R", "LEP", "LEPR", "PCSK9", "APOE",
                "LDLR", "HMGCR", "NPC1", "ABCA1", "TCF7L2", "PPARG", "IRS1",
                "KCNJ11", "SLC30A8", "HNF1A", "HNF4A", "GCK", "INS", "MTNR1B"),
  Phenotype_Hits = c(
    "abnormal DNA damage response|abnormal tumor incidence",
    "abnormal DNA damage response|embryonic lethality",
    "abnormal tumor incidence|premature death",
    "abnormal lung morphology|abnormal skin morphology",
    "embryonic lethality|abnormal embryo size",
    "embryonic lethality|abnormal cell proliferation",
    "abnormal tumor incidence|increased body weight",
    "embryonic lethality|abnormal eye morphology",
    "abnormal intestinal morphology|abnormal tumor incidence",
    "abnormal tumor incidence|abnormal melanocyte morphology",
    "increased body weight|increased fat mass",
    "increased body weight|abnormal food intake",
    "increased body weight|abnormal glucose homeostasis",
    "increased body weight|abnormal glucose homeostasis",
    "decreased circulating cholesterol level",
    "abnormal lipid homeostasis|abnormal learning",
    "increased circulating cholesterol level",
    "embryonic lethality",
    "abnormal lipid homeostasis|abnormal liver morphology",
    "abnormal lipid homeostasis",
    "abnormal glucose homeostasis|decreased body weight",
    "abnormal glucose homeostasis|abnormal fat morphology",
    "abnormal glucose homeostasis",
    "abnormal glucose homeostasis|abnormal insulin secretion",
    "abnormal glucose homeostasis|abnormal insulin secretion",
    "abnormal glucose homeostasis|abnormal liver morphology",
    "abnormal glucose homeostasis|abnormal liver morphology",
    "abnormal glucose homeostasis|postnatal lethality",
    "abnormal glucose homeostasis|decreased body weight",
    "abnormal glucose homeostasis|abnormal circadian rhythm"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 2. Sample OpenTargets Genetic Association Data
# ------------------------------------------------------------------------------
# Gene-disease associations from OpenTargets Platform
# Real data from: ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/

disease_target_genetic_association <- data.table(
  gene_name = c("FTO", "FTO", "MC4R", "MC4R", "LEP", "TCF7L2", "TCF7L2",
                "PPARG", "KCNJ11", "SLC30A8", "PCSK9", "LDLR", "APOE",
                "BRCA1", "BRCA2", "TP53", "EGFR", "KRAS"),
  gene_id = c("ENSG00000140718", "ENSG00000140718", "ENSG00000166603",
              "ENSG00000166603", "ENSG00000174697", "ENSG00000148737",
              "ENSG00000148737", "ENSG00000132170", "ENSG00000187486",
              "ENSG00000164400", "ENSG00000169174", "ENSG00000130164",
              "ENSG00000130203", "ENSG00000012048", "ENSG00000139618",
              "ENSG00000141510", "ENSG00000146648", "ENSG00000133703"),
  disease_id = c("EFO_0001073", "MONDO_0011122", "EFO_0001073", "MONDO_0011122",
                 "EFO_0001073", "EFO_0000400", "MONDO_0005148", "EFO_0000400",
                 "EFO_0000400", "EFO_0000400", "EFO_0000408", "EFO_0000408",
                 "EFO_0000249", "EFO_0000305", "EFO_0000305", "EFO_0000311",
                 "EFO_0000571", "EFO_0000571"),
  disease_name = c("obesity", "obesity disorder", "obesity", "obesity disorder",
                   "obesity", "type 2 diabetes", "diabetes mellitus",
                   "type 2 diabetes", "type 2 diabetes", "type 2 diabetes",
                   "hypercholesterolemia", "hypercholesterolemia",
                   "Alzheimer's disease", "breast carcinoma", "breast carcinoma",
                   "cancer", "lung carcinoma", "lung carcinoma"),
  score = c(0.85, 0.82, 0.78, 0.75, 0.72, 0.91, 0.88, 0.76, 0.81, 0.79,
            0.93, 0.89, 0.84, 0.95, 0.94, 0.88, 0.82, 0.79)
)

# ------------------------------------------------------------------------------
# 3. Sample OpenTargets Locus2Gene Data
# ------------------------------------------------------------------------------
# L2G (Locus to Gene) scores from OpenTargets Genetics
# Real data from: ftp://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/

l2g_annotated_full <- data.frame(
  study_id = c("GCST004773", "GCST004773", "GCST002783", "GCST002783",
               "GCST000679", "GCST000679", "GCST001633", "GCST001633",
               "GCST90002232", "GCST90002232", "GCST004988", "GCST004988"),
  variant_id = c("16_53786615_C_T", "16_53800954_G_A", "18_57829135_G_A",
                 "18_57838141_C_T", "10_114748339_G_A", "10_114758349_C_T",
                 "6_161006023_G_A", "6_161010276_C_T", "1_55039774_G_T",
                 "1_55051215_G_A", "11_92708710_C_T", "11_92673828_G_A"),
  gene_id = c("ENSG00000140718", "ENSG00000140718", "ENSG00000166603",
              "ENSG00000166603", "ENSG00000148737", "ENSG00000148737",
              "ENSG00000187486", "ENSG00000187486", "ENSG00000169174",
              "ENSG00000169174", "ENSG00000118046", "ENSG00000118046"),
  gene_name = c("FTO", "FTO", "MC4R", "MC4R", "TCF7L2", "TCF7L2",
                "KCNJ11", "KCNJ11", "PCSK9", "PCSK9", "STK33", "STK33"),
  y_proba_full_model = c(0.92, 0.88, 0.85, 0.82, 0.91, 0.87, 0.79, 0.76,
                         0.94, 0.90, 0.72, 0.68),
  y_proba_logi_distance = c(0.88, 0.85, 0.82, 0.79, 0.87, 0.84, 0.75, 0.72,
                            0.91, 0.88, 0.68, 0.65),
  y_proba_logi_interaction = c(0.75, 0.72, 0.68, 0.65, 0.78, 0.74, 0.62, 0.58,
                               0.82, 0.78, 0.55, 0.52),
  y_proba_logi_molecularQTL = c(0.85, 0.82, 0.78, 0.75, 0.84, 0.80, 0.72, 0.68,
                                0.88, 0.85, 0.65, 0.62),
  y_proba_logi_pathogenicity = c(0.45, 0.42, 0.52, 0.48, 0.38, 0.35, 0.55, 0.52,
                                 0.62, 0.58, 0.42, 0.38),
  trait_reported = c("Body mass index", "Body mass index", "Body mass index",
                     "Body mass index", "Type 2 diabetes", "Type 2 diabetes",
                     "Type 2 diabetes", "Type 2 diabetes",
                     "LDL cholesterol levels", "LDL cholesterol levels",
                     "Body mass index", "Body mass index"),
  trait_category = c("Body measurement", "Body measurement", "Body measurement",
                     "Body measurement", "Metabolic disorder", "Metabolic disorder",
                     "Metabolic disorder", "Metabolic disorder",
                     "Lipid measurement", "Lipid measurement",
                     "Body measurement", "Body measurement"),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Save internal datasets to data/ directory
# ------------------------------------------------------------------------------
cat("Saving internal datasets...\n")

# Create data directory if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

save(impc, file = "data/impc.rda", compress = "xz")
save(disease_target_genetic_association,
     file = "data/disease_target_genetic_association.rda", compress = "xz")
save(l2g_annotated_full, file = "data/l2g_annotated_full.rda", compress = "xz")

cat("Internal datasets saved to data/ directory\n")

# ------------------------------------------------------------------------------
# 4. Sample GWAS Summary Statistics
# ------------------------------------------------------------------------------
# Create example GWAS summary statistics file for testing

example_gwas <- data.frame(
  chr = c(rep("1", 5), rep("2", 5), rep("16", 5)),
  rs = c("rs1421085", "rs9939609", "rs1558902", "rs17817449", "rs8050136",
         "rs10938397", "rs6548238", "rs2867125", "rs543874", "rs7561317",
         "rs9930506", "rs7190492", "rs9941349", "rs1121980", "rs9939609"),
  ps = c(53786615, 53800954, 53803574, 53813367, 53816275,
         45182527, 622827, 624905, 622348, 637231,
         53786615, 53790000, 53795000, 53798000, 53800954),
  p_wald = c(2.5e-12, 5.1e-11, 8.3e-10, 1.2e-8, 3.4e-7,
             1.8e-9, 4.2e-8, 6.7e-7, 2.1e-6, 8.9e-6,
             3.2e-15, 7.8e-12, 1.5e-9, 4.3e-7, 9.1e-6),
  af = c(0.42, 0.41, 0.43, 0.40, 0.39,
         0.35, 0.28, 0.31, 0.33, 0.36,
         0.44, 0.42, 0.41, 0.38, 0.40),
  beta = c(0.35, 0.32, 0.28, 0.24, 0.19,
           0.27, 0.22, 0.18, 0.15, 0.12,
           0.41, 0.36, 0.29, 0.21, 0.16),
  se = c(0.05, 0.05, 0.05, 0.05, 0.05,
         0.05, 0.05, 0.05, 0.05, 0.05,
         0.05, 0.05, 0.05, 0.05, 0.05),
  stringsAsFactors = FALSE
)

# Create inst/extdata directory if it doesn't exist
if (!dir.exists("inst/extdata")) {
  dir.create("inst/extdata", recursive = TRUE)
}

# Save example GWAS data
write.table(example_gwas, file = "inst/extdata/example_gwas_sumstats.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Example GWAS summary statistics saved to inst/extdata/\n")

# ==============================================================================
# INSTRUCTIONS FOR PREPARING REAL DATA
# ==============================================================================
#
# The sample data above is for testing purposes only. For real analyses,
# you need to prepare data from the following sources:
#
# 1. IMPC Data:
#    Download: ftp://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-16.0/results/phenotypeHitsPerGene.csv.gz
#    Prepare using: IMPCprep() function
#
# 2. OpenTargets Genetic Association:
#    Download association and disease files from OpenTargets Platform
#    Prepare using: geneticAssocPrep() function
#
# 3. OpenTargets Locus2Gene:
#    Download L2G and study index files from OpenTargets Genetics
#    Prepare using: l2gPrep() function
#
# 4. GTF Files:
#    Human: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz
#    Cat: http://ftp.ensembl.org/pub/release-106/gtf/felis_catus/Felis_catus.Felis_catus_9.0.106.gtf.gz
#    Dog: http://ftp.ensembl.org/pub/release-106/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf.gz
#
# ==============================================================================

cat("\nData preparation complete!\n")
cat("To use the sample data, the package needs to be rebuilt.\n")
