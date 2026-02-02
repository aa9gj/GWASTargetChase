# GWASTargetChase

Prioritization of GWAS genes using OpenTargets, TWAS-hub and IMPC data.

## Purpose

Complex diseases are disorders that result from multiple genetic variants and genes coupled with influences from the environment. One way to study complex disease is Genome-wide association studies (GWAS). GWAS is an approach to study cohorts (populations) and associate loci (genomic regions containing multiple variants and genes) to complex diseases. However, a main challenge to interpreting GWAS results is 90% of these loci are found within intronic or intergenic regions suggesting they are involved in gene regulation rather than affecting the gene sequence. Therefore, GWAS follow-up (i.e.using multi-omics data and approaches to pin point causal genes) has been used as a tool in human genetics to identify causal genes in complex diseases and therefore hone in on potential therapeutic and intervention targets.

Research funding has been poured into generating resources of GWAS follow-up data in humans given the success this approach has shown in identifying gene targets. However, other species are lagging behind due to the lack of funding. Therefore, the idea to use humans as a model system for feline complex disease could potentially point to genes of interest to feline health. Thus, the purpose of this package is to prepare, and utilize data from human studies to follow up on feline and/or K9 GWAS summary statistics.

## Data Curation

The human GWAS follow-up data is obtained from the following sources:

- [OpenTargets](https://www.opentargets.org/)
- [IMPC](https://www.mousephenotype.org/)
- [PhenomeXcan](http://apps.hakyimlab.org/phenomexcan/)
- [Relevant WebTWAS traits](http://www.webtwas.net/#/diseases)

The data from OpenTargets and IMPC can be prepared using functions within this package. The data from PhenomeXcan and WebTWAS can also be downloaded and prepared using the provided functions.

## Installation

### Prerequisites

Requires R >= 3.5.0. Install Bioconductor dependencies first:

```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "rtracklayer", "S4Vectors"))
```

### Install from GitHub

```r
# Install remotes if not already installed
if (!require("remotes", quietly = TRUE))
  install.packages("remotes")

# Install GWASTargetChase
remotes::install_github("aa9gj/GWASTargetChase")

# Load the package
library(GWASTargetChase)
```

Alternatively, using devtools:

```r
devtools::install_github("aa9gj/GWASTargetChase")
```

## Sample Data

The package includes sample data files for testing and understanding the expected input formats:

```r
# Get path to example GWAS summary statistics
gwas_file <- system.file("extdata", "example_gwas_sumstats.tsv",
                         package = "GWASTargetChase")

# Get path to example GTF file (minimal example)
gtf_file <- system.file("extdata", "example_cat_gtf.gtf",
                        package = "GWASTargetChase")

# Read the example data to see the format
example_gwas <- data.table::fread(gwas_file)
head(example_gwas)
```

## GTF Files

Download the appropriate GTF files for your species. This package depends on Ensembl data structure for gene names and chromosome IDs.

### Cat (Felis catus)

```sh
wget http://ftp.ensembl.org/pub/release-106/gtf/felis_catus/Felis_catus.Felis_catus_9.0.106.gtf.gz
gunzip Felis_catus.Felis_catus_9.0.106.gtf.gz
```

### Dog (Canis lupus familiaris)

```sh
wget http://ftp.ensembl.org/pub/release-106/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf.gz
gunzip Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf.gz
```

### Human (for data preparation)

```sh
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz
gunzip gencode.v40.annotation.gtf.gz
```

## Quick Start

### Option A: One-command download and prepare (recommended)

The `downloadData()` function downloads and processes all required reference
data from IMPC, OpenTargets Platform, and OpenTargets Genetics:

```r
library(GWASTargetChase)

# Download and prepare all reference data (requires wget and internet)
paths <- downloadData(destdir = "GWASTargetChase_data")

# Run analysis using the downloaded data
gwasFollowupMan(sumStats = "your_gwas_sumstats.txt",
                felGTF = "Felis_catus.Felis_catus_9.0.106.gtf",
                pval = 0.00005,
                ResultsPath = "results",
                impc = paths$impc,
                assocOT = paths$assoc,
                l2gOT = paths$l2g)
```

Note: The downloads are large (several GB). Requires `wget` on your system.

### Option B: Manual step-by-step preparation

#### Step 1: Prepare the reference data

First, prepare the required reference datasets:

```r
library(GWASTargetChase)

# 1. Prepare IMPC data
# Download: ftp://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/phenotypeHitsPerGene.csv.gz
IMPCprep(impcData = "phenotypeHitsPerGene.csv",
         ResultsPath = "prepared_data")

# 2. Prepare OpenTargets locus2gene data
# Download l2g and study-index from OpenTargets Genetics
l2gPrep(l2gOT = "path/to/l2g/",
        studyOT = "path/to/study_index.json",
        humanGTF = "gencode.v40.annotation.gtf",
        ResultsPath = "prepared_data")

# 3. Prepare OpenTargets genetic association data
# Download from OpenTargets Platform
geneticAssocPrep(assocOT = "associationByDatasourceDirect.json",
                 diseaseOT = "diseases.json",
                 humanGTF = "gencode.v40.annotation.gtf",
                 ResultsPath = "prepared_data")
```

#### Step 2: Run the analysis

```r
# Use gwasFollowupMan with your prepared data
gwasFollowupMan(sumStats = "your_gwas_sumstats.txt",
                felGTF = "Felis_catus.Felis_catus_9.0.106.gtf",
                pval = 0.00005,
                ResultsPath = "results",
                impc = "prepared_data/impcData",
                assocOT = "prepared_data/disease_target_genetic_associations",
                l2gOT = "prepared_data/l2g_annotated_full")
```

## Functions

There are 9 exported functions within this package: `downloadData`, `gwasFollowup`, `gwasFollowupMan`, `gwasFollowupFull`, `gwasFollowuptest`, `IMPCprep`, `l2gPrep`, `geneticAssocPrep`, and `phenomePrep`

### Input Requirements

Your GWAS summary statistics file must have the following column names:
- `p_wald`: Corrected p-values
- `ps`: SNP position
- `chr`: Chromosome

### gwasFollowup

Main function that uses pre-loaded internal data (if available).

```r
gwasFollowup(sumStats = "/path/to/GWAS/Summary/stats",
             felGTF = "/path/to/feline/or/k9/gtf/file",
             pval = 0.00005,
             ResultsPath = "/path/to/results/directory/")
```

### gwasFollowupMan

Manual version that uses your own prepared data files. **Recommended for most users.**

```r
gwasFollowupMan(sumStats = "/path/to/sumstats",
                felGTF = "/path/to/Felis_catus.Felis_catus_9.0.106.gtf",
                pval = 0.00005,
                ResultsPath = "/path/to/results/",
                impc = "/path/to/impc/results",
                assocOT = "/path/to/disease_target_genetic_associations",
                l2gOT = "/path/to/l2g_annotated_full")
```

### gwasFollowupFull

Extended version with support for PhenomeXcan and TWAS data.

```r
gwasFollowupFull(sumStats = "/path/to/GWAS/summary/stats",
                 felGTF = "/path/to/feline/or/k9/gtf/file",
                 pval = 0.00000005,
                 ResultsPath = "/path/to/results/directory/",
                 phenomePath = "/path/to/PhenomeXcan/data",
                 twasPath = "/path/to/twas/data")
```

### Data Preparation Functions

#### IMPCprep

Prepares IMPC (International Mouse Phenotyping Consortium) data.

```sh
# Download the data first
wget ftp://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-16.0/results/phenotypeHitsPerGene.csv.gz
gunzip phenotypeHitsPerGene.csv.gz
```

```r
IMPCprep(impcData = "/path/to/phenotypeHitsPerGene.csv",
         ResultsPath = "/path/to/results/")
```

#### l2gPrep

Prepares OpenTargets locus2gene data.

```sh
# Download the data first
wget ftp://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/l2g/*
wget ftp://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/lut/study-index/*
cat *.json > study_index.json
```

```r
l2gPrep(l2gOT = "/path/to/l2g/dir/",
        studyOT = "/path/to/study_index.json",
        humanGTF = "/path/to/gencode.v40.annotation.gtf",
        ResultsPath = "/path/to/results/")
```

#### geneticAssocPrep

Prepares OpenTargets genetic association data.

```sh
# Download the data first
wget --recursive --no-parent --no-host-directories --cut-dirs 8 \
  ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/json/associationByDatasourceDirect
cat *.json > associationByDatasourceDirect.json

wget --recursive --no-parent --no-host-directories --cut-dirs 8 \
  ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/json/diseases
cat *.json > diseases.json
```

```r
geneticAssocPrep(assocOT = "/path/to/associationByDatasourceDirect.json",
                 diseaseOT = "/path/to/diseases.json",
                 humanGTF = "/path/to/gencode.v40.annotation.gtf",
                 ResultsPath = "/path/to/results/")
```

#### phenomePrep

Prepares PhenomeXcan data.

```sh
# Download the data first
wget https://zenodo.org/record/3911190/files/fastenloc-torus-rcp.tsv.gz
gunzip fastenloc-torus-rcp.tsv.gz
```

```r
phenomePrep(phenome = "/path/to/fastenloc-torus-rcp.tsv",
            gtf = "/path/to/human/gtf/file",
            ResultsPath = "/path/to/results/")
```

## Output

All main functions (`gwasFollowup`, `gwasFollowupMan`, `gwasFollowupFull`, `gwasFollowuptest`) generate:

- **g2d_results.txt**: Gene-to-disease association scores from OpenTargets
- **l2g_results.txt**: Full locus-to-gene results

Additionally, `gwasFollowup` generates:

- **l2g_filtered_results.txt**: Filtered results with key columns:
  1. `study_id`: GWAS study ID from GWAS Catalog
  2. `gene_id`: Ensembl gene ID
  3. `full_l2g_score`: L2G model probability score
  4. `y_proba_logi_*`: Component scores (distance, interaction, molecularQTL, pathogenicity)
  5. `trait_reported`: GWAS trait annotation
  6. `trait_category`: Type of trait
  7. `gene_name`: Gene symbol
  8. `IMPC_results`: Mouse knockout phenotypes
  9. `rs` and following: SNP information from input GWAS
- **plots.pdf**: LocusZoom-style plots for significant loci

## License

MIT License
