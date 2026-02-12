# GWASTargetChase

Gene prioritization for GWAS results using OpenTargets and IMPC data.

## Overview

GWASTargetChase helps you identify and prioritize candidate genes from GWAS results by integrating:

- **OpenTargets Platform**: Gene-disease associations and genetic evidence scores
- **IMPC** (International Mouse Phenotyping Consortium): Mouse knockout phenotypes

Supports **human**, **mouse**, **cat**, and **dog** GWAS analysis with automatic GTF annotation download.

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
# Install remotes if needed
if (!require("remotes", quietly = TRUE))
  install.packages("remotes")

# Install GWASTargetChase
remotes::install_github("aa9gj/GWASTargetChase")

# Install optional dependencies for API access
install.packages(c("httr", "rappdirs", "R.utils"))
```

## Quick Start

### One-command analysis

```r
library(GWASTargetChase)

# Run prioritization analysis
# GTF file is automatically downloaded and cached
results <- run_prioritization_analysis(
  sumStats = "my_gwas_results.txt",
  species = "cat",      # Options: "human", "mouse", "cat", "dog"
  pval = 5e-8,
  output_dir = "results/"
)

# View prioritized genes
head(results$gene_summary)
```

### Input Requirements

Your GWAS summary statistics file must have these columns:
- `chr`: Chromosome (e.g., "1", "X", or "chr1")
- `ps`: SNP position (base pairs)
- `p_wald`: P-value

Example:
```
chr	rs	ps	p_wald	af	beta	se
1	rs1421085	53786615	2.5e-12	0.42	0.35	0.05
1	rs9939609	53800954	5.1e-11	0.41	0.32	0.05
```

## What it does

`run_prioritization_analysis()` performs these steps:

1. **Downloads GTF** (if needed): Automatically fetches species-specific gene annotations from Ensembl and caches them locally
2. **Reads GWAS data**: Filters for significant SNPs based on your p-value threshold
3. **Identifies genes**: Finds protein-coding genes within 1Mb of significant loci
4. **Fetches OpenTargets data**: Queries the API for gene-disease associations
5. **Fetches IMPC data**: Queries the API for mouse knockout phenotypes
6. **Outputs results**: Creates integrated result files

## Output Files

| File | Description |
|------|-------------|
| `gene_summary.txt` | Main results with all genes, disease associations, and IMPC phenotypes |
| `gene_disease_associations.txt` | Full OpenTargets gene-disease data |
| `impc_phenotypes.txt` | IMPC mouse knockout phenotype data |

### gene_summary.txt columns

- `gene_name`: Gene symbol
- `MGI_Gene_id`: Mouse Gene ID (from IMPC)
- `#phenotype_hits`: Number of mouse knockout phenotypes
- `Phenotype_Hits`: List of observed phenotypes
- `n_disease_associations`: Number of disease associations in OpenTargets
- `top_diseases`: Top 3 associated diseases
- `max_genetic_assoc_score`: Highest genetic association score
- `genecards_url`: Link to GeneCards for more info

## Additional Functions

### Fetch data separately

```r
# Get GTF file for a species (downloads if needed)
gtf_path <- get_gtf("dog")

# Fetch OpenTargets gene-disease associations
assoc <- fetch_gene_disease_associations(
  gene_names = c("FTO", "MC4R", "BRCA1")
)

# Fetch IMPC mouse phenotypes
impc <- fetch_impc_data(
  gene_names = c("FTO", "MC4R", "BRCA1")
)

# Check OpenTargets API connection
check_opentargets_connection()
```

### Input validation

```r
# Validate GWAS file format
validate_gwas_input("my_gwas.txt")

# Check for significant SNPs
check_significant_snps("my_gwas.txt", pval = 5e-8)

# Validate GTF file
validate_gtf_file("my_annotation.gtf")
```

## Supported Species

| Species | GTF Source |
|---------|------------|
| Human | Ensembl GRCh38.110 |
| Mouse | Ensembl GRCm39.110 |
| Cat | Ensembl Felis_catus_9.0.110 |
| Dog | Ensembl ROS_Cfam_1.0.110 |

GTF files are automatically downloaded on first use and cached in `~/.cache/GWASTargetChase/gtf/`.

## Example with sample data

```r
library(GWASTargetChase)

# Use included sample data
gwas_file <- system.file("extdata", "example_gwas_sumstats.tsv",
                         package = "GWASTargetChase")

# Run with sample data
results <- run_prioritization_analysis(
  sumStats = gwas_file,
  species = "cat",
  pval = 1e-6,
  output_dir = tempdir()
)

# Explore results
print(results$genes)
View(results$gene_summary)
```

## License

MIT License
