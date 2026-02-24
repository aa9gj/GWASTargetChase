# GWASTargetChase

Gene prioritization for GWAS results using OpenTargets and IMPC data.

## Overview

GWASTargetChase helps you identify and prioritize candidate genes from GWAS results by integrating:

- **OpenTargets Platform**: Gene-disease associations and genetic evidence scores
- **IMPC** (International Mouse Phenotyping Consortium): Mouse knockout phenotypes
- **Zoonomia**: Cross-species orthology for translating genes between human, mouse, cat, and dog

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
```

## Quick Start

### TargetChase (API-based)

The main function queries OpenTargets and IMPC APIs directly:

```r
library(GWASTargetChase)

TargetChase(
  sumStats = "my_gwas_results.txt",
  felGTF = "my_species.gtf",
  species = "cat",      # Options: "human", "mouse", "cat", "dog"
  pval = 5e-8,
  ResultsPath = "results/"
)
```

### TargetChaseManual (local data)

For use with pre-downloaded OpenTargets and IMPC data:

```r
# Download data once
paths <- downloadData(destdir = "my_data")

# Run with local files
TargetChaseManual(
  sumStats = "my_gwas_results.txt",
  felGTF = "my_species.gtf",
  species = "cat",
  pval = 5e-8,
  ResultsPath = "results/",
  impc = paths$impc,
  assocOT = paths$assoc,
  l2gOT = paths$l2g
)
```

### Input Requirements

Your GWAS summary statistics file must have these columns:
- `chr`: Chromosome (e.g., "1", "X", or "chr1")
- `ps`: SNP position (base pairs)
- `p_wald`: P-value
- `gene`: Closest gene name per locus (required for `TargetChase`)

Example:
```
chr	rs	ps	p_wald	af	beta	se	gene
1	rs1421085	53786615	2.5e-12	0.42	0.35	0.05	FTO
1	rs9939609	53800954	5.1e-11	0.41	0.32	0.05	FTO
```

## Output Files

| File | Description |
|------|-------------|
| `g2d_results.txt` | Gene-disease associations from OpenTargets with IMPC phenotypes |
| `impc_results.txt` | IMPC mouse knockout phenotype data |
| `l2g_results.txt` | Locus-to-gene scores (TargetChaseManual only) |
| `plots.pdf` | LocusZoom plots per significant locus (TargetChase only) |

## Example with sample data

```r
library(GWASTargetChase)

# Use included sample data
gwas_file <- system.file("extdata", "example_gwas_sumstats.tsv",
                         package = "GWASTargetChase")
gtf_file <- system.file("extdata", "example_cat_gtf.gtf",
                        package = "GWASTargetChase")

# Run with sample data
TargetChase(
  sumStats = gwas_file,
  felGTF = gtf_file,
  species = "cat",
  pval = 5e-8,
  ResultsPath = tempdir()
)
```

## License

MIT License
