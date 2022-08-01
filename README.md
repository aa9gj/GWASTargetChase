---
title: "HillsGWASfollowup Package Documentation"
author: "Arby Abood"
date: "`r format(Sys.time(), '%B %d, %Y')`"
output:
  html_document:
    toc: true
    toc_depth: 3
    theme: united
    number_sections: true
    highlight: tango
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Purpose

Complex diseases are disorders that results from multiple genetic variants and genes coupled with influences from the environment. One way to study complex disease is Genome-wide association studies (GWAS). GWAS is an approach to study cohorts (populations) and associate loci (genomic regions containing multiple variants and genes) to complex diseases. However, a main challenge to interpreting GWAS results is 90% of these loci are found within intronic or intergenic regions suggesting they are involved in gene regulation rather than affecting the gene sequence. Therefore, GWAS follow-up (i.e.using multi-omics data and approaches to pin point causal genes) has been used as a tool in human genetics to identify causal genes in complex diseases and therefore hone in on potential therapeutic and intervention targets. \
Research funding has been poured into generating resources of GWAS follow-up data in humans given the success this approach has shown in identifying gene targets. However, other species are lagging behind due to the lack of funding. Therefore, the idea to use humans as a model system for feline complex disease could potentially point to genes of interest to feline health. Thus, the purpose of this package is to prepare, and utilize data from human studies to follow up on feline and/or K9 GWAS summary statistics. 

# Data Curation

The human GWAS follow-up data is obtained from the following sources: \

- OpenTargets <https://www.opentargets.org/>\
- IMPC <https://www.mousephenotype.org/>\
- PhenomeXcan <http://apps.hakyimlab.org/phenomexcan/>\
- Relevant WebTWAS traits <http://www.webtwas.net/#/diseases>\

The data from OpenTargets and IMPC has been wrangled and loaded onto the package (check functions section for reproducing the data). The data from PhenomeXcan and WebTWAS can be downloaded and prepared using functions within this package. 

# Installation and Loading

```{r, eval=FALSE }
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
install.packages("data.table")
install.packages("dplyr")
BiocManager::install("rtracklayer")
install.packages("arrow")
install.packages("jsonlite")
BiocManager::install("S4Vectors")
install.packages('tibble')
install.packages("magrittr")
install.packages("ggplot2")
install.packages("devtools")
library("devtools")
devtools::install_github("thomasp85/patchwork")
p <- c("GenomicRanges", "data.table", "dplyr", "rtracklayer", "arrow", "jsonlite", "S4Vectors", "tibble", "margrittr", "ggplot2", "patchwork")
lapply(p, require, character.only = T)
install.packages("/Users/Jeff Brockman/Documents/R/HillsGWASfollowup_1.0.0.0.tar.gz", repos = NULL, type = "source")
library(HillsGWASfollowup)
```

# Recommended GTF files

Download the following GTF files. Please note this package depends on Ensemble data structure for gene names and chromosome IDs.

```{sh, eval=FALSE}
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz
gunzip gencode.v40.annotation.gtf.gz
wget http://ftp.ensembl.org/pub/release-106/gtf/felis_catus/Felis_catus.Felis_catus_9.0.106.gtf.gz
gunzip Felis_catus.Felis_catus_9.0.106.gtf.gz
```

# Functions

There are 7 exported functions within this package  `gwasFollowup`, `gwasFollowupMan`, `gwasFollowupFull`, `IMPCprep`, `l2gPrepMan`, `geneticAssocPrepMan`, and `phenomePrep`

Please ensure that your sumStats file has the following column names exactly as such: (The rest can stay as is)\
- `p_wald`: This would be the corrected p-values\
- `ps`: SNP position\
- `chr`: chromosome\

## gwasFollouwp

This is the main function of this package. It takes GWAS summary stats and GTF file. Internally, it uses two OpenTarget files and an IMPC phenotype file per gene  already lazyloaded into the package by the author

```{r, eval=FALSE}
gwasFollowup(sumStats = "/path/to/GWAS/Summary/stats", 
             felGTF = "/path/to/feline/or/k9/gef/file", pval = 0.00005,
             ResultsPath = "/path/to/results/directory/")
```

## gwasFollowupFull

This function takes GWAS summary stats and GTF file. Internally, it uses two OpenTarget files and an IMPC phenotype file per gene  already lazyloaded onto the package by the author. The additional features here are providing the function with phenomeXcan data (default is NULL if no path is provided) and TWAS unprepared data (default is NULL if no path is provided)

```{r, eval= FALSE}
gwasFollowupFull(sumStats = "/path/to/GWAS/summary/stats",
             felGTF = "/path/to/feline/or/k9/gef/file", pval = 0.00000005,
             ResultsPath = "/path/to/results/directory/", 
             phenomePath = "/path/to/PhenomeXcan/data", 
             twasPath = "/path/to/twas/data")
```


## IMPCprep

You would need to download the relevant IMPC data from the command line
```{sh, eval=FALSE}
wget ftp://ftp.ebi.ac.uk/pub/databases/impc//all-data-releases/release-16.0/results/phenotypeHitsPerGene.csv.gz
gunzip phenotypeHitsPerGene.csv.gz
```

In R/Rstudio, use the following command to prepare the IMPC data for further analysis 
```{r, eval=FALSE}
IMPCprep(impcData = "/path/to/phenotypeHitsPerGene.csv", ResultsPath = "/path/to/results/")
```

## l2gPrepMan

If you don't want to use the lazyloaded data (in case of an updated release), you'll have to prepare the data yourself. This function prepares locus2gene data from OpenTargets (OT). You would need to download the l2g and study index data from OT using the following commands in the command line.

```{sh, eval=FALSE}
wget ftp://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/l2g/*
wget ftp://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/lut/study-index/*
cat *.json > study_index.json
```

In R/Rstudio, use the following command to prepare the locus2gene data for further analysis 

```{r, eval=FALSE}
l2gPrepMan(l2gOT = "/path/to/l2g/dir/l2g/",
        studyOT = "/path/to/study_index.json",
        humanGTF = "/path/to/humangtf/e.g.gencode.v40.annotation.gtf.gz",
        ResultsPath = "/path/to/results/dir/")
```

## geneticAssocPrep

If you don't want to use the lazyloaded data (in case of an updated release), you'll have to prepare the data yourself. This function prepares the target data from OpenTargets (OT).This information is not as reliable as their locus2gene model. You would need to download the target-disease associations (direct) and disease annotation data from OT using the following commands in the command line.

```{sh, eval=FALSE}
wget --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/json/associationByDatasourceDirect
cat *.json > associationByDatasourceDirect.json
wget --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/json/diseases
cat *.json > diseases.json
```

In R/Rstudio, use the following command to prepare the genetic association data for further analysis 
```{r, eval=FALSE}
geneticAssocPrep(assocOT = "/path/to/dir/associationByDatasourceDirect.json",
                 diseaseOT = "/path/to/dir/diseases.json",
                 humanGTF = "/path/to/e.ggencode.v40.annotation.gtf.gz",
                 ResultsPath = "/path/to/res/")
```

## gwasFollowupMan

Once the data has been prepared by the two functions above, you can use the function below

```{r, eval=FALSE}
gwasFollowupMan(sumStats = "/path/to/sumstats",
                felGTF = "/path/to/Felis_catus.Felis_catus_9.0.106.gtf",
                pval = e.g 0.00005, ResultsPath = "/path/to/resutls/",
                impc = "/path/to/impc/results",
                assocOT = "/path/to/disease_target_genetic_associations",
                l2gOT = "/path/to/l2g_annotated_full")
```

## phenomePrep

Columns containing the data from phenomeXcan analysis which uses a Bayesian colocalization framework through FastEnloc rather than Coloc. `data has been downloaded in bulk and provided in the package folder` from here: 

```{sh, eval=FALSE}
wget https://zenodo.org/record/3911190/files/fastenloc-torus-rcp.tsv.gz
gunzip fastenloc-torus-rcp.tsv.gz
```

```{r, eval=FALSE}
phenomePrep(phenome = "/path/to/downloaded/data",
            gtf = "/path/to/human/gtf/file", 
            ResultsPath = "/path/to/results/"
```

# Output 

If you use `gwasFollowup` function, here is the output you would get: \ 

- **g2d_results**: A score of how "targetable" the target is from OpenTargets\
- **l2g_results**: Full results section\
- **l2g_filtered**: Interesting results and includes the following columns:\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1. study_id: GWAS study id from GWAS catalog\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2. gene_id: Ensembl gene ID\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3. full_l2g_score: L2G model score from OpenTargets which is explained as the probability of the gene being causal at that locus. This score is calculated from the next four columns in the output (find out more about this here <https://genetics-docs.opentargets.org/our-approach/prioritising-causal-genes-at-gwas-loci-l2g>\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4. y_proba_logi_distance, y_proba_logi_interaction, y_proba_logi_molecularQTL, y_proba_logi_pathogenicity: The next four columns that are used to calculate the full model\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5. trait_reported: The annotation for the GWAS study id\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6. trait_category: The type of trait being studied\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;7. gene_name: Gene name from Ensembl\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8. IMPC_results: Functional impact of the gene knockout in mice\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;9. rs (and all following columns): SNP information from feline GWAS summary statistics\ 
- **plots**: LocusZoom plots for all significant loci (p_wald > 5)
