\name{example_data}
\alias{example_data}
\alias{example_gwas_sumstats}
\alias{example_cat_gtf}
\title{Example Data Files}
\description{
Example data files included with the GWASTargetChase package for testing and demonstration purposes.
}
\details{
The package includes the following example files in \code{inst/extdata/}:

\describe{
  \item{example_gwas_sumstats.tsv}{
    A sample GWAS summary statistics file with the required column format:
    \itemize{
      \item \code{chr}: Chromosome
      \item \code{rs}: SNP rsID
      \item \code{ps}: SNP position
      \item \code{p_wald}: P-value (corrected)
      \item \code{af}: Allele frequency
      \item \code{beta}: Effect size
      \item \code{se}: Standard error
    }
  }
  \item{example_cat_gtf.gtf}{
    A minimal GTF file excerpt demonstrating the expected format for feline/canine GTF files.
    For real analyses, download the full GTF files from Ensembl.
  }
}
}
\examples{
# Get path to example GWAS summary statistics
gwas_file <- system.file("extdata", "example_gwas_sumstats.tsv",
                         package = "GWASTargetChase")

# Get path to example GTF file
gtf_file <- system.file("extdata", "example_cat_gtf.gtf",
                        package = "GWASTargetChase")

# Read the example GWAS data
example_gwas <- data.table::fread(gwas_file)
head(example_gwas)
}
\seealso{
\code{\link{gwasFollowup}}, \code{\link{gwasFollowupMan}}
}
\keyword{datasets}
