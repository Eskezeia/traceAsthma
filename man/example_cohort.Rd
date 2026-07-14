\name{example_cohort}
\alias{example_cohort}
\docType{data}
\title{Synthetic example cohort for the TRACE-Asthma pipeline}
\format{
A list with elements:
\describe{
  \item{dosage}{Numeric matrix, 5 patients x 2 SNPs.}
  \item{methylation}{Numeric matrix, 2 CpGs x 5 patients.}
  \item{expression}{Numeric matrix, 4 genes x 5 patients.}
  \item{regulon_db}{Data frame, TF-target regulon table (3 TFs).}
  \item{clinical}{Data frame, one row per patient: \code{patient_id},
    \code{asthma}, \code{age}, \code{sex}, \code{bmi}, \code{smoking}.}
  \item{covariates}{Data frame, eQTM-model covariates (age, sex, PC1, PC2).}
  \item{eqtm_pairs}{Data frame, candidate CpG-gene pairs for
    \code{discover_eqtm()}.}
  \item{same_chr}{Logical vector flagging same-chromosome CpG-gene pairs.}
}
}
\source{
Simulated.
}
\description{
A fully synthetic, five-patient toy dataset used throughout the package
examples and vignette to demonstrate the end-to-end pipeline in seconds.
Contains two SNPs, two CpGs, four target genes, and a three-TF regulon --
deliberately tiny so every intermediate score can be hand-verified.
\strong{Not derived from any real cohort; for demonstration purposes only.}
}
\examples{
data(example_cohort)
str(example_cohort, max.level = 1)
}
\keyword{datasets}
