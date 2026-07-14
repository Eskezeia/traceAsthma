\name{validate_matrix}
\alias{validate_matrix}
\alias{validate_clinical}
\alias{validate_cohort_alignment}
\title{Validate real-world datasets before analysis}
\usage{
validate_matrix(x, type = c("methylation", "expression", "dosage"),
  allow_mvalues = FALSE, allow_imputed = TRUE, max_missing_frac = 0.2,
  min_samples = 3)

validate_clinical(clinical, required_cols = "asthma", id_col = "patient_id")

validate_cohort_alignment(..., id_col = "patient_id")
}
\arguments{
\item{x}{The matrix to validate.}
\item{type}{One of \code{"methylation"} (values expected in [0,1] unless
\code{allow_mvalues}), \code{"expression"} (non-negative), or \code{"dosage"} (values
in \{0,1,2\}, or [0,2] if \code{allow_imputed}).}
\item{allow_mvalues}{Logical; if \code{TRUE}, methylation values outside [0,1]
are permitted (treated as M-values). Default \code{FALSE}.}
\item{allow_imputed}{Logical; if \code{TRUE} (default), non-integer dosages in
[0,2] are permitted.}
\item{max_missing_frac}{Maximum allowed fraction of missing values. Default 0.2.}
\item{min_samples}{Minimum number of columns (subjects) required. Default 3.}
\item{clinical}{Data frame to validate.}
\item{required_cols}{Character vector of columns that must be present.
Default \code{"asthma"}.}
\item{id_col}{Character, subject-ID column name, used to check for
duplicates (and, in \code{validate_cohort_alignment()}, to extract subject
IDs from data-frame arguments). Default \code{"patient_id"}.}
\item{...}{Named matrices/data frames to check for consistent, identically
ordered subject sets.}
}
\value{
Invisibly \code{TRUE} if all checks pass; each function throws an informative,
specific error (naming the offending rows/columns/subjects, truncated to
the first 10) otherwise.
}
\description{
Structural and value-range validation for real-world datasets before any
analysis function runs, so a malformed dataset fails immediately with an
actionable message rather than producing a silent downstream error or an
incorrect result. \code{validate_matrix()} checks a methylation/expression/
dosage matrix's structure, missingness, and value ranges.
\code{validate_clinical()} checks a clinical data frame's required columns,
outcome coding, and plausible covariate ranges. \code{validate_cohort_alignment()}
checks that a set of matrices/tables (e.g. \code{methylation}, \code{expression},
\code{clinical}) describe an identical, identically-ordered subject set --
the single most common real-world data assembly bug.
}
\examples{
m <- matrix(runif(20, 0, 1), nrow = 4, dimnames = list(paste0("cg", 1:4), paste0("P", 1:5)))
validate_matrix(m, type = "methylation")

cl <- data.frame(patient_id = paste0("P", 1:5), asthma = c(0,1,0,1,1))
validate_clinical(cl)
}
