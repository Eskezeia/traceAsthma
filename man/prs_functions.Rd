\name{compute_prs}
\alias{compute_prs}
\alias{compute_mprs}
\title{Trait-specific and multi-trait polygenic risk scores}
\usage{
compute_prs(dosage_matrix, effect_sizes, ancestry = NULL)

compute_mprs(prs_matrix, asthma_status, train_idx = NULL, alpha = 0.5)
}
\arguments{
\item{dosage_matrix}{Numeric matrix, subjects x variants, of imputed
genotype dosages (0-2).}
\item{effect_sizes}{Named numeric vector of posterior SNP effect sizes;
names must match \code{colnames(dosage_matrix)}.}
\item{ancestry}{Optional factor/character vector (length = nrow), used to
standardize the score within stratum for cross-ancestry cohorts.}
\item{prs_matrix}{Numeric matrix or data frame, subjects x trait-specific
PRS.}
\item{asthma_status}{Binary (0/1) outcome vector.}
\item{train_idx}{Integer/logical vector indexing training subjects for
weight estimation, to avoid information leakage into later evaluation.
Defaults to all subjects.}
\item{alpha}{Elastic-net mixing parameter. Default 0.5.}
}
\value{
\code{compute_prs()} returns a named numeric vector of per-subject
standardized trait-specific PRS. \code{compute_mprs()} returns a list with
\code{mprs} (named numeric vector, standardized), \code{weights} (named
numeric vector of trait weights), and \code{fit} (the underlying
\code{glmnet::cv.glmnet} fit).
}
\description{
\code{compute_prs()} computes
\eqn{PRS_{i,t} = \sum_j \beta_{j,t} G_{ij}}, the dosage-weighted sum of
posterior SNP effect sizes (typically obtained upstream via PRS-CSx or
SBayesRC). \code{compute_mprs()} combines multiple trait-specific PRS
columns into a single multi-trait polygenic risk score via elastic-net
regression against asthma status,
\eqn{MPRS_i = \sum_t w_t PRS_{i,t}}.
}
\examples{
dosage <- matrix(c(0,1,2,1, 2,0,1,1), nrow = 4,
                  dimnames = list(paste0("S", 1:4), c("rs1","rs2")))
eff <- c(rs1 = 0.12, rs2 = -0.08)
compute_prs(dosage, eff)
}
