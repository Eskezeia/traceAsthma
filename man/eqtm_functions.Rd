\name{discover_eqtm}
\alias{discover_eqtm}
\alias{filter_asthma_relevant}
\alias{test_mediation}
\alias{compute_eqtm_score}
\title{Cis/trans eQTM discovery, asthma-relevant filtering, mediation testing, and scoring}
\usage{
discover_eqtm(methylation, expression, pairs, covariates, same_chr = NULL,
  cis_window_bp = 1e6, fdr_threshold = 0.05)

filter_asthma_relevant(eqtm_table, methylation, expression, asthma_status,
  covariates, p_threshold = 0.05)

test_mediation(eqtm_table, methylation, expression, asthma_status,
  covariates, n_boot = 10000, p_threshold = 0.05, seed = NULL)

compute_eqtm_score(methylation, asthma_status, train_idx = NULL,
  corr_threshold = 0.90, alpha = 0.5)
}
\arguments{
\item{methylation}{Numeric matrix, CpGs x subjects (M-values recommended
for \code{discover_eqtm()}; the CpGs retained after mediation filtering for
\code{compute_eqtm_score()}).}
\item{expression}{Numeric matrix, genes x subjects (normalized).}
\item{pairs}{Data frame with columns \code{cpg}, \code{gene}, \code{distance_bp}.}
\item{covariates}{Numeric matrix or data frame, subjects x covariates
(age, sex, ancestry PCs, smoking, cell-type proportions, surrogate
variables, ...). Row order must match the sample order of
\code{methylation}/\code{expression}.}
\item{same_chr}{Logical vector, one per row of \code{pairs}; \code{TRUE} if CpG and
gene are on the same chromosome. Defaults to \code{TRUE} for all pairs.}
\item{cis_window_bp}{Cis/trans boundary in base pairs. Default 1e6.}
\item{fdr_threshold}{FDR threshold for retaining an eQTM. Default 0.05.}
\item{eqtm_table}{Output of \code{discover_eqtm()} (or \code{filter_asthma_relevant()}).}
\item{asthma_status}{Binary (0/1) outcome vector.}
\item{p_threshold}{P-value threshold for retention. Default 0.05.}
\item{n_boot}{Number of bootstrap resamples for mediation CIs. Default 10000.}
\item{seed}{Optional integer seed.}
\item{train_idx}{Integer/logical index of training subjects. Defaults to all.}
\item{corr_threshold}{Correlation-pruning threshold. Default 0.90.}
\item{alpha}{Elastic-net mixing parameter. Default 0.5.}
}
\value{
\code{discover_eqtm()} returns a data frame of significant CpG-gene pairs
with columns \code{cpg}, \code{gene}, \code{type} ("cis"/"trans"), \code{alpha1}, \code{se},
\code{p_value}, \code{fdr}. \code{filter_asthma_relevant()} returns the input table
filtered to pairs with independent CpG- and gene-level asthma association,
with added columns \code{cpg_asthma_p}, \code{gene_asthma_p}. \code{test_mediation()}
returns the input table filtered to pairs with significant CpG -> expression
-> asthma mediation, with added columns \code{indirect_effect}, \code{ci_low},
\code{ci_high}, \code{p_mediation}. \code{compute_eqtm_score()} returns a list with
\code{eqtm_score} (standardized per-subject score), \code{weights}, \code{selected_cpgs},
and \code{fit}.
}
\description{
The core eQTM discovery pipeline. \code{discover_eqtm()} fits
\eqn{Y_g = \alpha_0 + \alpha_1 M_c + C\gamma + \epsilon} for each candidate
CpG-gene pair and classifies pairs as cis or trans, applying
Benjamini-Hochberg correction separately within each class.
\code{filter_asthma_relevant()} requires the CpG and target gene to each show
independent asthma association. \code{test_mediation()} tests the formal
CpG -> expression -> asthma mediation pathway (indirect effect \eqn{a \times b})
with bootstrap confidence intervals. \code{compute_eqtm_score()} applies
correlation pruning followed by elastic-net feature selection to construct
the final patient-level eQTM score,
\eqn{eQTM_i = \sum_{j=1}^k w_j M_{ij}}.
}
\examples{
data(example_cohort)
eqtm <- discover_eqtm(example_cohort$methylation, example_cohort$expression,
                       example_cohort$eqtm_pairs,
                       data.frame(age = example_cohort$covariates$age),
                       same_chr = example_cohort$same_chr, fdr_threshold = 1.0)
eqtm
}
