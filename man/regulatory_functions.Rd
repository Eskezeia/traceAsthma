\name{build_regulons}
\alias{build_regulons}
\alias{test_regulon_enrichment}
\alias{infer_tf_activity}
\title{TF regulon construction, enrichment testing, and activity inference}
\usage{
build_regulons(regulon_db, confidence_classes = c("A", "B", "C"),
  target_universe = NULL)

test_regulon_enrichment(regulons, target_genes, background_genes,
  fdr_threshold = 0.05)

infer_tf_activity(expression, regulons, method = c("viper", "decoupleR"),
  min_targets = 5)
}
\arguments{
\item{regulon_db}{Data frame with columns \code{tf}, \code{target}, \code{confidence}
(values in "A","B","C","D","E"), and \code{mor} (mode of regulation, +1/-1).}
\item{confidence_classes}{Character vector of confidence classes to retain.
Default \code{c("A","B","C")}.}
\item{target_universe}{Optional character vector restricting \code{target} to a
specific gene set (e.g. trans-eQTM targets).}
\item{regulons}{Output of \code{build_regulons()}.}
\item{target_genes}{Character vector of genes to test for regulon enrichment.}
\item{background_genes}{Character vector of all genes in the tested universe.}
\item{fdr_threshold}{FDR threshold for retaining a candidate TF. Default 0.05.}
\item{expression}{Numeric matrix, genes x subjects (normalized).}
\item{method}{Either \code{"viper"} or \code{"decoupleR"}. Default \code{"viper"}
if installed, otherwise automatically falls back to \code{"decoupleR"}
(via \code{run_ulm()}, which does not require the \code{viper} package).}
\item{min_targets}{Minimum regulon size to retain a TF. Default 5.}
}
\value{
\code{build_regulons()} returns a filtered regulon data frame.
\code{test_regulon_enrichment()} returns a data frame with columns \code{tf},
\code{n_overlap}, \code{regulon_size}, \code{p_value}, \code{fdr}, filtered to
significant TFs. \code{infer_tf_activity()} returns a numeric matrix, TFs x
subjects, of inferred activity scores (normalized enrichment scores).
}
\description{
\code{build_regulons()} filters a curated TF-target database to
high-confidence interactions, optionally restricted to a target gene
universe (e.g. trans-eQTM targets, for data-driven TF discovery).
\code{test_regulon_enrichment()} implements that data-driven discovery step:
tests each candidate TF's regulon for hypergeometric enrichment among a
target gene set (e.g. significant trans-eQTM targets), so candidate TFs
are nominated from the study's own data rather than picked a priori from
prior literature. \code{infer_tf_activity()} estimates per-subject TF
activity, \eqn{A_t = NES_t}, from the coordinated enrichment of a TF's
regulon among differentially expressed genes, using VIPER (\code{viper}
package) or decoupleR's \code{run_ulm()} method (which is genuinely
independent of the \code{viper} package, unlike \code{decoupleR::run_viper()}),
rather than the TF's own transcript abundance.
}
\examples{
data(example_cohort)
regs <- build_regulons(example_cohort$regulon_db, confidence_classes = c("A","B","C"))
test_regulon_enrichment(regs, target_genes = c("IL4","IL5","IL13"),
                         background_genes = rownames(example_cohort$expression),
                         fdr_threshold = 1.0)
}
