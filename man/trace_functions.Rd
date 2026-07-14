\name{compute_trace_score}
\alias{compute_trace_score}
\alias{summarize_trace_weights}
\title{TRACE score construction and weight summary}
\usage{
compute_trace_score(tf_activity, asthma_status, train_idx = NULL,
  fdr_threshold = 0.05, corr_threshold = 0.90, n_boot = 1000,
  stability_threshold = 0.80, alpha = 0.5, seed = NULL)

summarize_trace_weights(trace_fit, top_n = NULL)
}
\arguments{
\item{tf_activity}{Numeric matrix, TFs x subjects (output of \code{infer_tf_activity()}).}
\item{asthma_status}{Binary (0/1) outcome vector.}
\item{train_idx}{Integer/logical index of training subjects. Defaults to all.}
\item{fdr_threshold}{FDR threshold for the univariate screen. Default 0.05.}
\item{corr_threshold}{Correlation-pruning threshold. Default 0.90.}
\item{n_boot}{Number of stability-selection bootstrap resamples. Default 1000.}
\item{stability_threshold}{Minimum selection frequency to retain a TF. Default 0.80.}
\item{alpha}{Elastic-net mixing parameter. Default 0.5.}
\item{seed}{Optional integer seed.}
\item{trace_fit}{Output of \code{compute_trace_score()}.}
\item{top_n}{Number of top TFs to report. Default all.}
}
\value{
\code{compute_trace_score()} returns a list with \code{trace_score} (standardized
per-subject score), \code{weights}, and \code{stability_frequency}.
\code{summarize_trace_weights()} returns a data frame with columns \code{tf},
\code{weight}, \code{direction} ("activating"/"repressive"), sorted by \code{abs(weight)}.
}
\description{
Implements the full TF selection and TRACE scoring procedure: univariate
FDR screen, correlation pruning, bootstrap stability selection (retain
TFs selected in a high fraction of resamples), and final elastic-net
weighting on the stable TF set, producing
\eqn{TRACE_i = \sum_{t=1}^m \theta_t A_{it}}.
\code{summarize_trace_weights()} ranks contributing TFs for reporting in
figures/tables (see \code{plot_trace_weights()}).
}
\examples{
tf_activity <- matrix(rnorm(30), nrow = 3, dimnames = list(c("TF1","TF2","TF3"), paste0("P",1:10)))
asthma <- rbinom(10, 1, 0.4)
\donttest{
if (requireNamespace("glmnet", quietly = TRUE)) {
  fit <- compute_trace_score(tf_activity, asthma, fdr_threshold = 1, stability_threshold = 0)
  summarize_trace_weights(fit)
}
}
}
