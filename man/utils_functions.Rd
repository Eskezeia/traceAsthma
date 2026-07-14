\name{requirePkg}
\alias{requirePkg}
\alias{weighted_score}
\alias{z_score}
\alias{prune_correlated_features}
\alias{stability_select}
\title{Shared utility functions used across the pipeline}
\usage{
requirePkg(pkg, purpose = "this step", bioc = FALSE)

weighted_score(feature_matrix, weights, standardize = TRUE)

z_score(x, strata = NULL)

prune_correlated_features(feature_matrix, assoc_score, threshold = 0.90)

stability_select(x, y, n_boot = 1000, stability_threshold = 0.80,
  alpha = 0.5, seed = NULL)
}
\arguments{
\item{pkg}{Character package name.}
\item{purpose}{Short description of what the package is used for, inserted
into the error message if missing.}
\item{bioc}{Logical; if \code{TRUE}, install instructions reference
\code{BiocManager::install()} instead of \code{install.packages()}.}
\item{feature_matrix}{Numeric matrix or data frame, subjects x features.}
\item{weights}{Named numeric vector of weights.}
\item{standardize}{Logical; if \code{TRUE} (default) the raw weighted sum is
Z-scored across subjects.}
\item{x}{Numeric vector (\code{z_score()}) or feature matrix (\code{stability_select()}).}
\item{strata}{Optional factor/character vector for within-stratum standardization.}
\item{assoc_score}{Named numeric vector, one association statistic per feature.}
\item{threshold}{Absolute correlation above which two features are
considered redundant. Default 0.90.}
\item{y}{Binary outcome vector (0/1) for \code{stability_select()}.}
\item{n_boot}{Number of bootstrap resamples. Default 1000.}
\item{stability_threshold}{Minimum selection frequency required to retain
a feature. Default 0.80.}
\item{alpha}{Elastic-net mixing parameter. Default 0.5.}
\item{seed}{Optional integer seed.}
}
\value{
\code{requirePkg()} invisibly returns \code{TRUE} or throws an informative error.
\code{weighted_score()} returns a named numeric vector of per-subject scores.
\code{z_score()} returns a standardized numeric vector.
\code{prune_correlated_features()} returns a character vector of retained
feature names. \code{stability_select()} returns a list with \code{selected} and
\code{frequency}.
}
\description{
Internal engines shared across the eQTM, TRACE, and MPRS scoring
functions, exported for advanced/custom use: checking and requesting
installation of optional dependencies (\code{requirePkg()}), computing a
generic weighted linear score (\code{weighted_score()}), Z-standardizing a
vector overall or within strata (\code{z_score()}), clustering and pruning
highly correlated features before feature selection
(\code{prune_correlated_features()}), and bootstrap stability selection for
elastic-net-selected features (\code{stability_select()}).
}
\examples{
z_score(c(1, 2, 3, 4, 5))

fm <- matrix(c(0.72, 0.35, 0.61, 0.42), nrow = 2,
             dimnames = list(c("P001","P002"), c("cg_IL13","cg_GATA3")))
w <- c(cg_IL13 = 0.55, cg_GATA3 = 0.45)
weighted_score(fm, w)
}
