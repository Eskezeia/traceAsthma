#' Check for and load an optional (Suggested) package
#'
#' Many pipeline steps depend on specialized genomics packages (glmnet,
#' minfi, viper, lavaan, ...) that are declared in `Suggests` rather than
#' `Imports`, so that `traceAsthma` can be installed and its lightweight
#' functions (scoring, prediction) used without pulling in the full
#' upstream genomics stack. This helper checks availability and gives an
#' informative, actionable error if a required optional package is missing.
#'
#' @param pkg Character. Package name.
#' @param purpose Character. Short description of what the package is used
#'   for, inserted into the error message.
#' @param bioc Logical. If `TRUE`, install instructions reference
#'   `BiocManager::install()` instead of `install.packages()`.
#' @return Invisibly, `TRUE` if the package is available.
#' @keywords internal
#' @export
requirePkg <- function(pkg, purpose = "this step", bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install_cmd <- if (bioc) {
      sprintf('BiocManager::install("%s")', pkg)
    } else {
      sprintf('install.packages("%s")', pkg)
    }
    stop(sprintf(
      "Package '%s' is required for %s but is not installed.\nInstall it with:\n  %s",
      pkg, purpose, install_cmd
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Standardize a numeric vector to zero mean and unit variance
#'
#' Thin wrapper around [base::scale()] that always returns a plain numeric
#' vector (not a matrix with attributes), and optionally standardizes
#' within strata (e.g. within ancestry group or within sex) rather than
#' across the whole sample.
#'
#' @param x Numeric vector.
#' @param strata Optional factor/character vector of the same length as
#'   `x` giving the stratum for each observation. If supplied, `x` is
#'   centered and scaled separately within each stratum.
#' @return Numeric vector of standardized (Z-scored) values.
#' @export
#' @examples
#' z_score(c(1, 2, 3, 4, 5))
#' z_score(c(1, 2, 3, 10, 11, 12), strata = c("A","A","A","B","B","B"))
z_score <- function(x, strata = NULL) {
  if (is.null(strata)) {
    return(as.numeric(scale(x)[, 1]))
  }
  strata <- as.character(strata)
  out <- rep(NA_real_, length(x))
  for (s in unique(strata)) {
    idx <- strata == s
    out[idx] <- as.numeric(scale(x[idx])[, 1])
  }
  out
}

#' Compute a weighted linear score from a feature matrix and weight vector
#'
#' Generic engine behind the eQTM score, TRACE score, and MPRS: given a
#' matrix of per-subject feature values (methylation levels, TF activities,
#' trait-specific PRS, ...) and a named weight vector, returns the
#' per-subject weighted sum. Used internally by [compute_eqtm_score()],
#' [compute_trace_score()], and [compute_mprs()].
#'
#' @param feature_matrix Numeric matrix or data frame, subjects x features.
#'   Column names must match `names(weights)`.
#' @param weights Named numeric vector of weights, e.g. elastic-net
#'   coefficients.
#' @param standardize Logical; if `TRUE` (default) the raw weighted sum is
#'   Z-scored across subjects before being returned.
#' @return Named numeric vector of per-subject scores (raw or standardized).
#' @export
#' @examples
#' fm <- matrix(c(0.72, 0.35, 0.61, 0.42), nrow = 2,
#'              dimnames = list(c("P001","P002"), c("cg_IL13","cg_GATA3")))
#' w <- c(cg_IL13 = 0.55, cg_GATA3 = 0.45)
#' weighted_score(fm, w)
weighted_score <- function(feature_matrix, weights, standardize = TRUE) {
  feature_matrix <- as.matrix(feature_matrix)
  common <- intersect(colnames(feature_matrix), names(weights))
  if (length(common) == 0) {
    stop("No overlap between feature_matrix columns and names(weights).", call. = FALSE)
  }
  if (length(common) < length(weights)) {
    missing <- setdiff(names(weights), colnames(feature_matrix))
    warning(sprintf(
      "%d weighted feature(s) not found in feature_matrix and were dropped: %s",
      length(missing), paste(missing, collapse = ", ")
    ))
  }
  raw <- as.numeric(feature_matrix[, common, drop = FALSE] %*% weights[common])
  names(raw) <- rownames(feature_matrix)
  if (standardize) z_score(raw) else raw
}

#' Cluster and prune highly correlated features
#'
#' Implements the correlation-pruning step used before elastic-net feature
#' selection for both CpGs (eQTM layer) and transcription factors (TRACE
#' layer): features pairwise-correlated above `threshold` are grouped, and
#' only the feature with the strongest association with the outcome
#' (largest `|score|`) is retained from each group.
#'
#' @param feature_matrix Numeric matrix, subjects x features.
#' @param assoc_score Named numeric vector, one association statistic per
#'   feature (e.g. absolute t-statistic or -log10(P) from a univariate
#'   association test with the outcome). Names must match
#'   `colnames(feature_matrix)`.
#' @param threshold Numeric in (0, 1); absolute correlation above which two
#'   features are considered redundant. Default 0.90, matching the
#'   TRACE-Asthma specification.
#' @return Character vector of retained feature names.
#' @export
prune_correlated_features <- function(feature_matrix, assoc_score, threshold = 0.90) {
  feature_matrix <- as.matrix(feature_matrix)
  feats <- colnames(feature_matrix)
  if (is.null(feats)) stop("feature_matrix must have column names.", call. = FALSE)
  assoc_score <- assoc_score[feats]

  cor_mat <- suppressWarnings(stats::cor(feature_matrix, use = "pairwise.complete.obs"))
  cor_mat[is.na(cor_mat)] <- 0

  remaining <- feats
  retained <- character(0)

  # order by strength of association, strongest first, so the "winner"
  # of each correlated cluster is decided by outcome relevance
  ordered <- names(sort(abs(assoc_score), decreasing = TRUE))

  dropped <- character(0)
  for (f in ordered) {
    if (f %in% dropped) next
    retained <- c(retained, f)
    correlated_with_f <- feats[abs(cor_mat[f, ]) > threshold & feats != f]
    dropped <- union(dropped, correlated_with_f)
  }
  retained
}

#' Bootstrap stability selection for elastic-net-selected features
#'
#' Repeatedly re-fits an elastic-net model on bootstrap resamples of the
#' data and records, for each candidate feature, the proportion of
#' resamples in which it received a non-zero coefficient. Used for TF
#' selection in the TRACE layer, as specified in the TRACE-Asthma design
#' (1,000 resamples, retain features selected in > 80% of resamples), but
#' generic enough to reuse for the eQTM layer as well.
#'
#' @param x Numeric feature matrix, subjects x features.
#' @param y Binary outcome vector (0/1), same length as `nrow(x)`.
#' @param n_boot Integer, number of bootstrap resamples. Default 1000.
#' @param stability_threshold Numeric in (0, 1); minimum selection
#'   frequency required to retain a feature. Default 0.80.
#' @param alpha Elastic-net mixing parameter passed to `glmnet::cv.glmnet`.
#'   Default 0.5.
#' @param seed Optional integer seed for reproducibility.
#' @return A list with `selected` (character vector of stable feature
#'   names) and `frequency` (named numeric vector of per-feature selection
#'   frequencies across all resamples).
#' @export
stability_select <- function(x, y, n_boot = 1000, stability_threshold = 0.80,
                              alpha = 0.5, seed = NULL) {
  requirePkg("glmnet", "elastic-net feature selection / stability selection")
  x <- as.matrix(x)
  n <- nrow(x)
  feats <- colnames(x)
  if (!is.null(seed)) set.seed(seed)

  counts <- stats::setNames(rep(0L, length(feats)), feats)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    xb <- x[idx, , drop = FALSE]
    yb <- y[idx]
    if (length(unique(yb)) < 2) next  # skip degenerate resamples
    fit <- tryCatch(
      glmnet::cv.glmnet(xb, yb, family = "binomial", alpha = alpha, nfolds = 5),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    coefs <- as.matrix(stats::coef(fit, s = "lambda.min"))
    nz <- rownames(coefs)[abs(coefs[, 1]) > 0 & rownames(coefs) != "(Intercept)"]
    counts[nz] <- counts[nz] + 1L
  }

  freq <- counts / n_boot
  list(
    selected = names(freq)[freq > stability_threshold],
    frequency = freq
  )
}
