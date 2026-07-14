#' Construct the TRACE score from inferred TF activities
#'
#' Implements the full TRACE-Asthma TF selection and scoring procedure:
#' (1) univariate association of each TF's activity with asthma status,
#' FDR-corrected; (2) correlation pruning at |r| > 0.90; (3) bootstrap
#' stability selection (default 1000 resamples, retain TFs selected in
#' > 80% of resamples); (4) final elastic-net weighting on the stability-
#' selected TF set. Produces
#' \deqn{TRACE_i = \sum_{t=1}^m \theta_t A_{it}}
#'
#' @param tf_activity Numeric matrix, TFs x subjects (output of
#'   [infer_tf_activity()]).
#' @param asthma_status Binary (0/1) outcome vector.
#' @param train_idx Integer/logical index of training subjects. Defaults to all.
#' @param fdr_threshold FDR threshold for the univariate screen. Default 0.05.
#' @param corr_threshold Correlation-pruning threshold. Default 0.90.
#' @param n_boot Number of stability-selection bootstrap resamples. Default 1000.
#' @param stability_threshold Minimum selection frequency to retain a TF.
#'   Default 0.80.
#' @param alpha Elastic-net mixing parameter. Default 0.5.
#' @param seed Optional integer seed.
#' @return A list with `trace_score` (named numeric vector, all subjects,
#'   standardized), `weights` (named numeric vector, selected TFs), and
#'   `stability_frequency` (named numeric vector, all TFs entering
#'   stability selection).
#' @export
compute_trace_score <- function(tf_activity, asthma_status, train_idx = NULL,
                                 fdr_threshold = 0.05, corr_threshold = 0.90,
                                 n_boot = 1000, stability_threshold = 0.80,
                                 alpha = 0.5, seed = NULL) {
  requirePkg("glmnet", "TRACE score elastic-net weighting")
  if (is.null(train_idx)) train_idx <- seq_len(ncol(tf_activity))

  feature_matrix <- t(tf_activity)  # subjects x TFs
  train_x <- feature_matrix[train_idx, , drop = FALSE]
  train_y <- asthma_status[train_idx]

  # Step 1: univariate association + FDR
  univ <- apply(train_x, 2, function(col) {
    fit <- tryCatch(stats::glm(train_y ~ col, family = stats::binomial()), error = function(e) NULL)
    if (is.null(fit)) return(c(z = 0, p = 1))
    smry <- summary(fit)$coefficients
    if (!("col" %in% rownames(smry))) return(c(z = 0, p = 1))
    c(z = smry["col", "z value"], p = smry["col", "Pr(>|z|)"])
  })
  fdr <- stats::p.adjust(univ["p", ], method = "BH")
  screened <- names(fdr)[fdr <= fdr_threshold]
  if (length(screened) == 0) {
    stop("No TFs survived the univariate FDR screen at threshold ", fdr_threshold, call. = FALSE)
  }

  # Step 2: correlation pruning
  assoc_score <- univ["z", screened]
  retained <- prune_correlated_features(train_x[, screened, drop = FALSE], assoc_score,
                                         threshold = corr_threshold)

  # Step 3: bootstrap stability selection
  stab <- stability_select(train_x[, retained, drop = FALSE], train_y,
                            n_boot = n_boot, stability_threshold = stability_threshold,
                            alpha = alpha, seed = seed)
  stable_tfs <- stab$selected
  if (length(stable_tfs) == 0) {
    stop("No TFs survived stability selection at threshold ", stability_threshold,
         ". Consider lowering `stability_threshold` or `n_boot`.", call. = FALSE)
  }

  # Step 4: final elastic-net weighting on the stable TF set
  fit <- glmnet::cv.glmnet(train_x[, stable_tfs, drop = FALSE], train_y,
                            family = "binomial", alpha = alpha)
  coefs <- as.matrix(stats::coef(fit, s = "lambda.min"))
  weights <- coefs[rownames(coefs) != "(Intercept)", 1]
  weights <- weights[weights != 0]

  if (length(weights) == 0) {
    stop("Final elastic-net weighting shrank all stability-selected TFs to zero.", call. = FALSE)
  }

  trace_score <- weighted_score(feature_matrix, weights, standardize = TRUE)

  list(
    trace_score = trace_score,
    weights = weights,
    stability_frequency = stab$frequency
  )
}

#' Summarize the biological interpretation of a fitted TRACE score
#'
#' Small convenience function that ranks the TFs contributing to a fitted
#' TRACE score by absolute weight, for reporting in Results tables /
#' figures (e.g. Table X in the manuscript template).
#'
#' @param trace_fit Output of [compute_trace_score()].
#' @param top_n Number of top TFs to report. Default all.
#' @return Data frame with columns `tf`, `weight`, `direction`
#'   (`"activating"` / `"repressive"` of asthma risk), sorted by
#'   `abs(weight)` descending.
#' @export
summarize_trace_weights <- function(trace_fit, top_n = NULL) {
  w <- trace_fit$weights
  df <- data.frame(
    tf = names(w),
    weight = as.numeric(w),
    direction = ifelse(w > 0, "activating", "repressive")
  )
  df <- df[order(-abs(df$weight)), ]
  if (!is.null(top_n)) df <- utils::head(df, top_n)
  rownames(df) <- NULL
  df
}
