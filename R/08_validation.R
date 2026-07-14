#' Compute discrimination, calibration, and reclassification metrics
#'
#' One-stop validation function producing the full metric set specified
#' in the TRACE-Asthma Methods (AUROC, AUPRC, calibration slope/intercept,
#' Hosmer-Lemeshow test) for a single fitted model, and NRI/IDI relative
#' to a reference (simpler, nested) model when one is supplied.
#'
#' @param observed Binary (0/1) observed outcome vector.
#' @param predicted Numeric vector of predicted probabilities from the
#'   model being evaluated.
#' @param predicted_reference Optional numeric vector of predicted
#'   probabilities from a simpler reference model, for NRI/IDI.
#' @param n_calibration_bins Number of bins for the calibration plot data
#'   and Hosmer-Lemeshow test. Default 10.
#' @return A list with `auroc`, `auprc`, `calibration` (list with `slope`,
#'   `intercept`, `hl_p`, and `plot_data` -- a data frame of
#'   observed-vs-predicted risk per bin), and, if `predicted_reference` is
#'   supplied, `nri` and `idi`.
#' @export
validate_model <- function(observed, predicted, predicted_reference = NULL,
                            n_calibration_bins = 10) {
  requirePkg("pROC", "AUROC computation")

  roc_obj <- pROC::roc(observed, predicted, quiet = TRUE)
  auroc <- as.numeric(pROC::auc(roc_obj))

  auprc <- if (requireNamespace("PRROC", quietly = TRUE)) {
    pr <- PRROC::pr.curve(scores.class0 = predicted[observed == 1],
                           scores.class1 = predicted[observed == 0], curve = FALSE)
    pr$auc.integral
  } else {
    warning("PRROC not installed; AUPRC not computed. install.packages('PRROC').")
    NA_real_
  }

  calibration <- .calibration_metrics(observed, predicted, n_bins = n_calibration_bins)

  result <- list(auroc = auroc, auprc = auprc, calibration = calibration)

  if (!is.null(predicted_reference)) {
    result$nri <- .compute_nri(observed, predicted, predicted_reference)
    result$idi <- .compute_idi(observed, predicted, predicted_reference)
  }

  structure(result, class = "traceAsthma_validation")
}

#' @keywords internal
.calibration_metrics <- function(observed, predicted, n_bins = 10) {
  eps <- 1e-6
  p <- pmin(pmax(predicted, eps), 1 - eps)
  logit_p <- log(p / (1 - p))

  cal_fit <- stats::glm(observed ~ logit_p, family = stats::binomial())
  slope <- stats::coef(cal_fit)[["logit_p"]]
  intercept <- stats::coef(cal_fit)[["(Intercept)"]]

  bins <- cut(p, breaks = stats::quantile(p, probs = seq(0, 1, length.out = n_bins + 1)),
              include.lowest = TRUE, labels = FALSE)
  plot_data <- do.call(rbind, lapply(sort(unique(bins)), function(b) {
    idx <- bins == b
    data.frame(bin = b, mean_predicted = mean(p[idx]), mean_observed = mean(observed[idx]),
               n = sum(idx))
  }))

  expected <- tapply(p, bins, sum)
  observed_sum <- tapply(observed, bins, sum)
  n_bin <- tapply(p, bins, length)
  hl_stat <- sum((observed_sum - expected)^2 / (expected * (1 - expected / n_bin)), na.rm = TRUE)
  hl_p <- stats::pchisq(hl_stat, df = length(unique(bins)) - 2, lower.tail = FALSE)

  list(slope = slope, intercept = intercept, hl_statistic = hl_stat, hl_p = hl_p,
       plot_data = plot_data)
}

#' @keywords internal
.compute_nri <- function(observed, predicted, predicted_reference) {
  up_events   <- mean(predicted[observed == 1] > predicted_reference[observed == 1]) -
                 mean(predicted[observed == 1] < predicted_reference[observed == 1])
  up_nonevents <- mean(predicted[observed == 0] < predicted_reference[observed == 0]) -
                  mean(predicted[observed == 0] > predicted_reference[observed == 0])
  up_events + up_nonevents
}

#' @keywords internal
.compute_idi <- function(observed, predicted, predicted_reference) {
  (mean(predicted[observed == 1]) - mean(predicted_reference[observed == 1])) -
    (mean(predicted[observed == 0]) - mean(predicted_reference[observed == 0]))
}

#' @export
print.traceAsthma_validation <- function(x, ...) {
  cat(sprintf("AUROC: %.3f\n", x$auroc))
  cat(sprintf("AUPRC: %.3f\n", x$auprc))
  cat(sprintf("Calibration slope: %.3f, intercept: %.3f, Hosmer-Lemeshow P = %.3f\n",
              x$calibration$slope, x$calibration$intercept, x$calibration$hl_p))
  if (!is.null(x$nri)) cat(sprintf("NRI: %.3f\n", x$nri))
  if (!is.null(x$idi)) cat(sprintf("IDI: %.3f\n", x$idi))
  invisible(x)
}

#' Decision curve analysis
#'
#' Thin wrapper around `dcurves::dca()` computing net benefit of the model
#' relative to treat-all / treat-none strategies across a range of risk
#' thresholds, for clinical-utility assessment as specified in the
#' TRACE-Asthma Methods.
#'
#' @param data Data frame with the outcome and predicted-probability
#'   columns.
#' @param outcome_col Character, name of the 0/1 outcome column.
#' @param model_cols Character vector of column name(s) holding predicted
#'   probabilities for one or more models to compare.
#' @param thresholds Numeric vector of risk thresholds to evaluate.
#'   Default `seq(0.01, 0.5, by = 0.01)`.
#' @return The `dcurves::dca()` result object (has its own `plot()` method).
#' @export
decision_curve <- function(data, outcome_col, model_cols,
                            thresholds = seq(0.01, 0.5, by = 0.01)) {
  requirePkg("dcurves", "decision curve analysis")
  form <- stats::as.formula(paste(outcome_col, "~", paste(model_cols, collapse = " + ")))
  dcurves::dca(form, data = data, thresholds = thresholds)
}
