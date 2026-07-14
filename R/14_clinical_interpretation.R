#' Translate a fitted risk model into plain-language clinical relevance
#'
#' Statistical output (odds ratios, AUC, mediated proportion) answers "is
#' this associated" and "does this discriminate" -- this function answers
#' the question a clinician actually asks: "how much does this change the
#' estimated risk for a patient, can I trust that number, and is there a
#' modifiable step in the pathway." It computes absolute predicted risk
#' at low/median/high values of a chosen score (holding other covariates
#' at their sample median/mode), attaches a calibration statement if a
#' [validate_model()] result is supplied, and attaches a mechanistic
#' mediation statement if a [fit_mediation_sem()] result is supplied,
#' then assembles all of it into one narrative paragraph alongside the
#' underlying table.
#'
#' @param risk_model A fitted `glm` object (e.g. from [fit_risk_model()]).
#' @param data The data frame used to fit `risk_model` (or an equivalent
#'   reference sample), used to determine percentile values of `score_var`
#'   and median/mode values at which to hold other covariates fixed.
#' @param score_var Character, the model term to profile across
#'   percentiles (e.g. `"trace"`, `"eqtm"`, `"mprs"`). If `NULL`
#'   (default), the first of `"trace"`, `"eqtm"`, `"mprs"` found among
#'   `risk_model`'s terms is used.
#' @param percentiles Numeric vector of percentiles (0-1) at which to
#'   profile `score_var`. Default `c(0.10, 0.50, 0.90)` (low/median/high).
#' @param validation Optional output of [validate_model()]; if supplied,
#'   its calibration slope/intercept and Hosmer-Lemeshow P are translated
#'   into a plain-language trustworthiness statement.
#' @param sem_result Optional output of [fit_mediation_sem()]; if
#'   supplied, the proportion-mediated estimate is translated into a
#'   plain-language mechanistic statement naming the modifiable pathway.
#' @param outcome_label Character, used in the narrative in place of the
#'   generic word "outcome" (e.g. `"asthma"`, the default).
#' @return A `traceAsthma_interpretation` object: a list with
#'   `absolute_risk` (data frame: percentile label, `score_var` value,
#'   predicted probability, and, for all rows after the first, the
#'   absolute risk difference from the lowest profile),
#'   `calibration_statement` (character or `NULL`),
#'   `mediation_statement` (character or `NULL`), and `narrative`
#'   (character, the full assembled paragraph). Printing the object shows
#'   the narrative followed by the table.
#' @export
#' @examples
#' set.seed(1)
#' n <- 60
#' d <- data.frame(asthma = rbinom(n, 1, 0.4), age = round(rnorm(n, 40, 10)),
#'                  sex = rbinom(n, 1, 0.5), trace = rnorm(n))
#' fit <- fit_risk_model(d, predictors = c("age", "sex", "trace"))
#' interp <- interpret_clinical_relevance(fit, d)
#' print(interp)
interpret_clinical_relevance <- function(risk_model, data, score_var = NULL,
                                          percentiles = c(0.10, 0.50, 0.90),
                                          validation = NULL, sem_result = NULL,
                                          outcome_label = "asthma") {
  terms <- names(stats::coef(risk_model))[-1]
  if (is.null(score_var)) {
    priority <- c("trace", "eqtm", "mprs")
    hit <- priority[priority %in% terms]
    if (length(hit) == 0) {
      stop("Could not auto-detect a score variable among risk_model's terms (",
           paste(terms, collapse = ", "), "). Supply `score_var` explicitly.",
           call. = FALSE)
    }
    score_var <- hit[1]
  }
  if (!(score_var %in% colnames(data))) {
    stop("`score_var` ('", score_var, "') not found in `data`.", call. = FALSE)
  }

  other_vars <- setdiff(terms, score_var)
  baseline <- lapply(other_vars, function(v) {
    x <- data[[v]]
    if (all(x %in% c(0, 1), na.rm = TRUE)) {
      round(stats::median(x, na.rm = TRUE))
    } else {
      stats::median(x, na.rm = TRUE)
    }
  })
  names(baseline) <- other_vars

  score_vals <- stats::quantile(data[[score_var]], probs = percentiles, na.rm = TRUE)
  pct_labels <- .percentile_labels(percentiles)

  newdata <- do.call(rbind, lapply(seq_along(score_vals), function(i) {
    row <- as.data.frame(baseline, stringsAsFactors = FALSE)
    row[[score_var]] <- unname(score_vals[i])
    row
  }))
  pred <- stats::predict(risk_model, newdata = newdata, type = "response")

  absolute_risk <- data.frame(
    Profile = pct_labels,
    Score_value = round(unname(score_vals), 3),
    Predicted_probability_pct = round(100 * pred, 1),
    Absolute_difference_vs_low_pp = round(100 * (pred - pred[1]), 1)
  )
  colnames(absolute_risk) <- c(
    "Profile", paste0(score_var, " value"),
    paste0("Predicted probability of ", outcome_label, " (%)"),
    "Absolute difference vs. lowest profile (pp)"
  )

  calibration_statement <- NULL
  if (!is.null(validation)) {
    slope <- validation$calibration$slope
    hl_p <- validation$calibration$hl_p
    well_calibrated <- !is.na(slope) && abs(slope - 1) < 0.2 && (!is.na(hl_p) && hl_p > 0.05)
    calibration_statement <- if (well_calibrated) {
      sprintf(paste("Predicted probabilities are well calibrated (calibration slope = %.2f;",
                     "Hosmer-Lemeshow P = %.2f), meaning the percentages above can be read",
                     "as face-value probabilities rather than a relative score."),
              slope, hl_p)
    } else {
      sprintf(paste("Predicted probabilities show evidence of miscalibration (calibration",
                     "slope = %.2f; Hosmer-Lemeshow P = %.3f); the relative ordering of",
                     "patients by risk is likely still informative, but the absolute",
                     "percentages above should be interpreted cautiously and the model",
                     "should be recalibrated on local data before clinical use."),
              slope, hl_p)
    }
  }

  mediation_statement <- NULL
  if (!is.null(sem_result)) {
    pm <- sem_result$proportion_mediated
    mediation_statement <- sprintf(
      paste("Approximately %.0f%% of the total estimated genetic effect on %s risk operates",
            "through the measured methylation and transcriptional regulatory pathway",
            "(eQTM -> TRACE) rather than through fixed genetic sequence alone, identifying",
            "a potentially modifiable component of an otherwise fixed inherited risk."),
      100 * pm, outcome_label
    )
  }

  risk_range <- absolute_risk[[3]]
  narrative <- sprintf(
    paste("Moving from the %s to the %s of the %s distribution is associated with an estimated",
          "predicted %s probability of %.1f%% versus %.1f%%, holding other model covariates at",
          "their sample median/typical value -- an absolute difference of %.1f percentage points",
          "across the observed range of this score."),
    tolower(pct_labels[1]), tolower(pct_labels[length(pct_labels)]), score_var,
    outcome_label, risk_range[1], risk_range[length(risk_range)],
    risk_range[length(risk_range)] - risk_range[1]
  )
  if (!is.null(calibration_statement)) narrative <- paste(narrative, calibration_statement)
  if (!is.null(mediation_statement)) narrative <- paste(narrative, mediation_statement)

  structure(
    list(
      score_var = score_var,
      absolute_risk = absolute_risk,
      calibration_statement = calibration_statement,
      mediation_statement = mediation_statement,
      narrative = narrative
    ),
    class = "traceAsthma_interpretation"
  )
}

#' @keywords internal
.percentile_labels <- function(percentiles) {
  n <- length(percentiles)
  if (n == 3 && isTRUE(all.equal(percentiles, c(0.10, 0.50, 0.90)))) {
    return(c("10th percentile (low)", "50th percentile (median)", "90th percentile (high)"))
  }
  sprintf("%gth percentile", 100 * percentiles)
}

#' @export
print.traceAsthma_interpretation <- function(x, ...) {
  cat(strwrap(x$narrative, width = 80), sep = "\n")
  cat("\n\n")
  print.data.frame(x$absolute_risk, row.names = FALSE)
  invisible(x)
}
