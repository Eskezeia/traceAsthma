#' Fit the integrated TRACE-Asthma logistic risk model
#'
#' \deqn{logit(P(Asthma)) = \beta_0 + \beta_1 MPRS + \beta_2 eQTM +
#'   \beta_3 TRACE + \beta_4 Age + \beta_5 Sex + \beta_6 BMI + \beta_7 Smoking}
#'
#' A thin, explicit wrapper around [stats::glm()] so the final model
#' formula, and every coefficient in it, is transparent and auditable --
#' appropriate for a model intended for clinical use. Any subset of
#' predictors may be supplied (e.g. omit `mprs`/`bmi`/`smoking` if
#' genotype or clinical data are unavailable at a given site; see the
#' package vignette for the reduced eQTM + TRACE + age + sex
#' specification).
#'
#' @param data Data frame containing the outcome column `asthma` (0/1) and
#'   any subset of predictor columns among `mprs`, `eqtm`, `trace`, `age`,
#'   `sex`, `bmi`, `smoking`.
#' @param predictors Character vector of column names in `data` to include
#'   as predictors. Default: all of `mprs`, `eqtm`, `trace`, `age`, `sex`,
#'   `bmi`, `smoking` that are present in `data`.
#' @return A fitted `glm` object (family = binomial).
#' @export
fit_risk_model <- function(data, predictors = NULL) {
  default_predictors <- c("mprs", "eqtm", "trace", "age", "sex", "bmi", "smoking")
  if (is.null(predictors)) predictors <- intersect(default_predictors, colnames(data))
  if (!("asthma" %in% colnames(data))) {
    stop("`data` must contain an `asthma` (0/1) outcome column.", call. = FALSE)
  }
  form <- stats::as.formula(paste("asthma ~", paste(predictors, collapse = " + ")))
  stats::glm(form, data = data, family = stats::binomial())
}

#' Fit and compare nested TRACE-Asthma models (M1-M5)
#'
#' Fits the five nested models specified in the TRACE-Asthma design --
#' M1 (asthma-only PRS), M2 (MPRS), M3 (MPRS + eQTM), M4 (MPRS + eQTM +
#' TRACE), M5 (M4 + clinical covariates) -- and compares each to its
#' immediate predecessor via likelihood-ratio test and AUC change, exactly
#' as required for the incremental-value claims in the manuscript.
#' Predictor names are adapted automatically to whichever of `mprs`/
#' `prs_asthma`/`eqtm`/`trace`/clinical columns are actually present in
#' `data`, so this also works with the reduced (no-genotype) specification.
#'
#' @param data Data frame with outcome column `asthma` and any of
#'   `prs_asthma`, `mprs`, `eqtm`, `trace`, `age`, `sex`, `bmi`, `smoking`.
#' @return A list of class `traceAsthma_nested` with elements `models`
#'   (named list of fitted `glm` objects), `comparison` (data frame with
#'   one row per model, columns `model`, `predictors`, `auc`,
#'   `delta_auc`, `lrt_p`).
#' @export
compare_nested_models <- function(data) {
  requirePkg("pROC", "AUC computation for nested model comparison")

  spec <- list(
    M1 = "prs_asthma",
    M2 = "mprs",
    M3 = c("mprs", "eqtm"),
    M4 = c("mprs", "eqtm", "trace"),
    M5 = c("mprs", "eqtm", "trace", "age", "sex", "bmi", "smoking")
  )
  spec <- lapply(spec, function(v) intersect(v, colnames(data)))
  spec <- Filter(function(v) length(v) > 0, spec)

  models <- list()
  rows <- list()
  prev_fit <- NULL
  prev_name <- NULL

  for (nm in names(spec)) {
    preds <- spec[[nm]]
    fit <- fit_risk_model(data, predictors = preds)
    models[[nm]] <- fit

    pred_prob <- stats::predict(fit, type = "response")
    roc_obj <- pROC::roc(data$asthma, pred_prob, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))

    lrt_p <- NA_real_
    delta_auc <- NA_real_
    if (!is.null(prev_fit)) {
      lrt <- stats::anova(prev_fit, fit, test = "LRT")
      lrt_p <- lrt[["Pr(>Chi)"]][2]
      prev_pred <- stats::predict(prev_fit, type = "response")
      prev_roc <- pROC::roc(data$asthma, prev_pred, quiet = TRUE)
      delta_auc <- auc_val - as.numeric(pROC::auc(prev_roc))
    }

    rows[[nm]] <- data.frame(
      model = nm, predictors = paste(preds, collapse = " + "),
      auc = auc_val, delta_auc = delta_auc, lrt_p = lrt_p
    )
    prev_fit <- fit
    prev_name <- nm
  }

  structure(
    list(models = models, comparison = do.call(rbind, rows)),
    class = "traceAsthma_nested"
  )
}

#' @export
print.traceAsthma_nested <- function(x, ...) {
  cat("TRACE-Asthma nested model comparison\n")
  print(x$comparison, row.names = FALSE)
  invisible(x)
}

#' Nested cross-validation for TRACE-Asthma model development
#'
#' Runs outer k-fold cross-validation; within each outer training fold,
#' re-executes the full feature-construction pipeline (correlation
#' pruning, elastic-net feature selection for eQTM, TF stability selection
#' for TRACE, elastic-net weighting for MPRS) via user-supplied
#' `build_fun`, so that no information from the held-out outer test fold
#' leaks into feature selection -- required for an unbiased performance
#' estimate. Inner-loop hyperparameter tuning is delegated to
#' `glmnet::cv.glmnet()` inside each `*_score()` function.
#'
#' @param subject_ids Character/integer vector of subject identifiers
#'   defining the cross-validation partition.
#' @param build_fun A function taking `train_idx` (indices into
#'   `subject_ids`) and returning a fitted model object exposing a
#'   `predict(newdata_idx)` interface via `predict_fun` (see below). In
#'   practice this will call [compute_mprs()], [compute_eqtm_score()],
#'   [compute_trace_score()], and [fit_risk_model()] in sequence using
#'   only `train_idx`, and return everything needed to score the held-out
#'   fold.
#' @param predict_fun A function taking the object returned by
#'   `build_fun` and a `test_idx`, returning a numeric vector of predicted
#'   probabilities for the held-out subjects.
#' @param outcome Binary (0/1) outcome vector, same length/order as
#'   `subject_ids`.
#' @param k Number of outer folds. Default 10.
#' @param seed Optional integer seed for fold assignment.
#' @return A list with `fold_auc` (numeric vector, length `k`),
#'   `mean_auc`, `sd_auc`, and `oof_predictions` (named numeric vector of
#'   out-of-fold predicted probabilities for every subject, suitable for a
#'   single pooled ROC/calibration curve).
#' @export
nested_cv <- function(subject_ids, build_fun, predict_fun, outcome, k = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(subject_ids)
  fold_id <- sample(rep(seq_len(k), length.out = n))

  fold_auc <- rep(NA_real_, k)
  oof_pred <- rep(NA_real_, n)
  names(oof_pred) <- subject_ids

  for (fold in seq_len(k)) {
    train_idx <- which(fold_id != fold)
    test_idx  <- which(fold_id == fold)

    model_obj <- build_fun(train_idx)
    pred <- predict_fun(model_obj, test_idx)
    oof_pred[test_idx] <- pred

    if (requireNamespace("pROC", quietly = TRUE) && length(unique(outcome[test_idx])) == 2) {
      roc_obj <- pROC::roc(outcome[test_idx], pred, quiet = TRUE)
      fold_auc[fold] <- as.numeric(pROC::auc(roc_obj))
    }
  }

  list(
    fold_auc = fold_auc,
    mean_auc = mean(fold_auc, na.rm = TRUE),
    sd_auc = stats::sd(fold_auc, na.rm = TRUE),
    oof_predictions = oof_pred
  )
}
