#' Format the nested-model comparison as a publication table
#'
#' Produces a table matching Nature-family journal conventions: fully
#' spelled-out column headers (no bare "AUC" without units/definition in
#' a footnote), the specific predictors added at each step named directly
#' rather than "M1"-"M5", and 95% CIs / P-values formatted to a fixed,
#' consistent number of decimal places.
#'
#' @param nested_comparison Output of [compare_nested_models()].
#' @param auc_ci Optional data frame with columns `model`, `ci_low`,
#'   `ci_high` (see [plot_nested_auc()]).
#' @return A `traceAsthma_table` object (data frame with a `footnote`
#'   attribute); print it directly, or pass to [export_table()].
#' @export
format_nested_model_table <- function(nested_comparison, auc_ci = NULL) {
  cmp <- if (inherits(nested_comparison, "traceAsthma_nested")) nested_comparison$comparison else nested_comparison
  label_map <- c(
    M1 = "Asthma-only polygenic risk score",
    M2 = "Multi-trait polygenic risk score (MPRS)",
    M3 = "MPRS + methylation eQTM score",
    M4 = "MPRS + eQTM score + TRACE score",
    M5 = "MPRS + eQTM score + TRACE score + clinical covariates"
  )
  model_desc <- ifelse(cmp$model %in% names(label_map), label_map[cmp$model], cmp$predictors)

  auc_fmt <- if (!is.null(auc_ci)) {
    ord <- match(cmp$model, auc_ci$model)
    sprintf("%.3f (%.3f\u2013%.3f)", cmp$auc, auc_ci$ci_low[ord], auc_ci$ci_high[ord])
  } else {
    sprintf("%.3f", cmp$auc)
  }

  out <- data.frame(
    Model = model_desc,
    Predictors_included = cmp$predictors,
    AUC_95CI = auc_fmt,
    Delta_AUC = ifelse(is.na(cmp$delta_auc), "\u2014", sprintf("%+.3f", cmp$delta_auc)),
    LRT_P = ifelse(is.na(cmp$lrt_p), "\u2014", .format_p(cmp$lrt_p)),
    check.names = FALSE
  )
  colnames(out) <- c("Model", "Predictors included", "AUC (95% CI)",
                      "\u0394AUC vs. prior model", "Likelihood-ratio test P")

  structure(out, class = c("traceAsthma_table", "data.frame"),
            footnote = paste(
              "AUC, area under the receiver operating characteristic curve.",
              "\u0394AUC, change in AUC relative to the immediately preceding",
              "(one layer simpler) nested model. Likelihood-ratio test P-values",
              "compare each model to its immediate predecessor; the first model",
              "listed has no predecessor and is shown with \u2014."
            ))
}

#' Format a fitted risk model's coefficients as a publication table
#'
#' @param risk_model A fitted `glm` object (e.g. from [fit_risk_model()]).
#' @param variable_labels Optional named character vector mapping raw
#'   model term names (e.g. `"eqtm"`) to publication-ready labels (e.g.
#'   `"Methylation eQTM score (per SD)"`). Terms not in this vector are
#'   shown as-is; supplying this is strongly recommended so the exported
#'   table never shows a bare, non-descriptive variable name.
#' @return A `traceAsthma_table` object.
#' @export
format_coefficient_table <- function(risk_model, variable_labels = NULL) {
  smry <- summary(risk_model)$coefficients
  terms <- rownames(smry)
  labels <- if (!is.null(variable_labels)) {
    ifelse(terms %in% names(variable_labels), variable_labels[terms], terms)
  } else terms

  est <- smry[, "Estimate"]
  se  <- smry[, "Std. Error"]
  pv  <- smry[, "Pr(>|z|)"]
  or  <- exp(est)
  or_lo <- exp(est - 1.96 * se)
  or_hi <- exp(est + 1.96 * se)

  out <- data.frame(
    Predictor = labels,
    Beta_SE = sprintf("%.3f (%.3f)", est, se),
    OR_95CI = ifelse(terms == "(Intercept)", "\u2014",
                      sprintf("%.2f (%.2f\u2013%.2f)", or, or_lo, or_hi)),
    P_value = .format_p(pv),
    check.names = FALSE
  )
  colnames(out) <- c("Predictor", "\u03b2 (SE)", "Odds ratio (95% CI)", "P value")
  structure(out, class = c("traceAsthma_table", "data.frame"),
            footnote = paste(
              "\u03b2, logistic regression coefficient on the log-odds scale; SE,",
              "standard error. Odds ratios are per 1-SD increase for standardized",
              "continuous predictors (eQTM score, TRACE score, MPRS). Intercept",
              "is shown on the log-odds scale only, as an odds ratio for the",
              "intercept is not interpretable."
            ))
}

#' Format the mediation SEM effect decomposition as a publication table
#'
#' @param sem_result Output of [fit_mediation_sem()].
#' @return A `traceAsthma_table` object.
#' @export
format_mediation_table <- function(sem_result) {
  eff <- sem_result$effects
  label_map <- c(
    direct = "Direct effect (MPRS \u2192 Asthma)",
    indirect = "Indirect effect (MPRS \u2192 eQTM \u2192 TRACE \u2192 Asthma)",
    total = "Total effect (MPRS \u2192 Asthma)"
  )
  out <- data.frame(
    Effect = ifelse(eff$effect %in% names(label_map), label_map[eff$effect], eff$effect),
    Estimate = sprintf("%.3f", eff$estimate),
    CI_95 = sprintf("%.3f\u2013%.3f", eff$ci_low, eff$ci_high),
    P_value = .format_p(eff$p_value),
    check.names = FALSE
  )
  colnames(out) <- c("Effect", "Estimate", "95% CI", "P value")
  structure(out, class = c("traceAsthma_table", "data.frame"),
            footnote = sprintf(
              paste("Effects estimated by structural equation modeling with",
                    "bootstrap standard errors. Approximately %.1f%% of the total",
                    "effect of MPRS on asthma risk is statistically mediated",
                    "through the eQTM \u2192 TRACE pathway."),
              100 * sem_result$proportion_mediated
            ))
}

#' @keywords internal
.format_p <- function(p) {
  ifelse(is.na(p), "\u2014",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

#' @export
print.traceAsthma_table <- function(x, ...) {
  print.data.frame(x, row.names = FALSE, ...)
  fn <- attr(x, "footnote")
  if (!is.null(fn)) {
    cat("\n", strwrap(fn, width = 80, prefix = "  "), sep = "\n")
  }
  invisible(x)
}

#' Export a formatted table to CSV (with footnote) or Excel
#'
#' Writes the table as a CSV alongside a companion `_footnote.txt` file
#' (the default, dependency-free option, matching Nature's requirement
#' that source-data tables be plain, editable, and separated from
#' narrative formatting), or as a single-sheet `.xlsx` file with the
#' footnote written into the cell below the table if the `writexl`
#' package is installed.
#'
#' @param table A `traceAsthma_table` object (from `format_*_table()`).
#' @param path Output file path (`.csv` or `.xlsx`).
#' @return Invisibly, `path`.
#' @export
export_table <- function(table, path) {
  ext <- tolower(tools::file_ext(path))
  fn <- attr(table, "footnote")

  if (ext == "xlsx") {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      message("writexl not installed; falling back to CSV. install.packages('writexl') for native .xlsx export.")
      path <- sub("\\.xlsx$", ".csv", path)
      ext <- "csv"
    } else {
      out <- rbind(table, stats::setNames(as.list(rep("", ncol(table))), colnames(table)))
      out <- rbind(out, stats::setNames(c(fn, rep("", ncol(table) - 1)), colnames(table)))
      writexl::write_xlsx(list(Table = out), path)
      return(invisible(path))
    }
  }

  df <- as.data.frame(table)
  esc <- function(x) paste0('"', gsub('"', '""', x), '"')
  header <- paste(esc(colnames(df)), collapse = ",")
  rows <- apply(df, 1, function(r) paste(esc(r), collapse = ","))
  content <- enc2utf8(c(header, rows))
  con <- file(path, open = "wb")
  writeLines(content, con, useBytes = TRUE)
  close(con)

  if (!is.null(fn)) {
    fn_con <- file(sub("\\.csv$", "_footnote.txt", path), open = "wb")
    writeLines(enc2utf8(fn), fn_con, useBytes = TRUE)
    close(fn_con)
  }
  invisible(path)
}
