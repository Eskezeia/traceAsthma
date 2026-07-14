#' Fit the mechanistic mediation pathway MPRS -> eQTM -> TRACE -> Asthma
#'
#' Fits the structural equation model
#' \deqn{eQTM = a_0 + a_1 MPRS + \epsilon_1}
#' \deqn{TRACE = b_0 + b_1 eQTM + \epsilon_2}
#' \deqn{Asthma = c_0 + c_1 MPRS + c_2 eQTM + c_3 TRACE + \epsilon_3}
#' using `lavaan`, and decomposes the total effect of MPRS on asthma risk
#' into a direct effect (`c1`) and an indirect effect mediated through
#' eQTM and TRACE (`a1 * b1 * c3`). Bootstrap confidence intervals for
#' direct, indirect, and total effects are obtained via `lavaan`'s
#' bootstrap facility.
#'
#' @param data Data frame with columns `mprs`, `eqtm`, `trace`, `asthma`
#'   (all subjects; standardized scores recommended for interpretable
#'   path coefficients).
#' @param n_boot Number of bootstrap resamples for CI estimation. Default 10000.
#' @param seed Optional integer seed.
#' @return A list with `fit` (the `lavaan` fit object), `effects` (data
#'   frame with rows `direct`, `indirect`, `total`, columns `estimate`,
#'   `ci_low`, `ci_high`, `p_value`), and `proportion_mediated` (indirect
#'   effect as a fraction of the total effect).
#' @export
fit_mediation_sem <- function(data, n_boot = 10000, seed = NULL) {
  requirePkg("lavaan", "structural equation modeling of the mediation pathway")
  if (!is.null(seed)) set.seed(seed)

  required_cols <- c("mprs", "eqtm", "trace", "asthma")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("`data` is missing required column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  model <- '
    eqtm  ~ a1 * mprs
    trace ~ b1 * eqtm
    asthma ~ c1 * mprs + c2 * eqtm + c3 * trace

    indirect := a1 * b1 * c3
    direct   := c1
    total    := direct + indirect
  '

  fit <- lavaan::sem(model, data = data, se = "bootstrap", bootstrap = n_boot,
                      ordered = "asthma")

  pe <- lavaan::parameterEstimates(fit, boot.ci.type = "bca.simple", level = 0.95)
  eff <- pe[pe$label %in% c("direct", "indirect", "total"), ]

  effects <- data.frame(
    effect = eff$label,
    estimate = eff$est,
    ci_low = eff$ci.lower,
    ci_high = eff$ci.upper,
    p_value = eff$pvalue,
    row.names = NULL
  )

  total_est <- effects$estimate[effects$effect == "total"]
  indirect_est <- effects$estimate[effects$effect == "indirect"]
  prop_mediated <- if (length(total_est) && total_est != 0) indirect_est / total_est else NA_real_

  list(fit = fit, effects = effects, proportion_mediated = prop_mediated)
}

#' Print a plain-language summary of a fitted mediation SEM
#'
#' @param sem_result Output of [fit_mediation_sem()].
#' @return Invisibly returns `sem_result`; called for its printed side effect.
#' @export
summarize_mediation <- function(sem_result) {
  eff <- sem_result$effects
  fmt <- function(row) {
    sprintf("%s effect: %.3f (95%% CI %.3f to %.3f, P = %.4g)",
            tools::toTitleCase(row["effect"]),
            as.numeric(row["estimate"]), as.numeric(row["ci_low"]),
            as.numeric(row["ci_high"]), as.numeric(row["p_value"]))
  }
  for (i in seq_len(nrow(eff))) cat(fmt(eff[i, ]), "\n")
  cat(sprintf(
    "\nApproximately %.1f%% of the total MPRS effect on asthma risk is mediated through the eQTM -> TRACE pathway.\n",
    100 * sem_result$proportion_mediated
  ))
  invisible(sem_result)
}
