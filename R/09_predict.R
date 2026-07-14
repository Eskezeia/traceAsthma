#' Construct a fitted TRACE-Asthma model object for deployment
#'
#' Bundles everything needed to score a new patient -- eQTM weights, TF
#' regulons + TRACE weights, and the final risk-model coefficients -- into
#' a single portable object that can be saved with [saveRDS()] and shipped
#' to a clinical site, without needing to re-run discovery-stage analyses
#' (mediation testing, stability selection, etc.) at the point of care.
#'
#' @param eqtm_weights Named numeric vector of CpG weights (from
#'   [compute_eqtm_score()]$weights).
#' @param tf_weights Named numeric vector of TF weights (from
#'   [compute_trace_score()]$weights).
#' @param regulons Regulon table used to compute TF activity at deployment
#'   time (output of [build_regulons()], restricted to the TFs in
#'   `tf_weights`).
#' @param risk_model A fitted `glm` object from [fit_risk_model()].
#' @param mprs_weights Optional named numeric vector of trait-PRS weights
#'   (from [compute_mprs()]$weights), if the deployed model includes MPRS.
#' @param score_center_scale Optional named list with elements `eqtm`,
#'   `trace`, `mprs`, each a length-2 numeric vector `c(mean, sd)` from the
#'   discovery cohort, used to standardize new patients' raw scores onto
#'   the same scale the risk model was trained on (critical for a
#'   single-patient prediction, where within-sample standardization is not
#'   possible).
#' @param metadata Optional named list of free-form metadata (cohort name,
#'   date frozen, package version used to train, contact) stored with the
#'   object for provenance / audit purposes.
#' @return An object of class `traceAsthma_model`.
#' @export
build_deployable_model <- function(eqtm_weights, tf_weights, regulons, risk_model,
                                    mprs_weights = NULL, score_center_scale = NULL,
                                    metadata = list()) {
  metadata$traceAsthma_version <- tryCatch(
    as.character(utils::packageVersion("traceAsthma")), error = function(e) "dev"
  )
  metadata$created <- as.character(Sys.time())

  structure(
    list(
      eqtm_weights = eqtm_weights,
      tf_weights = tf_weights,
      regulons = regulons,
      mprs_weights = mprs_weights,
      risk_model = risk_model,
      score_center_scale = score_center_scale,
      metadata = metadata
    ),
    class = "traceAsthma_model"
  )
}

#' @export
print.traceAsthma_model <- function(x, ...) {
  cat("<traceAsthma_model>\n")
  cat(" eQTM CpGs      :", length(x$eqtm_weights), "\n")
  cat(" TRACE TFs      :", length(x$tf_weights), "\n")
  cat(" MPRS traits    :", if (is.null(x$mprs_weights)) 0 else length(x$mprs_weights), "\n")
  cat(" Risk model terms:", paste(names(stats::coef(x$risk_model))[-1], collapse = ", "), "\n")
  cat(" Created        :", x$metadata$created %||% "unknown", "\n")
  cat(" Package version:", x$metadata$traceAsthma_version %||% "unknown", "\n")
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Predict asthma risk for one or more new patients
#'
#' The primary clinician-facing entry point. Takes raw per-patient
#' molecular measurements (methylation at the model's CpGs, expression at
#' the genes needed for TF activity inference, and optionally genotype
#' dosages / clinical covariates) and returns a predicted asthma
#' probability together with the intermediate eQTM / TF-activity / TRACE
#' scores, so a clinician or genetic counselor can see *why* a patient was
#' scored as high or low risk rather than only the final number.
#'
#' @param model A `traceAsthma_model` object from
#'   [build_deployable_model()].
#' @param methylation Numeric matrix, CpGs x patients, containing at least
#'   the CpGs in `model$eqtm_weights`.
#' @param expression Numeric matrix, genes x patients, containing at least
#'   the target genes of the TFs in `model$tf_weights` (per
#'   `model$regulons`).
#' @param clinical Optional data frame, one row per patient (same order as
#'   the columns of `methylation`/`expression`), with any of `age`, `sex`,
#'   `bmi`, `smoking` needed by `model$risk_model`.
#' @param dosage Optional numeric matrix, patients x variants, for MPRS
#'   computation if `model$mprs_weights` is set.
#' @param mprs_effect_sizes Optional named list of per-trait effect-size
#'   vectors, required if `dosage` is supplied (see [compute_prs()]).
#' @return A data frame, one row per patient, with columns `patient_id`,
#'   `eqtm_score`, `trace_score`, `mprs_score` (if applicable),
#'   `predicted_probability`, and `risk_category` (a simple three-tier
#'   Low/Moderate/High banding at 10%/25% predicted probability, provided
#'   as a starting point only -- **thresholds should be locally validated
#'   and clinically calibrated before use in practice; see Limitations in
#'   the package vignette**).
#' @export
predict_asthma_risk <- function(model, methylation, expression, clinical = NULL,
                                 dosage = NULL, mprs_effect_sizes = NULL) {
  if (!inherits(model, "traceAsthma_model")) {
    stop("`model` must be a traceAsthma_model object; see build_deployable_model().", call. = FALSE)
  }

  patient_ids <- colnames(expression)
  if (is.null(patient_ids)) patient_ids <- colnames(methylation)
  if (is.null(patient_ids)) stop("methylation/expression matrices must have column names (patient IDs).", call. = FALSE)
  n <- length(patient_ids)

  # --- eQTM score ---
  missing_cpgs <- setdiff(names(model$eqtm_weights), rownames(methylation))
  if (length(missing_cpgs) > 0) {
    warning(sprintf("%d of %d model CpGs not found in supplied methylation data: %s%s",
                     length(missing_cpgs), length(model$eqtm_weights),
                     paste(utils::head(missing_cpgs, 5), collapse = ", "),
                     if (length(missing_cpgs) > 5) ", ..." else ""))
  }
  eqtm_raw <- weighted_score(t(methylation), model$eqtm_weights, standardize = FALSE)
  eqtm_score <- .apply_center_scale(eqtm_raw, model$score_center_scale$eqtm)

  # --- TF activity + TRACE score ---
  regulons_used <- model$regulons[model$regulons$tf %in% names(model$tf_weights), ]
  tf_activity <- infer_tf_activity(expression, regulons_used)
  trace_raw <- weighted_score(t(tf_activity), model$tf_weights, standardize = FALSE)
  trace_score <- .apply_center_scale(trace_raw, model$score_center_scale$trace)

  out <- data.frame(patient_id = patient_ids, eqtm_score = eqtm_score, trace_score = trace_score)

  # --- MPRS (optional) ---
  if (!is.null(model$mprs_weights)) {
    if (is.null(dosage) || is.null(mprs_effect_sizes)) {
      stop("model includes MPRS but `dosage` / `mprs_effect_sizes` were not supplied.", call. = FALSE)
    }
    trait_prs <- sapply(names(mprs_effect_sizes), function(trait) {
      compute_prs(dosage, mprs_effect_sizes[[trait]])
    })
    mprs_raw <- weighted_score(trait_prs, model$mprs_weights, standardize = FALSE)
    out$mprs_score <- .apply_center_scale(mprs_raw, model$score_center_scale$mprs)
  }

  # --- clinical covariates ---
  if (!is.null(clinical)) {
    stopifnot(nrow(clinical) == n)
    out <- cbind(out, clinical)
  }

  # --- final prediction ---
  newdata <- out
  pred <- stats::predict(model$risk_model, newdata = newdata, type = "response")
  out$predicted_probability <- as.numeric(pred)
  out$risk_category <- cut(
    out$predicted_probability,
    breaks = c(-Inf, 0.10, 0.25, Inf),
    labels = c("Low", "Moderate", "High")
  )

  out
}

#' @keywords internal
.apply_center_scale <- function(raw, center_scale) {
  if (is.null(center_scale)) return(raw)  # no discovery-cohort reference available; return raw
  (raw - center_scale[1]) / center_scale[2]
}
