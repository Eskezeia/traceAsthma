#' traceAsthma: A Mechanistic Multi-Omics Asthma Risk Model
#'
#' `traceAsthma` implements the TRACE-Asthma pipeline end to end: genotype
#' / methylation / RNA-seq QC, multi-trait polygenic scoring, cis/trans
#' eQTM discovery with mediation-based filtering, data-driven transcription
#' factor discovery and activity inference, TRACE score construction with
#' stability selection, structural equation modeling of the
#' MPRS -> eQTM -> TRACE -> Asthma mediation pathway, nested
#' cross-validated model training, external validation, and a
#' clinician-facing prediction interface.
#'
#' @section Getting started:
#' See `vignette("trace-asthma-pipeline", package = "traceAsthma")` for a
#' full worked example on simulated data, and `?run_trace_asthma_pipeline`
#' for the end-to-end orchestrator. For applying an already-fitted model
#' to a new patient, see `?predict_asthma_risk`.
#'
#' @section Important clinical disclaimer:
#' TRACE-Asthma is a research risk-modeling framework. Any model built
#' with this package requires site-specific validation, calibration, and
#' appropriate regulatory clearance before use in clinical decision-making.
#' See `vignette("clinical-deployment-notes")`.
#'
#' @keywords internal
"_PACKAGE"
