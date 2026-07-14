\name{traceAsthma-package}
\alias{traceAsthma-package}
\alias{traceAsthma}
\docType{package}
\title{traceAsthma: A Mechanistic Multi-Omics Asthma Risk Model}
\description{
Implements the TRACE-Asthma pipeline end to end: genotype / methylation /
RNA-seq QC, multi-trait polygenic scoring, cis/trans eQTM discovery with
mediation-based filtering, data-driven transcription factor discovery and
activity inference, TRACE score construction with stability selection,
structural equation modeling of the MPRS -> eQTM -> TRACE -> Asthma
mediation pathway, nested cross-validated model training, external
validation, and a clinician-facing prediction interface.
}
\details{
See \code{vignette("trace-asthma-pipeline", package = "traceAsthma")} for
a full worked example on simulated data, and
\code{\link{run_trace_asthma_pipeline}} for the end-to-end orchestrator.
For applying an already-fitted model to a new patient, see
\code{\link{predict_asthma_risk}}.

\strong{Important clinical disclaimer:} TRACE-Asthma is a research
risk-modeling framework. Any model built with this package requires
site-specific validation, calibration, and appropriate regulatory
clearance before use in clinical decision-making. See
\code{vignette("clinical-deployment-notes", package = "traceAsthma")}.
}
\keyword{internal}
