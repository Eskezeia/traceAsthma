\name{predict_asthma_risk}
\alias{predict_asthma_risk}
\title{Predict asthma risk for one or more new patients}
\usage{
predict_asthma_risk(model, methylation, expression, clinical = NULL,
  dosage = NULL, mprs_effect_sizes = NULL)
}
\arguments{
\item{model}{A \code{traceAsthma_model} object from
\code{build_deployable_model()}.}
\item{methylation}{Numeric matrix, CpGs x patients, containing at least
the CpGs in \code{model$eqtm_weights}.}
\item{expression}{Numeric matrix, genes x patients, containing at least
the target genes of the TFs in \code{model$tf_weights}.}
\item{clinical}{Optional data frame, one row per patient (same order as
the columns of \code{methylation}/\code{expression}), with any of
\code{age}, \code{sex}, \code{bmi}, \code{smoking} needed by
\code{model$risk_model}.}
\item{dosage}{Optional numeric matrix, patients x variants, for MPRS
computation if \code{model$mprs_weights} is set.}
\item{mprs_effect_sizes}{Optional named list of per-trait effect-size
vectors, required if \code{dosage} is supplied.}
}
\value{
A data frame, one row per patient, with columns \code{patient_id},
\code{eqtm_score}, \code{trace_score}, \code{mprs_score} (if applicable),
\code{predicted_probability}, and \code{risk_category} (a simple
three-tier Low/Moderate/High banding at 10\%/25\% predicted probability,
provided as a starting point only -- thresholds should be locally
validated and clinically calibrated before use in practice).
}
\description{
The primary clinician-facing entry point. Takes raw per-patient molecular
measurements and returns a predicted asthma probability together with the
intermediate eQTM / TF-activity / TRACE scores, so a clinician or genetic
counselor can see \emph{why} a patient was scored as high or low risk
rather than only the final number.
}
\section{Clinical use warning}{
This function returns a research-grade risk estimate. See
\code{vignette("clinical-deployment-notes", package = "traceAsthma")}
before using its output to inform any clinical decision.
}
\examples{
\dontrun{
predict_asthma_risk(
  model = deployable_model,
  methylation = new_patient_methylation,
  expression  = new_patient_expression,
  clinical    = data.frame(age = 34, sex = "F", bmi = 24.1, smoking = "never")
)
}
}
