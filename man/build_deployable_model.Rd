\name{build_deployable_model}
\alias{build_deployable_model}
\alias{print.traceAsthma_model}
\title{Construct a fitted TRACE-Asthma model object for deployment}
\usage{
build_deployable_model(eqtm_weights, tf_weights, regulons, risk_model,
  mprs_weights = NULL, score_center_scale = NULL, metadata = list())
}
\arguments{
\item{eqtm_weights}{Named numeric vector of CpG weights (from
\code{compute_eqtm_score()}$weights).}
\item{tf_weights}{Named numeric vector of TF weights (from
\code{compute_trace_score()}$weights).}
\item{regulons}{Regulon table used to compute TF activity at deployment
time, restricted to the TFs in \code{tf_weights}.}
\item{risk_model}{A fitted \code{glm} object from \code{fit_risk_model()}.}
\item{mprs_weights}{Optional named numeric vector of trait-PRS weights, if
the deployed model includes MPRS.}
\item{score_center_scale}{Optional named list with elements \code{eqtm},
\code{trace}, \code{mprs}, each a length-2 numeric vector \code{c(mean, sd)} from the
discovery cohort, used to standardize new patients' raw scores onto the
scale the risk model was trained on.}
\item{metadata}{Optional named list of free-form metadata (cohort name,
date frozen, contact) stored with the object for provenance/audit purposes.}
}
\value{
An object of class \code{traceAsthma_model}, suitable for \code{saveRDS()} and
later use with \code{predict_asthma_risk()}.
}
\description{
Bundles everything needed to score a new patient -- eQTM weights, TF
regulons + TRACE weights, and the final risk-model coefficients -- into a
single portable object that can be shipped to a clinical site without
needing to re-run discovery-stage analyses at the point of care.
}
\examples{
data(example_cohort)
fit <- fit_risk_model(
  data.frame(asthma = example_cohort$clinical$asthma, age = example_cohort$clinical$age),
  predictors = "age"
)
model <- build_deployable_model(
  eqtm_weights = c(cg_IL13 = 0.5),
  tf_weights = c(GATA3 = 0.5),
  regulons = example_cohort$regulon_db[example_cohort$regulon_db$tf == "GATA3", ],
  risk_model = fit
)
print(model)
}
