\name{validate_model}
\alias{validate_model}
\alias{print.traceAsthma_validation}
\alias{decision_curve}
\title{Model validation: discrimination, calibration, reclassification, and decision curves}
\usage{
validate_model(observed, predicted, predicted_reference = NULL,
  n_calibration_bins = 10)

decision_curve(data, outcome_col, model_cols, thresholds = seq(0.01, 0.5, by = 0.01))
}
\arguments{
\item{observed}{Binary (0/1) observed outcome vector.}
\item{predicted}{Numeric vector of predicted probabilities from the model
being evaluated.}
\item{predicted_reference}{Optional numeric vector of predicted
probabilities from a simpler reference model, for NRI/IDI.}
\item{n_calibration_bins}{Number of bins for the calibration plot data and
Hosmer-Lemeshow test. Default 10.}
\item{data}{Data frame with the outcome and predicted-probability columns.}
\item{outcome_col}{Character, name of the 0/1 outcome column.}
\item{model_cols}{Character vector of column name(s) holding predicted
probabilities for one or more models to compare.}
\item{thresholds}{Numeric vector of risk thresholds to evaluate. Default
\code{seq(0.01, 0.5, by = 0.01)}.}
}
\value{
\code{validate_model()} returns a \code{traceAsthma_validation} object (list
with \code{auroc}, \code{auprc}, \code{calibration}, and, if \code{predicted_reference} is
supplied, \code{nri}/\code{idi}); printing it shows a formatted summary.
\code{decision_curve()} returns the \code{dcurves::dca()} result object.
}
\description{
\code{validate_model()} computes the full validation metric set specified in
the TRACE-Asthma design: AUROC, AUPRC, calibration slope/intercept and
Hosmer-Lemeshow test, and (when a reference model's predictions are
supplied) net reclassification improvement (NRI) and integrated
discrimination improvement (IDI). \code{decision_curve()} wraps
\code{dcurves::dca()} to compute net benefit relative to treat-all/treat-none
strategies across a range of risk thresholds, for clinical-utility
assessment.
}
\examples{
\donttest{
if (requireNamespace("pROC", quietly = TRUE)) {
  set.seed(1)
  obs <- rbinom(100, 1, 0.4)
  pred <- plogis(rnorm(100, mean = obs - 0.5))
  val <- validate_model(obs, pred)
  print(val)
}
}
}
