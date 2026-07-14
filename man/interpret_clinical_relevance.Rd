\name{interpret_clinical_relevance}
\alias{interpret_clinical_relevance}
\alias{print.traceAsthma_interpretation}
\title{Translate a fitted risk model into plain-language clinical relevance}
\usage{
interpret_clinical_relevance(risk_model, data, score_var = NULL,
  percentiles = c(0.10, 0.50, 0.90), validation = NULL, sem_result = NULL,
  outcome_label = "asthma")
}
\arguments{
\item{risk_model}{A fitted \code{glm} object (e.g. from \code{fit_risk_model()}).}
\item{data}{The data frame used to fit \code{risk_model} (or an equivalent
reference sample), used to determine percentile values of \code{score_var}
and median/mode values at which to hold other covariates fixed.}
\item{score_var}{Character, the model term to profile across percentiles
(e.g. \code{"trace"}, \code{"eqtm"}, \code{"mprs"}). If \code{NULL}
(default), the first of \code{"trace"}, \code{"eqtm"}, \code{"mprs"}
found among \code{risk_model}'s terms is used.}
\item{percentiles}{Numeric vector of percentiles (0-1) at which to profile
\code{score_var}. Default \code{c(0.10, 0.50, 0.90)} (low/median/high).}
\item{validation}{Optional output of \code{validate_model()}; if supplied,
its calibration slope/intercept and Hosmer-Lemeshow P are translated into
a plain-language trustworthiness statement.}
\item{sem_result}{Optional output of \code{fit_mediation_sem()}; if
supplied, the proportion-mediated estimate is translated into a
plain-language mechanistic statement naming the modifiable pathway.}
\item{outcome_label}{Character, used in the narrative in place of the
generic word "outcome" (default \code{"asthma"}).}
}
\value{
A \code{traceAsthma_interpretation} object: a list with \code{absolute_risk}
(data frame of predicted probability at each profile, with absolute risk
difference from the lowest profile), \code{calibration_statement},
\code{mediation_statement} (each character or \code{NULL}), and
\code{narrative} (the full assembled paragraph). Printing the object shows
the narrative followed by the table.
}
\description{
Statistical output (odds ratios, AUC, mediated proportion) answers "is
this associated" and "does this discriminate" -- this function answers
the question a clinician actually asks: how much does this change the
estimated risk for a patient, can that number be trusted, and is there a
modifiable step in the pathway. It computes absolute predicted risk at
low/median/high values of a chosen score (holding other covariates at
their sample median/mode), attaches a calibration statement if a
\code{validate_model()} result is supplied, and attaches a mechanistic
mediation statement if a \code{fit_mediation_sem()} result is supplied.
}
\examples{
set.seed(1)
n <- 60
d <- data.frame(asthma = rbinom(n, 1, 0.4), age = round(rnorm(n, 40, 10)),
                 sex = rbinom(n, 1, 0.5), trace = rnorm(n))
fit <- fit_risk_model(d, predictors = c("age", "sex", "trace"))
interp <- interpret_clinical_relevance(fit, d)
print(interp)
}
