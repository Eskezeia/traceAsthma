\name{fit_risk_model}
\alias{fit_risk_model}
\alias{compare_nested_models}
\alias{print.traceAsthma_nested}
\alias{nested_cv}
\title{Risk model fitting, nested model comparison, and nested cross-validation}
\usage{
fit_risk_model(data, predictors = NULL)

compare_nested_models(data)

nested_cv(subject_ids, build_fun, predict_fun, outcome, k = 10, seed = NULL)
}
\arguments{
\item{data}{Data frame containing the outcome column \code{asthma} (0/1) and
any subset of predictor columns among \code{mprs}, \code{eqtm}, \code{trace}, \code{age},
\code{sex}, \code{bmi}, \code{smoking} (for \code{fit_risk_model()}), or additionally
\code{prs_asthma} (for \code{compare_nested_models()}).}
\item{predictors}{Character vector of column names in \code{data} to include as
predictors. Defaults to whichever of the standard predictor names are present.}
\item{subject_ids}{Character/integer vector of subject identifiers defining
the cross-validation partition.}
\item{build_fun}{Function taking \code{train_idx} and returning a fitted model
object.}
\item{predict_fun}{Function taking the object returned by \code{build_fun} and a
\code{test_idx}, returning predicted probabilities for the held-out subjects.}
\item{outcome}{Binary (0/1) outcome vector, same length/order as \code{subject_ids}.}
\item{k}{Number of outer folds. Default 10.}
\item{seed}{Optional integer seed for fold assignment.}
}
\value{
\code{fit_risk_model()} returns a fitted \code{glm} object.
\code{compare_nested_models()} returns a \code{traceAsthma_nested} object (list
with \code{models} and \code{comparison}); printing it shows the comparison table.
\code{nested_cv()} returns a list with \code{fold_auc}, \code{mean_auc}, \code{sd_auc}, and
\code{oof_predictions}.
}
\description{
\code{fit_risk_model()} is a thin, explicit wrapper around \code{\link[stats]{glm}}
implementing
\eqn{logit(P(Asthma)) = \beta_0 + \beta_1 MPRS + \beta_2 eQTM + \beta_3 TRACE + \beta_4 Age + \beta_5 Sex + \beta_6 BMI + \beta_7 Smoking}
(any subset of predictors may be supplied). \code{compare_nested_models()} fits
the five nested models M1 (asthma-only PRS) through M5 (full model +
clinical covariates) and compares each to its immediate predecessor via
likelihood-ratio test and AUC change. \code{nested_cv()} runs outer k-fold
cross-validation, re-executing feature construction within each training
fold via user-supplied \code{build_fun}/\code{predict_fun} to avoid information
leakage from the held-out fold.
}
\examples{
d <- data.frame(asthma = rbinom(50, 1, 0.4), age = rnorm(50, 40, 10),
                 sex = rbinom(50, 1, 0.5))
fit <- fit_risk_model(d, predictors = c("age", "sex"))
summary(fit)
}
