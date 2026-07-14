\name{fit_mediation_sem}
\alias{fit_mediation_sem}
\alias{summarize_mediation}
\title{Structural equation modeling of the mediation pathway}
\usage{
fit_mediation_sem(data, n_boot = 10000, seed = NULL)

summarize_mediation(sem_result)
}
\arguments{
\item{data}{Data frame with columns \code{mprs}, \code{eqtm}, \code{trace}, \code{asthma}
(standardized scores recommended).}
\item{n_boot}{Number of bootstrap resamples for CI estimation. Default 10000.}
\item{seed}{Optional integer seed.}
\item{sem_result}{Output of \code{fit_mediation_sem()}.}
}
\value{
\code{fit_mediation_sem()} returns a list with \code{fit} (the \code{lavaan} fit
object), \code{effects} (data frame with rows direct/indirect/total), and
\code{proportion_mediated}. \code{summarize_mediation()} prints a plain-language
summary and invisibly returns \code{sem_result}.
}
\description{
Fits \eqn{eQTM = a_0 + a_1 MPRS + \epsilon_1}, \eqn{TRACE = b_0 + b_1 eQTM + \epsilon_2},
\eqn{Asthma = c_0 + c_1 MPRS + c_2 eQTM + c_3 TRACE + \epsilon_3} via \code{lavaan},
decomposing the total effect of MPRS on asthma risk into a direct effect
and an indirect effect mediated through eQTM and TRACE
(\eqn{a_1 \times b_1 \times c_3}), with bootstrap confidence intervals.
}
\examples{
\donttest{
if (requireNamespace("lavaan", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(mprs = rnorm(200))
  d$eqtm <- 0.4 * d$mprs + rnorm(200)
  d$trace <- 0.5 * d$eqtm + rnorm(200)
  d$asthma <- as.integer(plogis(0.3 * d$mprs + 0.5 * d$trace) > 0.5)
  sem <- fit_mediation_sem(d, n_boot = 200)
  summarize_mediation(sem)
}
}
}
