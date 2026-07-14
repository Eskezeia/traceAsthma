\name{run_package_diagnostics}
\alias{run_package_diagnostics}
\alias{print.traceAsthma_diagnostics}
\title{Run an end-to-end self-diagnostic of the installed package}
\usage{
run_package_diagnostics(cohort = NULL, verbose = TRUE)
}
\arguments{
\item{cohort}{A cohort list in the format returned by
\code{\link{simulate_trace_asthma_cohort}} (or \code{\link{example_cohort}}/
\code{\link{simulated_benchmark_cohort}}). Defaults to
\code{simulated_benchmark_cohort} -- the bundled, fixed-seed dataset -- so
that every user's diagnostic run is directly comparable.}
\item{verbose}{Logical; if \code{TRUE} (default), prints a running
stage-by-stage report as it goes, not just the final summary.}
}
\value{
Invisibly, a \code{traceAsthma_diagnostics} object (data frame with
columns \code{stage}, \code{status} ("PASS"/"FAIL"/"SKIPPED"), \code{detail});
printing it shows a formatted report.
}
\description{
The single command a new user (or a community member verifying an
installation, a new release, or a bug report) should run first. Exercises
every pipeline stage -- QC/validation, eQTM discovery with ground-truth
recovery checking, mediation, eQTM scoring, TF regulon discovery and
activity inference, TRACE scoring, the integrated risk model, and
validation metrics -- against a realistic simulated cohort with a known
causal architecture, and prints a clear PASS/FAIL/SKIPPED report per
stage. Stages requiring an optional package that isn't installed are
reported as \code{SKIPPED} (not \code{FAIL}), with the exact install
command needed -- see \code{\link{list_dependencies}}.

This function is the corrected, hardened version of a community-reported
diagnostic pattern that surfaced two real bugs (fixed in v0.3.2-0.3.3;
see \code{NEWS.md}): it always uses each candidate TF's \emph{full}
regulon for activity inference, never the \code{target_universe}-restricted
regulon used only for TF \emph{discovery} -- reusing the restricted
regulon for inference is exactly the mistake that produces a
\code{"Network is empty after intersecting..."} error from
viper/decoupleR.
}
\examples{
\donttest{
diag <- run_package_diagnostics()
print(diag)
}
}
