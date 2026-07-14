\name{run_trace_asthma_pipeline}
\alias{run_trace_asthma_pipeline}
\title{Run the full TRACE-Asthma discovery pipeline end to end}
\usage{
run_trace_asthma_pipeline(genotype_dosage = NULL, mprs_effect_sizes = NULL,
  methylation, expression, eqtm_pairs, regulon_db, clinical, covariates,
  train_idx = NULL, seed = NULL)
}
\arguments{
\item{genotype_dosage}{Optional numeric matrix, subjects x variants. Omit
to skip the MPRS layer entirely.}
\item{mprs_effect_sizes}{Optional named list of per-trait effect-size
vectors, required if \code{genotype_dosage} is supplied.}
\item{methylation}{Numeric matrix, CpGs x subjects.}
\item{expression}{Numeric matrix, genes x subjects.}
\item{eqtm_pairs}{Data frame of candidate CpG-gene pairs, as expected by
\code{discover_eqtm()}.}
\item{regulon_db}{Regulon database, as expected by \code{build_regulons()}.}
\item{clinical}{Data frame, subjects x clinical covariates (must include
\code{asthma}; may include \code{age}, \code{sex}, \code{bmi}, \code{smoking}).}
\item{covariates}{Numeric matrix/data frame, subjects x eQTM-model
covariates (age, sex, ancestry PCs, smoking, cell composition, surrogate
variables).}
\item{train_idx}{Integer/logical index of training subjects, for use
inside \code{nested_cv()}. Defaults to all subjects.}
\item{seed}{Optional integer seed for reproducibility of stochastic steps.}
}
\value{
A list with all intermediate objects (\code{eqtm_table},
\code{mediation_table}, \code{eqtm_fit}, \code{tf_activity},
\code{trace_fit}, \code{mprs_fit}, \code{risk_model},
\code{nested_comparison}) plus a ready-to-deploy \code{deployable_model}
(see \code{build_deployable_model()}).
}
\description{
Orchestrates QC through final model fitting on a single (training) cohort,
calling each pipeline stage in sequence with sensible defaults. Every
stage it calls is independently exported and documented; advanced users
are encouraged to call the stages directly for full control over
intermediate QC and filtering decisions.
}
\seealso{
\code{\link{predict_asthma_risk}} to apply a fitted model to new patients;
\code{nested_cv()} for unbiased performance estimation.
}
