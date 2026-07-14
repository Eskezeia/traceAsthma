\name{run_trace_asthma_auto}
\alias{run_trace_asthma_auto}
\alias{detect_available_data}
\title{Automated pipeline that adapts to whatever data you have}
\usage{
detect_available_data(methylation, expression, clinical,
  genotype_dosage = NULL, mprs_effect_sizes = NULL,
  regulon_db = NULL, eqtm_pairs = NULL)

run_trace_asthma_auto(methylation, expression, clinical,
  genotype_dosage = NULL, mprs_effect_sizes = NULL,
  regulon_db = NULL, eqtm_pairs = NULL, train_idx = NULL, seed = NULL)
}
\arguments{
\item{methylation}{Numeric matrix, CpGs x subjects. Required.}
\item{expression}{Numeric matrix, genes x subjects. Required.}
\item{clinical}{Data frame with an \code{asthma} (0/1) column and, if
available, any of \code{age}, \code{sex}, \code{bmi}, \code{smoking},
\code{ancestry_pc1}, \code{ancestry_pc2}, .... Required (at minimum,
\code{asthma}).}
\item{genotype_dosage}{Optional numeric matrix, subjects x variants.}
\item{mprs_effect_sizes}{Optional named list of per-trait effect-size vectors.}
\item{regulon_db}{Optional TF-target regulon database; if omitted, the
TRACE layer cannot be built.}
\item{eqtm_pairs}{Optional data frame of candidate CpG-gene pairs; if
omitted, eQTM discovery (and everything downstream of it) is skipped.}
\item{train_idx}{Integer/logical index of training subjects.}
\item{seed}{Optional integer seed.}
}
\value{
\code{detect_available_data()} returns a \code{traceAsthma_plan} object
describing which layers can be built from the supplied data.
\code{run_trace_asthma_auto()} returns a list with \code{plan} plus
whichever of \code{eqtm_table}, \code{mediation_table}, \code{eqtm_fit},
\code{trace_fit}, \code{mprs_fit}, \code{risk_model}, and
\code{deployable_model} the available data supported.
}
\description{
The recommended entry point for most users. Real datasets rarely include
every modality in the full TRACE-Asthma design (genotype, methylation,
expression, and a complete clinical panel). \code{run_trace_asthma_auto()}
inspects whatever you \emph{do} have -- at minimum, methylation,
expression, and an asthma outcome for the same subjects -- and builds the
best-supported model, automatically skipping the MPRS layer if no
genotype data is supplied, skipping the TRACE layer if no regulon
database is supplied, and using only the clinical covariates
(\code{age}/\code{sex}/\code{bmi}/\code{smoking}/ancestry PCs) actually
present in \code{clinical}.
}
\examples{
data(example_cohort)
plan <- detect_available_data(
  methylation = example_cohort$methylation,
  expression  = example_cohort$expression,
  clinical    = example_cohort$clinical,
  eqtm_pairs  = example_cohort$eqtm_pairs
)
print(plan)
}
