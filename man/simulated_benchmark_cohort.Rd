\name{simulated_benchmark_cohort}
\alias{simulated_benchmark_cohort}
\docType{data}
\title{Realistic, fixed-seed simulated benchmark cohort}
\format{
A list with the same structure as \code{\link{example_cohort}}, plus a
\code{ground_truth} element:
\describe{
  \item{dosage}{Numeric matrix, 300 subjects x 30 SNPs.}
  \item{methylation}{Numeric matrix, 500 CpGs x 300 subjects.}
  \item{expression}{Numeric matrix, 200 genes x 300 subjects.}
  \item{regulon_db}{Data frame, TF-target regulon table (10 TFs).}
  \item{clinical}{Data frame: \code{patient_id}, \code{asthma}, \code{age},
    \code{sex}, \code{bmi}, \code{smoking}, \code{ancestry_pc1}, \code{ancestry_pc2}.}
  \item{covariates}{Data frame, eQTM-model covariates.}
  \item{eqtm_pairs}{Data frame, candidate CpG-gene pairs (true causal pairs
    plus null pairs).}
  \item{same_chr}{Logical vector flagging same-chromosome CpG-gene pairs.}
  \item{ground_truth}{List: \code{causal_cpgs}, \code{causal_genes},
    \code{active_tfs}, \code{true_effects}, \code{maf} -- the features/effect
    sizes actually used to generate the data, for evaluating discovery
    precision/recall and effect-size recovery.}
}
}
\source{
Simulated via \code{simulate_trace_asthma_cohort(seed = 20260714)}.
}
\description{
A larger, more realistic companion to \code{example_cohort}, generated once
with a fixed random seed via \code{\link{simulate_trace_asthma_cohort}} and
bundled with the package so that every user gets the \strong{exact same}
dataset -- useful as a shared community benchmark for testing pipeline
behavior, comparing results across systems/versions, and reporting
reproducible issues (e.g. "on \code{simulated_benchmark_cohort},
\code{discover_eqtm()} gives me X"). 300 subjects, 500 CpGs, 200 genes, 10
transcription factors, mild two-cluster ancestry structure, and an injected
causal architecture with known ground truth. \strong{Fully synthetic; not
derived from any real cohort.}

Regenerate (if you change the simulator and want an updated benchmark) with
\code{data-raw/make_simulated_benchmark.R}; the seed there (\code{20260714})
must stay fixed for the bundled object to remain a stable community reference.
}
\seealso{
\code{\link{simulate_trace_asthma_cohort}} to generate your own cohorts at
any size/seed; \code{\link{example_cohort}} for a smaller instant-run toy dataset.
}
\examples{
data(simulated_benchmark_cohort)
str(simulated_benchmark_cohort, max.level = 1)
table(simulated_benchmark_cohort$clinical$asthma)
}
\keyword{datasets}
