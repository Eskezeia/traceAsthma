\name{simulate_trace_asthma_cohort}
\alias{simulate_trace_asthma_cohort}
\alias{simulate_genotypes}
\alias{simulate_methylation}
\alias{simulate_expression}
\alias{simulate_regulon}
\alias{simulate_clinical}
\title{Simulate a realistic TRACE-Asthma cohort with a known ground truth}
\usage{
simulate_trace_asthma_cohort(n_subjects = 300, n_snps = 30, n_cpgs = 500,
  n_genes = 200, n_tfs = 8, n_causal_cpgs = 15, n_causal_genes = 10,
  n_active_tfs = 2, true_effects = list(), prevalence = 0.35,
  ancestry_groups = 1, seed = NULL)

simulate_genotypes(n_subjects, n_snps, maf_range = c(0.05, 0.5), seed = NULL)

simulate_methylation(n_subjects, n_cpgs, genotype_dosage = NULL,
  meqtl_idx = NULL, meqtl_beta = 0.4, case_control_idx = integer(0),
  case_control_status = NULL, case_control_delta = 0.5,
  intermediate_frac = 0.2, seed = NULL)

simulate_expression(n_subjects, n_genes, methylation_beta = NULL,
  eqtm_pairs = NULL, regulon = NULL, active_tfs = NULL,
  active_tf_effect = 0.6, case_control_status = NULL, dispersion = 0.15,
  baseline_mean_range = c(20, 500), seed = NULL)

simulate_regulon(gene_ids, n_tfs, min_targets = 10, max_targets = 40, seed = NULL)

simulate_clinical(n_subjects, ancestry_groups = 1, seed = NULL)
}
\arguments{
\item{n_subjects}{Integer, cohort size.}
\item{n_snps, n_cpgs, n_genes, n_tfs}{Integers, feature-space sizes.}
\item{n_causal_cpgs, n_causal_genes, n_active_tfs}{Integers, how many
features carry true injected signal.}
\item{true_effects}{Named list of true log-odds effect sizes used to
generate the outcome, e.g. \code{list(mprs = 0.4, eqtm = 0.5, trace = 0.7,
age = 0.01, sex = 0.3)}. Override any subset.}
\item{prevalence}{Numeric in (0,1), approximate target asthma prevalence;
the intercept is calibrated to hit this via root-finding.}
\item{ancestry_groups}{Integer, number of discrete ancestry clusters.}
\item{seed}{Optional integer seed for full reproducibility.}
\item{maf_range}{Numeric length-2, per-SNP minor allele frequency range.}
\item{genotype_dosage}{Optional dosage matrix for injecting meQTL effects.}
\item{meqtl_idx, meqtl_beta}{CpG-to-SNP pairing and per-allele effect size
for synthetic meQTL signal.}
\item{case_control_idx, case_control_status, case_control_delta}{Which
CpGs/genes carry a case-control mean shift, the case/control label vector,
and the shift magnitude.}
\item{intermediate_frac}{Fraction of CpGs drawn from the variable/
intermediate component of the bimodal methylation mixture.}
\item{methylation_beta, eqtm_pairs}{Inputs for injecting eQTM
(CpG -> expression) effects into simulated counts.}
\item{regulon, active_tfs, active_tf_effect}{Inputs for injecting
TF-activity-driven, case-coordinated expression shifts.}
\item{dispersion, baseline_mean_range}{Negative-binomial RNA-seq count
generation parameters.}
\item{gene_ids}{Character vector of gene IDs to draw regulon targets from.}
\item{min_targets, max_targets}{Regulon size range per TF.}
}
\value{
\code{simulate_trace_asthma_cohort()} returns a list in the same format as
\code{example_cohort} (\code{dosage}, \code{methylation}, \code{expression},
\code{regulon_db}, \code{clinical}, \code{covariates}, \code{eqtm_pairs},
\code{same_chr}), ready to pass directly to \code{run_trace_asthma_pipeline()}
or \code{run_trace_asthma_auto()}, plus a \code{ground_truth} element
listing the causal CpGs/genes/TFs and true effect sizes used to generate
the data -- for evaluating discovery precision/recall and effect-size
recovery. The individual \code{simulate_*()} functions are exported for
building custom simulation scenarios and each return the specific
data type plus its own ground-truth vector where applicable.
}
\description{
Generates a synthetic cohort whose individual data types match the
distributional properties of real genetic, epigenetic, transcriptomic,
and clinical data (Hardy-Weinberg genotype sampling with a realistic MAF
spectrum; bimodal Beta-mixture methylation values as seen on real
Illumina arrays; negative-binomial RNA-seq counts with library-size
factors, matching the DESeq2/edgeR generative model; realistic clinical
covariate distributions and ancestry structure) and are internally
consistent with a specified causal architecture
(SNP -> methylation -> eQTM -> TF activity -> TRACE -> asthma), with
asthma case/control status simulated from a logistic model using
user-specified true effect sizes. Intended for end-to-end pipeline
testing (statistical power, discovery calibration, code correctness)
before real data is available -- \strong{not} intended to produce numbers
to report as biological findings.
}
\examples{
sim <- simulate_trace_asthma_cohort(n_subjects = 100, n_cpgs = 100,
                                     n_genes = 60, n_tfs = 5, seed = 1)
str(sim$ground_truth)
table(sim$clinical$asthma)
}
