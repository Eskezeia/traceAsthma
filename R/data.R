#' Synthetic example cohort for the TRACE-Asthma pipeline
#'
#' A fully synthetic, five-patient toy dataset used throughout the package
#' examples and vignette to demonstrate the end-to-end pipeline in seconds.
#' Contains two SNPs, two CpGs, four target genes, and a three-TF regulon
#' -- deliberately tiny so every intermediate score can be hand-verified.
#' **Not derived from any real cohort; for demonstration purposes only.**
#'
#' @format A list with elements:
#' \describe{
#'   \item{dosage}{Numeric matrix, 5 patients x 2 SNPs.}
#'   \item{methylation}{Numeric matrix, 2 CpGs x 5 patients.}
#'   \item{expression}{Numeric matrix, 4 genes x 5 patients.}
#'   \item{regulon_db}{Data frame, TF-target regulon table (3 TFs).}
#'   \item{clinical}{Data frame, one row per patient: `patient_id`,
#'     `asthma`, `age`, `sex`, `bmi`, `smoking`.}
#'   \item{covariates}{Data frame, eQTM-model covariates (age, sex, PC1, PC2).}
#'   \item{eqtm_pairs}{Data frame, candidate CpG-gene pairs for
#'     [discover_eqtm()].}
#'   \item{same_chr}{Logical vector flagging same-chromosome CpG-gene pairs.}
#' }
#' @source Simulated.
#' @examples
#' data(example_cohort)
#' str(example_cohort, max.level = 1)
"example_cohort"

#' Realistic, fixed-seed simulated benchmark cohort
#'
#' A larger, more realistic companion to `example_cohort`, generated once
#' with a fixed random seed via [simulate_trace_asthma_cohort()] and
#' bundled with the package so that every user gets the *exact same*
#' dataset -- useful as a shared community benchmark for testing pipeline
#' behavior, comparing results across systems/versions, and reporting
#' reproducible issues (e.g. "on `simulated_benchmark_cohort`,
#' `discover_eqtm()` gives me X"). 300 subjects, 500 CpGs, 200 genes, 10
#' transcription factors, mild two-cluster ancestry structure, and an
#' injected causal architecture with known ground truth. **Fully
#' synthetic; not derived from any real cohort.**
#'
#' Regenerate (if you change the simulator and want an updated benchmark)
#' with `data-raw/make_simulated_benchmark.R`; the seed there
#' (`20260714`) must stay fixed for the bundled object to remain a stable
#' community reference.
#'
#' @format A list with the same structure as [example_cohort], plus a
#'   `ground_truth` element:
#' \describe{
#'   \item{dosage}{Numeric matrix, 300 subjects x 30 SNPs.}
#'   \item{methylation}{Numeric matrix, 500 CpGs x 300 subjects.}
#'   \item{expression}{Numeric matrix, 200 genes x 300 subjects.}
#'   \item{regulon_db}{Data frame, TF-target regulon table (10 TFs).}
#'   \item{clinical}{Data frame: `patient_id`, `asthma`, `age`, `sex`,
#'     `bmi`, `smoking`, `ancestry_pc1`, `ancestry_pc2`.}
#'   \item{covariates}{Data frame, eQTM-model covariates.}
#'   \item{eqtm_pairs}{Data frame, candidate CpG-gene pairs (true causal
#'     pairs plus null pairs).}
#'   \item{same_chr}{Logical vector flagging same-chromosome CpG-gene pairs.}
#'   \item{ground_truth}{List: `causal_cpgs`, `causal_genes`, `active_tfs`,
#'     `true_effects`, `maf` -- the features/effect sizes actually used to
#'     generate the data, for evaluating discovery precision/recall and
#'     effect-size recovery.}
#' }
#' @source Simulated via `simulate_trace_asthma_cohort(seed = 20260714)`.
#' @seealso [simulate_trace_asthma_cohort()] to generate your own cohorts
#'   at any size/seed; [example_cohort] for a smaller instant-run toy dataset.
#' @examples
#' data(simulated_benchmark_cohort)
#' str(simulated_benchmark_cohort, max.level = 1)
#' table(simulated_benchmark_cohort$clinical$asthma)
"simulated_benchmark_cohort"
