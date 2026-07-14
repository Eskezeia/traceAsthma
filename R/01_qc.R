#' Quality control for genome-wide genotype data
#'
#' Applies standard variant- and sample-level genotype QC thresholds
#' (call rate, Hardy-Weinberg equilibrium, minor allele frequency,
#' sample missingness / heterozygosity). This is a thin, documented wrapper
#' intended to be run on PLINK-derived summary statistics; heavy lifting
#' (phasing, imputation) is expected to be done upstream with PLINK2 /
#' Eagle2 / minimac4 and is out of scope for this R package -- see the
#' package vignette for a recommended shell pipeline.
#'
#' @param variant_stats Data frame with columns `variant_id`, `call_rate`,
#'   `hwe_p`, `maf`.
#' @param sample_stats Data frame with columns `sample_id`, `missingness`,
#'   `heterozygosity_z`.
#' @param call_rate_min Minimum variant call rate. Default 0.95.
#' @param hwe_p_min Minimum HWE P-value (in controls). Default 1e-6.
#' @param maf_min Minimum minor allele frequency. Default 0.01.
#' @param sample_missingness_max Maximum per-sample missingness. Default 0.02.
#' @param het_z_max Maximum absolute heterozygosity Z-score. Default 3.
#' @return A list with `keep_variants`, `keep_samples` (character vectors of
#'   IDs passing QC) and `summary` (a one-row data frame of counts).
#' @export
qc_genotype <- function(variant_stats, sample_stats,
                         call_rate_min = 0.95, hwe_p_min = 1e-6, maf_min = 0.01,
                         sample_missingness_max = 0.02, het_z_max = 3) {
  keep_variants <- variant_stats$variant_id[
    variant_stats$call_rate >= call_rate_min &
      variant_stats$hwe_p >= hwe_p_min &
      variant_stats$maf >= maf_min
  ]
  keep_samples <- sample_stats$sample_id[
    sample_stats$missingness <= sample_missingness_max &
      abs(sample_stats$heterozygosity_z) <= het_z_max
  ]
  summary_df <- data.frame(
    variants_input = nrow(variant_stats),
    variants_retained = length(keep_variants),
    samples_input = nrow(sample_stats),
    samples_retained = length(keep_samples)
  )
  list(keep_variants = keep_variants, keep_samples = keep_samples, summary = summary_df)
}

#' Quality control for EPIC/450K DNA methylation data
#'
#' Flags low-quality probes (poor detection, SNP-overlapping,
#' cross-reactive, sex-chromosome) and low-quality samples, mirroring the
#' filtering steps described in the TRACE-Asthma Methods. Detection
#' P-values and cross-reactive/SNP-overlap probe lists are expected to be
#' precomputed upstream (e.g. via `minfi::detectionP()` and standard
#' EPIC/450K manifest annotation resources); this function operates on
#' those precomputed summaries so it has no hard Bioconductor dependency.
#'
#' @param detection_p Matrix of detection P-values, probes x samples.
#' @param flagged_probes Character vector of probe IDs to exclude a priori
#'   (SNP-overlapping, cross-reactive, sex-chromosome probes).
#' @param detection_p_max Detection P-value threshold. Default 0.01.
#' @param probe_fail_frac_max Maximum fraction of samples in which a probe
#'   may fail detection before the probe is dropped. Default 0.05.
#' @return A list with `keep_probes` and `keep_samples` character vectors,
#'   plus `summary`.
#' @export
qc_methylation <- function(detection_p, flagged_probes = character(0),
                            detection_p_max = 0.01, probe_fail_frac_max = 0.05) {
  probe_fail_frac <- rowMeans(detection_p > detection_p_max, na.rm = TRUE)
  keep_probes <- rownames(detection_p)[probe_fail_frac <= probe_fail_frac_max]
  keep_probes <- setdiff(keep_probes, flagged_probes)

  sample_fail_frac <- colMeans(detection_p > detection_p_max, na.rm = TRUE)
  keep_samples <- colnames(detection_p)[sample_fail_frac <= probe_fail_frac_max]

  summary_df <- data.frame(
    probes_input = nrow(detection_p),
    probes_retained = length(keep_probes),
    samples_input = ncol(detection_p),
    samples_retained = length(keep_samples)
  )
  list(keep_probes = keep_probes, keep_samples = keep_samples, summary = summary_df)
}

#' Quality control and low-expression filtering for RNA-seq count data
#'
#' @param counts Integer/numeric matrix of raw counts, genes x samples.
#' @param min_count Minimum count. Default 6.
#' @param min_samples Minimum number of samples in which `min_count` must
#'   be met for a gene to be retained. Default 3.
#' @return Character vector of retained gene IDs.
#' @export
qc_rnaseq <- function(counts, min_count = 6, min_samples = 3) {
  counts <- as.matrix(counts)
  keep <- rowSums(counts >= min_count) >= min_samples
  rownames(counts)[keep]
}

#' Estimate cell-type proportions by reference-based deconvolution
#'
#' Documented pass-through wrapper: reference-based deconvolution for
#' methylation (e.g. Houseman-style constrained projection, as implemented
#' in `minfi::estimateCellCounts2()`) requires a tissue-appropriate
#' reference panel and is intentionally not reimplemented here. Supply a
#' `deconvolution_fun` that returns a samples x cell-type matrix; this
#' wrapper simply validates and standardizes its output for downstream use
#' as covariates in [discover_eqtm()].
#'
#' @param methylation_matrix Numeric matrix, probes x samples.
#' @param deconvolution_fun A function taking `methylation_matrix` and
#'   returning a samples x cell-type numeric matrix (e.g. a thin wrapper
#'   around `minfi::estimateCellCounts2()` or `EpiDISH::epidish()`).
#' @return Samples x cell-type numeric matrix.
#' @export
estimate_cell_proportions <- function(methylation_matrix, deconvolution_fun) {
  if (missing(deconvolution_fun)) {
    stop(
      "estimate_cell_proportions() requires a `deconvolution_fun` implementing ",
      "reference-based deconvolution for your tissue (e.g. via minfi::estimateCellCounts2() ",
      "or EpiDISH::epidish()). This is intentionally not bundled in traceAsthma since the ",
      "correct reference panel is tissue- and platform-specific.",
      call. = FALSE
    )
  }
  props <- deconvolution_fun(methylation_matrix)
  if (!is.matrix(props)) props <- as.matrix(props)
  props
}
