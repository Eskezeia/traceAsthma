\name{qc_genotype}
\alias{qc_genotype}
\alias{qc_methylation}
\alias{qc_rnaseq}
\alias{estimate_cell_proportions}
\title{Quality control for genotype, methylation, and RNA-seq data}
\usage{
qc_genotype(variant_stats, sample_stats, call_rate_min = 0.95,
  hwe_p_min = 1e-6, maf_min = 0.01, sample_missingness_max = 0.02,
  het_z_max = 3)

qc_methylation(detection_p, flagged_probes = character(0),
  detection_p_max = 0.01, probe_fail_frac_max = 0.05)

qc_rnaseq(counts, min_count = 6, min_samples = 3)

estimate_cell_proportions(methylation_matrix, deconvolution_fun)
}
\arguments{
\item{variant_stats}{Data frame with columns \code{variant_id}, \code{call_rate},
\code{hwe_p}, \code{maf}.}
\item{sample_stats}{Data frame with columns \code{sample_id}, \code{missingness},
\code{heterozygosity_z}.}
\item{call_rate_min}{Minimum variant call rate. Default 0.95.}
\item{hwe_p_min}{Minimum Hardy-Weinberg equilibrium P-value. Default 1e-6.}
\item{maf_min}{Minimum minor allele frequency. Default 0.01.}
\item{sample_missingness_max}{Maximum per-sample missingness. Default 0.02.}
\item{het_z_max}{Maximum absolute heterozygosity Z-score. Default 3.}
\item{detection_p}{Matrix of detection P-values, probes x samples.}
\item{flagged_probes}{Character vector of probe IDs to exclude a priori
(SNP-overlapping, cross-reactive, sex-chromosome probes).}
\item{detection_p_max}{Detection P-value threshold. Default 0.01.}
\item{probe_fail_frac_max}{Maximum fraction of samples in which a probe/sample
may fail detection. Default 0.05.}
\item{counts}{Integer/numeric matrix of raw RNA-seq counts, genes x samples.}
\item{min_count}{Minimum count. Default 6.}
\item{min_samples}{Minimum number of samples meeting \code{min_count}. Default 3.}
\item{methylation_matrix}{Numeric matrix, probes x samples.}
\item{deconvolution_fun}{A function taking \code{methylation_matrix} and
returning a samples x cell-type numeric matrix.}
}
\value{
\code{qc_genotype()} and \code{qc_methylation()} return a list with
\code{keep_variants}/\code{keep_probes} and \code{keep_samples} (character
vectors of IDs passing QC) and \code{summary}. \code{qc_rnaseq()} returns a
character vector of retained gene IDs. \code{estimate_cell_proportions()}
returns a samples x cell-type numeric matrix.
}
\description{
Applies standard variant-, probe-, and sample-level quality control
thresholds for genome-wide genotype, DNA methylation, and RNA-seq count
data, operating on precomputed summary statistics (call rates, detection
P-values) rather than raw platform files, so this package has no hard
dependency on the upstream genomics toolchain (PLINK2, minfi, STAR).
\code{estimate_cell_proportions()} is a documented pass-through: reference-based
deconvolution requires a tissue-appropriate reference panel and is
intentionally not reimplemented here -- supply a \code{deconvolution_fun}
wrapping e.g. \code{minfi::estimateCellCounts2()} or \code{EpiDISH::epidish()}.
}
\examples{
vstats <- data.frame(variant_id = c("rs1","rs2","rs3"), call_rate = c(0.99,0.80,0.97),
                      hwe_p = c(0.5,0.5,1e-8), maf = c(0.2,0.3,0.02))
sstats <- data.frame(sample_id = c("P1","P2"), missingness = c(0.01,0.03),
                      heterozygosity_z = c(0.5,-0.2))
qc_genotype(vstats, sstats)

counts <- matrix(rpois(40, 10), nrow = 4, dimnames = list(paste0("g",1:4), paste0("P",1:10)))
qc_rnaseq(counts, min_count = 5, min_samples = 3)
}
