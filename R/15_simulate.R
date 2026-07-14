#' Simulate genome-wide genotype dosages
#'
#' Samples independent SNPs under Hardy-Weinberg equilibrium with a
#' realistic minor allele frequency (MAF) spectrum, rather than a flat or
#' arbitrary distribution -- MAF is drawn from Uniform(0.05, 0.5), and
#' genotype dosage per subject-SNP is `rbinom(1, 2, maf)`, matching the
#' standard generative model used in PRS/QTL simulation literature (e.g.
#' HAPGEN-style marginal genotype sampling, without explicit LD -- see
#' Details).
#'
#' @details
#' This generator does **not** simulate linkage disequilibrium between
#' SNPs (each SNP's genotype is independent). This is a deliberate
#' simplification appropriate for testing pipeline mechanics (PRS
#' construction, downstream score combination) rather than for methods
#' development that specifically depends on LD structure (e.g. fine-
#' mapping). If your use case needs LD, simulate genotypes upstream with
#' a dedicated tool (e.g. `msprime`, HAPGEN2) and pass the resulting
#' dosage matrix directly to [simulate_trace_asthma_cohort()] via the
#' `genotype_dosage` override argument.
#'
#' @param n_subjects Integer, number of subjects.
#' @param n_snps Integer, number of SNPs.
#' @param maf_range Numeric length-2, range for per-SNP MAF. Default `c(0.05, 0.5)`.
#' @param seed Optional integer seed.
#' @return A list with `dosage` (numeric matrix, subjects x SNPs, values
#'   in {0,1,2}) and `maf` (named numeric vector of the realized MAF per SNP).
#' @export
simulate_genotypes <- function(n_subjects, n_snps, maf_range = c(0.05, 0.5), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  maf <- stats::runif(n_snps, maf_range[1], maf_range[2])
  dosage <- vapply(maf, function(p) stats::rbinom(n_subjects, 2, p), numeric(n_subjects))
  colnames(dosage) <- sprintf("rs%06d", seq_len(n_snps))
  rownames(dosage) <- sprintf("P%04d", seq_len(n_subjects))
  names(maf) <- colnames(dosage)
  list(dosage = dosage, maf = maf)
}

#' Simulate DNA methylation beta values with realistic marginal structure
#'
#' Real Illumina array methylation data is strongly bimodal at the
#' population level: most CpGs are either uniformly hypomethylated or
#' uniformly hypermethylated across samples (Beta distributions concentrated
#' near 0 or 1), with a minority of "intermediate" or variably-methylated
#' CpGs that carry most of the biologically interesting signal. This
#' function reproduces that mixture rather than sampling from a single
#' unimodal distribution (e.g. plain Uniform or Normal beta values), which
#' would make downstream QC/normalization code paths behave unrealistically
#' leniently compared to real array data.
#'
#' A subset of CpGs (`causal_cis_idx`, `causal_trans_idx`) are additionally
#' given a genotype-dependent shift (a synthetic meQTL effect) and/or a
#' case/control mean shift, so that [discover_eqtm()] and [test_mediation()]
#' have real signal to recover when run on the simulated data.
#'
#' @param n_subjects,n_cpgs Integers.
#' @param genotype_dosage Optional numeric matrix, subjects x SNPs (e.g.
#'   from [simulate_genotypes()]), used to inject meQTL effects.
#' @param meqtl_idx Optional integer vector of length `min(n_cpgs, ncol(genotype_dosage))`
#'   pairing CpG index `i` with SNP column `meqtl_idx[i]` for a synthetic
#'   meQTL effect on CpGs `seq_along(meqtl_idx)`. If `NULL`, no meQTL
#'   effects are injected.
#' @param meqtl_beta Numeric, per-allele effect of the paired SNP on the
#'   underlying methylation M-value for meQTL CpGs. Default 0.4.
#' @param case_control_idx Integer vector of CpG indices that additionally
#'   carry a case/control mean shift (the "differentially methylated CpG"
#'   ground truth set feeding the eQTM/mediation layers). Default: none
#'   (`integer(0)`) -- supply this together with `case_control_status` to
#'   inject a case/control effect.
#' @param case_control_status Optional 0/1 vector, length `n_subjects`;
#'   required if `case_control_idx` is non-empty.
#' @param case_control_delta Numeric, mean M-value shift in cases for
#'   `case_control_idx` CpGs. Default 0.5.
#' @param intermediate_frac Numeric in (0,1), fraction of CpGs drawn from
#'   the variable/intermediate component of the mixture (the rest split
#'   evenly between the low- and high-methylation components). Default 0.2,
#'   consistent with typical EPIC/450K array-wide marginal distributions.
#' @param seed Optional integer seed.
#' @return A list with `beta` (methylation beta values, CpGs x subjects,
#'   in \[0,1\]) and `causal_cpgs` (character vector of CpG IDs carrying an
#'   injected meQTL and/or case-control effect -- the discovery ground truth).
#' @export
simulate_methylation <- function(n_subjects, n_cpgs, genotype_dosage = NULL,
                                  meqtl_idx = NULL, meqtl_beta = 0.4,
                                  case_control_idx = integer(0),
                                  case_control_status = NULL, case_control_delta = 0.5,
                                  intermediate_frac = 0.2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cpg_ids <- sprintf("cg%08d", seq_len(n_cpgs))
  subj_ids <- if (!is.null(genotype_dosage)) rownames(genotype_dosage) else sprintf("P%04d", seq_len(n_subjects))

  n_inter <- round(intermediate_frac * n_cpgs)
  n_each_extreme <- floor((n_cpgs - n_inter) / 2)
  component <- sample(c(rep("low", n_each_extreme), rep("high", n_each_extreme),
                         rep("intermediate", n_cpgs - 2 * n_each_extreme)))

  beta <- matrix(NA_real_, nrow = n_cpgs, ncol = n_subjects, dimnames = list(cpg_ids, subj_ids))
  for (i in seq_len(n_cpgs)) {
    beta[i, ] <- switch(component[i],
      low          = stats::rbeta(n_subjects, 2, 25),
      high         = stats::rbeta(n_subjects, 25, 2),
      intermediate = stats::rbeta(n_subjects, 3, 3)
    )
  }

  # inject meQTL effects on the M-value scale, then convert back to beta
  mval <- log2(beta / (1 - beta))
  if (!is.null(meqtl_idx) && !is.null(genotype_dosage)) {
    if (length(meqtl_idx) > n_cpgs) {
      warning(sprintf(
        "length(meqtl_idx) = %d exceeds n_cpgs = %d; truncating meqtl_idx to the first %d entries.",
        length(meqtl_idx), n_cpgs, n_cpgs
      ))
      meqtl_idx <- meqtl_idx[seq_len(n_cpgs)]
    }
    for (i in seq_along(meqtl_idx)) {
      snp_col <- meqtl_idx[i]
      if (is.na(snp_col) || snp_col > ncol(genotype_dosage)) next
      mval[i, ] <- mval[i, ] + meqtl_beta * genotype_dosage[, snp_col]
    }
  }
  if (length(case_control_idx) > 0) {
    if (is.null(case_control_status)) {
      stop("case_control_status is required when case_control_idx is non-empty.", call. = FALSE)
    }
    if (any(case_control_idx > n_cpgs | case_control_idx < 1)) {
      stop(sprintf("case_control_idx contains value(s) outside [1, n_cpgs = %d].", n_cpgs), call. = FALSE)
    }
    for (i in case_control_idx) {
      mval[i, ] <- mval[i, ] + case_control_delta * case_control_status
    }
  }
  beta <- 2^mval / (1 + 2^mval)  # back-transform, keeps values in (0,1)

  list(beta = beta, causal_cpgs = cpg_ids[union(meqtl_idx %||% integer(0), case_control_idx)])
}

#' Simulate a TF-target regulon
#'
#' Assigns each transcription factor a random target-gene set of realistic
#' size (10-40 targets, right-skewed so a few "hub" TFs have larger
#' regulons, matching curated databases like DoRothEA), with mode-of-
#' regulation (activation/repression) roughly 80/20 activating, and
#' confidence classes distributed across A-C.
#'
#' @param gene_ids Character vector of gene IDs to draw targets from.
#' @param n_tfs Integer, number of TFs.
#' @param min_targets,max_targets Integers, regulon size range. Default
#'   `c(10, 40)`.
#' @param seed Optional integer seed.
#' @return Data frame with columns `tf`, `target`, `mor` (+1/-1),
#'   `confidence` (A/B/C), matching the format expected by [build_regulons()].
#' @export
simulate_regulon <- function(gene_ids, n_tfs, min_targets = 10, max_targets = 40, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  tf_ids <- sprintf("TF%02d", seq_len(n_tfs))
  hub_weight <- stats::rgamma(n_tfs, shape = 2)  # right-skewed -> a few hub TFs
  sizes <- pmin(max_targets, pmax(min_targets,
                round(min_targets + hub_weight / max(hub_weight) * (max_targets - min_targets))))

  do.call(rbind, lapply(seq_len(n_tfs), function(i) {
    targets <- sample(gene_ids, min(sizes[i], length(gene_ids)))
    data.frame(
      tf = tf_ids[i], target = targets,
      mor = sample(c(1, -1), length(targets), replace = TRUE, prob = c(0.8, 0.2)),
      confidence = sample(c("A", "B", "C"), length(targets), replace = TRUE, prob = c(0.3, 0.4, 0.3)),
      stringsAsFactors = FALSE
    )
  }))
}

#' Simulate RNA-seq gene expression counts
#'
#' Generates negative-binomial count data with sample-level library size
#' factors and gene-level dispersion, matching the generative model
#' underlying standard RNA-seq differential expression tools (DESeq2,
#' edgeR) rather than a Normal/log-Normal approximation -- so downstream
#' QC (`qc_rnaseq()`) and normalization behave the way they would on real
#' count data. A subset of genes are made downstream targets of causal
#' CpGs (the eQTM ground truth: `mean expression ~ f(paired CpG beta)`)
#' and/or driven by TF regulon membership (genes sharing an "active" TF in
#' cases get a coordinated expression shift), so [infer_tf_activity()] and
#' the TRACE layer have real, recoverable structure.
#'
#' @param n_subjects,n_genes Integers.
#' @param methylation_beta Optional numeric matrix, CpGs x subjects (e.g.
#'   from [simulate_methylation()]), used to inject eQTM effects.
#' @param eqtm_pairs Optional data frame with columns `cpg_idx`
#'   (row index into `methylation_beta`) and `gene_idx` (index into the
#'   simulated gene set), defining which CpG drives which gene's mean
#'   expression, and `effect` (numeric, effect of CpG beta on log-mean
#'   expression).
#' @param regulon Optional regulon data frame (from [simulate_regulon()]);
#'   if supplied along with `active_tfs` and `case_control_status`, target
#'   genes of the active TF(s) receive a coordinated expression shift in
#'   cases, giving [infer_tf_activity()] real signal to detect.
#' @param active_tfs Character vector of TF IDs (from `regulon$tf`) that
#'   are "biologically active" in cases -- the TRACE-layer ground truth.
#' @param active_tf_effect Numeric, log-fold-change applied to active TFs'
#'   activating targets in cases (and the opposite sign for repressive
#'   targets). Default 0.6.
#' @param case_control_status Optional 0/1 vector, required if `regulon`/
#'   `active_tfs` are supplied.
#' @param dispersion Numeric, negative-binomial dispersion (1/size).
#'   Default 0.15 (moderate, typical of bulk RNA-seq).
#' @param baseline_mean_range Numeric length-2, range for gene-level
#'   baseline mean count (library-size-normalized). Default `c(20, 500)`.
#' @param seed Optional integer seed.
#' @return A list with `counts` (integer matrix, genes x subjects) and
#'   `causal_genes` (character vector of gene IDs carrying an injected
#'   eQTM and/or TF-activity effect -- the discovery ground truth).
#' @export
simulate_expression <- function(n_subjects, n_genes, methylation_beta = NULL,
                                 eqtm_pairs = NULL, regulon = NULL, active_tfs = NULL,
                                 active_tf_effect = 0.6, case_control_status = NULL,
                                 dispersion = 0.15, baseline_mean_range = c(20, 500),
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  gene_ids <- sprintf("GENE%05d", seq_len(n_genes))
  subj_ids <- if (!is.null(methylation_beta)) colnames(methylation_beta) else sprintf("P%04d", seq_len(n_subjects))

  baseline_mean <- stats::runif(n_genes, baseline_mean_range[1], baseline_mean_range[2])
  size_factor <- stats::rlnorm(length(subj_ids), meanlog = 0, sdlog = 0.2)  # library size variation

  log_mu <- matrix(log(baseline_mean), nrow = n_genes, ncol = length(subj_ids))

  causal_genes <- character(0)
  if (!is.null(eqtm_pairs) && !is.null(methylation_beta)) {
    for (k in seq_len(nrow(eqtm_pairs))) {
      gi <- eqtm_pairs$gene_idx[k]; ci <- eqtm_pairs$cpg_idx[k]; eff <- eqtm_pairs$effect[k]
      if (gi > n_genes || ci > nrow(methylation_beta)) next
      log_mu[gi, ] <- log_mu[gi, ] + eff * (methylation_beta[ci, ] - mean(methylation_beta[ci, ]))
      causal_genes <- union(causal_genes, gene_ids[gi])
    }
  }
  if (!is.null(regulon) && !is.null(active_tfs) && !is.null(case_control_status)) {
    for (tf in active_tfs) {
      targets <- regulon[regulon$tf == tf, ]
      idx <- match(targets$target, gene_ids)
      valid <- !is.na(idx)
      idx <- idx[valid]; mor <- targets$mor[valid]
      if (length(idx) == 0) next
      shift <- outer(mor * active_tf_effect, case_control_status)
      log_mu[idx, ] <- log_mu[idx, ] + shift
      causal_genes <- union(causal_genes, gene_ids[idx])
    }
  }

  log_mu <- sweep(log_mu, 2, log(size_factor), "+")
  mu <- exp(log_mu)
  size <- 1 / dispersion
  counts <- matrix(stats::rnbinom(n_genes * length(subj_ids), mu = mu, size = size),
                    nrow = n_genes, dimnames = list(gene_ids, subj_ids))

  list(counts = counts, causal_genes = causal_genes)
}

#' Simulate a clinical covariate panel
#'
#' Generates age, sex, BMI, smoking status, and genetic-ancestry principal
#' components with realistic marginal distributions and cross-correlations
#' (BMI right-skewed and weakly age-associated; smoking status with
#' realistic never/former/current proportions), rather than independent
#' uniform/normal draws for every field.
#'
#' @param n_subjects Integer.
#' @param ancestry_groups Integer, number of discrete ancestry clusters to
#'   simulate PCs for (each cluster is a Gaussian blob in PC space).
#'   Default 1 (single-ancestry cohort; set > 1 for a multi-ancestry design).
#' @param seed Optional integer seed.
#' @return Data frame, one row per subject: `age`, `sex` (0/1), `bmi`,
#'   `smoking` (factor: never/former/current), `ancestry_pc1`, `ancestry_pc2`.
#' @export
simulate_clinical <- function(n_subjects, ancestry_groups = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  age <- round(pmin(90, pmax(5, stats::rnorm(n_subjects, 38, 15))))
  sex <- stats::rbinom(n_subjects, 1, 0.52)
  bmi <- round(pmin(55, pmax(14, stats::rgamma(n_subjects, shape = 8, scale = 3) + 12)), 1)
  smoking <- factor(sample(c("never", "former", "current"), n_subjects, replace = TRUE,
                            prob = c(0.55, 0.30, 0.15)), levels = c("never", "former", "current"))

  group <- sample(seq_len(ancestry_groups), n_subjects, replace = TRUE)
  centers <- if (ancestry_groups > 1) {
    matrix(stats::rnorm(ancestry_groups * 2, sd = 3), ncol = 2)
  } else matrix(0, nrow = 1, ncol = 2)
  pcs <- centers[group, , drop = FALSE] + matrix(stats::rnorm(n_subjects * 2, sd = 0.5), ncol = 2)

  data.frame(
    age = age, sex = sex, bmi = bmi, smoking = smoking,
    ancestry_pc1 = round(pcs[, 1], 3), ancestry_pc2 = round(pcs[, 2], 3),
    row.names = sprintf("P%04d", seq_len(n_subjects))
  )
}

#' Simulate a full TRACE-Asthma-format cohort with a known ground truth
#'
#' The top-level simulator: generates genotypes, methylation, expression,
#' a TF regulon, and a clinical panel that are internally consistent with
#' each other and with a specified causal architecture
#' (SNP -> methylation -> eQTM -> TF activity -> TRACE -> asthma), then
#' simulates asthma case/control status from a logistic model combining
#' those true effects -- so that running the real `traceAsthma` pipeline
#' on the output should recover approximately the specified effect sizes
#' and flag approximately the correct causal CpGs/genes/TFs. Intended for
#' end-to-end pipeline testing (statistical power, discovery calibration,
#' code correctness) before real data is available; **not** intended to
#' produce numbers to report as biological findings.
#'
#' @param n_subjects Integer, cohort size. Default 300.
#' @param n_snps,n_cpgs,n_genes,n_tfs Integers, feature-space sizes.
#'   Defaults 30, 500, 200, 8 -- large enough for realistic multiple-testing
#'   behavior while remaining fast to simulate and to run through discovery.
#' @param n_causal_cpgs,n_causal_genes,n_active_tfs Integers, how many
#'   features carry true injected signal. Defaults 15, 10, 2.
#' @param true_effects Named list of true log-odds effect sizes for the
#'   final asthma-risk model, used to generate the outcome:
#'   `list(mprs = 0.4, eqtm = 0.5, trace = 0.7, age = 0.01, sex = 0.3)`
#'   by default. Override any subset; unspecified terms default as shown.
#' @param prevalence Numeric in (0,1), approximate target asthma
#'   prevalence; the intercept is calibrated (via a short root-finding
#'   step) to hit this prevalence given the other true effects. Default 0.35.
#' @param ancestry_groups Integer, passed to [simulate_clinical()]. Default 1.
#' @param seed Optional integer seed for full reproducibility.
#' @return A list in the same format as `example_cohort`, with
#'   \describe{
#'     \item{dosage, methylation, expression, regulon_db, clinical, covariates, eqtm_pairs, same_chr}{
#'       As in `example_cohort`, ready to pass directly to
#'       [run_trace_asthma_pipeline()] or [run_trace_asthma_auto()].}
#'     \item{ground_truth}{A list with `causal_cpgs`, `causal_genes`,
#'       `active_tfs`, and `true_effects` -- the features/effect sizes
#'       actually used to generate the data, for evaluating discovery
#'       precision/recall and effect-size recovery once you run the
#'       pipeline on this simulated cohort.}
#'   }
#' @export
#' @examples
#' sim <- simulate_trace_asthma_cohort(n_subjects = 100, n_cpgs = 100,
#'                                      n_genes = 60, n_tfs = 5, seed = 1)
#' str(sim$ground_truth)
#' table(sim$clinical$asthma)
simulate_trace_asthma_cohort <- function(n_subjects = 300, n_snps = 30, n_cpgs = 500,
                                          n_genes = 200, n_tfs = 8,
                                          n_causal_cpgs = 15, n_causal_genes = 10,
                                          n_active_tfs = 2, true_effects = list(),
                                          prevalence = 0.35, ancestry_groups = 1,
                                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  defaults <- list(mprs = 0.4, eqtm = 0.5, trace = 0.7, age = 0.01, sex = 0.3)
  true_effects <- utils::modifyList(defaults, true_effects)

  clinical <- simulate_clinical(n_subjects, ancestry_groups = ancestry_groups)
  geno <- simulate_genotypes(n_subjects, n_snps)
  rownames(geno$dosage) <- rownames(clinical)

  # provisional case/control labels for injecting methylation/expression
  # effects, refined below once the true logistic model is fit
  provisional_status <- stats::rbinom(n_subjects, 1, prevalence)

  meqtl_idx <- if (n_snps > 0) sample.int(n_snps, min(n_causal_cpgs, n_snps, n_cpgs), replace = TRUE) else NULL
  cc_idx <- sample.int(n_cpgs, min(n_causal_cpgs, n_cpgs))
  meth <- simulate_methylation(
    n_subjects, n_cpgs, genotype_dosage = geno$dosage, meqtl_idx = meqtl_idx,
    case_control_idx = cc_idx, case_control_status = provisional_status
  )
  colnames(meth$beta) <- rownames(clinical)

  gene_ids <- sprintf("GENE%05d", seq_len(n_genes))
  eqtm_gene_idx <- sample.int(n_genes, min(n_causal_genes, n_genes))
  eqtm_pair_spec <- data.frame(
    cpg_idx = sample(cc_idx, length(eqtm_gene_idx), replace = TRUE),
    gene_idx = eqtm_gene_idx,
    effect = stats::runif(length(eqtm_gene_idx), 1.5, 3.5)
  )

  regulon <- simulate_regulon(gene_ids, n_tfs)
  active_tfs <- sample(unique(regulon$tf), min(n_active_tfs, length(unique(regulon$tf))))

  expr <- simulate_expression(
    n_subjects, n_genes, methylation_beta = meth$beta, eqtm_pairs = eqtm_pair_spec,
    regulon = regulon, active_tfs = active_tfs, case_control_status = provisional_status
  )
  colnames(expr$counts) <- rownames(clinical)

  # --- derive approximate true MPRS / eQTM / TRACE scores for outcome generation ---
  mprs_true <- z_score(rowSums(geno$dosage[, seq_len(min(5, n_snps)), drop = FALSE]))
  eqtm_true <- z_score(colMeans(meth$beta[cc_idx, , drop = FALSE]))
  active_targets <- unique(regulon$target[regulon$tf %in% active_tfs])
  active_target_idx <- match(active_targets, gene_ids)
  active_target_idx <- active_target_idx[!is.na(active_target_idx)]
  trace_true <- if (length(active_target_idx) > 0) {
    z_score(colMeans(log1p(expr$counts[active_target_idx, , drop = FALSE])))
  } else rep(0, n_subjects)

  lin <- true_effects$mprs * mprs_true + true_effects$eqtm * eqtm_true +
    true_effects$trace * trace_true + true_effects$age * (clinical$age - mean(clinical$age)) +
    true_effects$sex * clinical$sex

  # calibrate intercept to hit target prevalence
  target_logit_offset <- stats::uniroot(
    function(b0) mean(stats::plogis(b0 + lin)) - prevalence,
    interval = c(-10, 10)
  )$root
  p_true <- stats::plogis(target_logit_offset + lin)
  asthma <- stats::rbinom(n_subjects, 1, p_true)
  clinical$asthma <- asthma
  clinical$patient_id <- rownames(clinical)
  clinical <- clinical[, c("patient_id", setdiff(colnames(clinical), "patient_id"))]

  # candidate eQTM pairs for discover_eqtm(): the true causal pairs plus a
  # random sample of null pairs, mirroring a realistic pre-screened
  # candidate list rather than the full CpG x gene cross product
  n_null_pairs <- min(50, n_cpgs * n_genes - nrow(eqtm_pair_spec))
  null_cpg <- sample.int(n_cpgs, n_null_pairs, replace = TRUE)
  null_gene <- sample.int(n_genes, n_null_pairs, replace = TRUE)
  eqtm_pairs <- data.frame(
    cpg = c(rownames(meth$beta)[eqtm_pair_spec$cpg_idx], rownames(meth$beta)[null_cpg]),
    gene = c(gene_ids[eqtm_pair_spec$gene_idx], gene_ids[null_gene]),
    distance_bp = c(rep(50000, nrow(eqtm_pair_spec)), rep(NA, n_null_pairs)),
    stringsAsFactors = FALSE
  )
  same_chr <- c(rep(TRUE, nrow(eqtm_pair_spec)), rep(FALSE, n_null_pairs))

  covariates <- data.frame(
    age = clinical$age, sex = clinical$sex,
    ancestry_pc1 = clinical$ancestry_pc1, ancestry_pc2 = clinical$ancestry_pc2,
    row.names = rownames(clinical)
  )

  list(
    dosage = geno$dosage,
    methylation = meth$beta,
    expression = expr$counts,
    regulon_db = regulon,
    clinical = clinical,
    covariates = covariates,
    eqtm_pairs = eqtm_pairs,
    same_chr = same_chr,
    ground_truth = list(
      causal_cpgs = meth$causal_cpgs,
      causal_genes = expr$causal_genes,
      active_tfs = active_tfs,
      true_effects = true_effects,
      maf = geno$maf
    )
  )
}
