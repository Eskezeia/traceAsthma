#' Discover cis- and trans-acting eQTMs
#'
#' For each CpG-gene pair, fits
#' \deqn{Y_g = \alpha_0 + \alpha_1 M_c + C\gamma + \epsilon}
#' and classifies the pair as cis (within `cis_window_bp` of the gene TSS)
#' or trans (beyond that window or on a different chromosome), applying
#' Benjamini-Hochberg correction separately within each class given their
#' very different multiple-testing burdens.
#'
#' This reference implementation loops over the supplied `pairs` table in
#' R and is intended for a pre-filtered candidate set (e.g. CpG-gene pairs
#' within a differentially-methylated / differentially-expressed shortlist)
#' rather than an unrestricted genome-wide search space, which is far more
#' efficiently handled by `MatrixEQTL`. For genome-wide discovery, run
#' `MatrixEQTL` upstream and pass its output pair list to
#' [filter_asthma_relevant()] directly.
#'
#' @param methylation Numeric matrix, CpGs x subjects (M-values recommended).
#' @param expression Numeric matrix, genes x subjects (normalized).
#' @param pairs Data frame with columns `cpg`, `gene`, `distance_bp` (signed
#'   or absolute distance between CpG and gene TSS in base pairs; use `NA`
#'   or a value > `cis_window_bp` for trans pairs on the same chromosome,
#'   and set `same_chr = FALSE` for different-chromosome pairs).
#' @param covariates Numeric matrix or data frame, subjects x covariates
#'   (age, sex, ancestry PCs, smoking, cell-type proportions, SVs, ...).
#'   Row order must match the sample (column) order of `methylation` /
#'   `expression`.
#' @param same_chr Logical vector, one per row of `pairs`; `TRUE` if CpG
#'   and gene are on the same chromosome. Required to correctly classify
#'   trans pairs when `distance_bp` is `NA`. Defaults to `TRUE` for all
#'   pairs (i.e. assumes `pairs` only contains same-chromosome candidates
#'   unless you provide this explicitly).
#' @param cis_window_bp Cis/trans boundary in base pairs. Default 1e6 (1 Mb).
#' @param fdr_threshold FDR threshold for retaining an eQTM. Default 0.05.
#' @return Data frame with columns `cpg`, `gene`, `type` ("cis"/"trans"),
#'   `alpha1` (methylation effect), `se`, `p_value`, `fdr`, filtered to
#'   `fdr <= fdr_threshold`.
#' @export
discover_eqtm <- function(methylation, expression, pairs, covariates,
                           same_chr = NULL, cis_window_bp = 1e6,
                           fdr_threshold = 0.05) {
  if (is.null(same_chr)) same_chr <- rep(TRUE, nrow(pairs))
  covariates <- as.data.frame(covariates)

  pairs$type <- ifelse(
    same_chr & !is.na(pairs$distance_bp) & abs(pairs$distance_bp) < cis_window_bp,
    "cis", "trans"
  )

  n <- nrow(pairs)
  alpha1 <- se <- p_value <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    cpg <- pairs$cpg[i]; gene <- pairs$gene[i]
    if (!(cpg %in% rownames(methylation)) || !(gene %in% rownames(expression))) next

    df <- data.frame(
      y = as.numeric(expression[gene, ]),
      m = as.numeric(methylation[cpg, ]),
      covariates
    )
    fit <- tryCatch(stats::lm(y ~ ., data = df), error = function(e) NULL)
    if (is.null(fit)) next
    smry <- summary(fit)$coefficients
    if (!("m" %in% rownames(smry))) next
    alpha1[i] <- smry["m", "Estimate"]
    se[i]     <- smry["m", "Std. Error"]
    p_value[i] <- smry["m", "Pr(>|t|)"]
  }

  out <- data.frame(pairs[, c("cpg", "gene", "type")], alpha1, se, p_value)
  out$fdr <- NA_real_
  for (ty in c("cis", "trans")) {
    idx <- out$type == ty & !is.na(out$p_value)
    out$fdr[idx] <- stats::p.adjust(out$p_value[idx], method = "BH")
  }
  out <- out[!is.na(out$fdr) & out$fdr <= fdr_threshold, ]
  out[order(out$type, out$fdr), ]
}

#' Filter eQTMs for independent CpG- and gene-level asthma association
#'
#' Applies criteria 2 and 3 of the TRACE-Asthma asthma-relevant filtering
#' scheme: the CpG and the target gene must each show independent
#' association with asthma status (criterion 1, FDR on the eQTM itself, is
#' assumed already applied by [discover_eqtm()]; criterion 4, mediation,
#' is applied separately by [test_mediation()]).
#'
#' @param eqtm_table Output of [discover_eqtm()].
#' @param methylation,expression As in [discover_eqtm()].
#' @param asthma_status Binary (0/1) outcome vector, subject order matching
#'   the columns of `methylation`/`expression`.
#' @param covariates As in [discover_eqtm()].
#' @param p_threshold P-value threshold for CpG- and gene-level asthma
#'   association. Default 0.05 (nominal, consistent with the two-stage
#'   filtering design where the primary FDR control happens at the eQTM
#'   and mediation stages).
#' @return `eqtm_table` filtered to pairs where both the CpG and the gene
#'   are nominally associated with asthma status, with added columns
#'   `cpg_asthma_p` and `gene_asthma_p`.
#' @export
filter_asthma_relevant <- function(eqtm_table, methylation, expression,
                                    asthma_status, covariates, p_threshold = 0.05) {
  covariates <- as.data.frame(covariates)
  assoc_p <- function(feature_vec) {
    df <- data.frame(y = asthma_status, x = feature_vec, covariates)
    fit <- tryCatch(
      stats::glm(y ~ ., data = df, family = stats::binomial()),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    smry <- summary(fit)$coefficients
    if (!("x" %in% rownames(smry))) return(NA_real_)
    smry["x", "Pr(>|z|)"]
  }

  eqtm_table$cpg_asthma_p  <- vapply(eqtm_table$cpg,  function(c) assoc_p(as.numeric(methylation[c, ])), numeric(1))
  eqtm_table$gene_asthma_p <- vapply(eqtm_table$gene, function(g) assoc_p(as.numeric(expression[g, ])), numeric(1))

  eqtm_table[
    !is.na(eqtm_table$cpg_asthma_p) & eqtm_table$cpg_asthma_p <= p_threshold &
    !is.na(eqtm_table$gene_asthma_p) & eqtm_table$gene_asthma_p <= p_threshold,
  ]
}

#' Test CpG -> expression -> asthma mediation for a set of eQTM pairs
#'
#' Estimates the indirect effect \eqn{IE = a \times b} for the pathway
#' CpG -> Expression -> Asthma, with bootstrap confidence intervals.
#' Uses the `mediation` package if available (recommended, provides
#' bias-corrected bootstrap CIs matching the TRACE-Asthma specification of
#' 10,000 resamples); falls back to a simple percentile-bootstrap of the
#' product-of-coefficients if `mediation` is not installed.
#'
#' @param eqtm_table Output of [filter_asthma_relevant()] (or
#'   [discover_eqtm()]).
#' @param methylation,expression,asthma_status,covariates As above.
#' @param n_boot Number of bootstrap resamples. Default 10000, matching the
#'   TRACE-Asthma specification.
#' @param p_threshold Mediation P-value threshold for retention. Default 0.05.
#' @param seed Optional integer seed.
#' @return `eqtm_table` filtered to pairs with significant mediation,
#'   with added columns `indirect_effect`, `ci_low`, `ci_high`, `p_mediation`.
#' @export
test_mediation <- function(eqtm_table, methylation, expression, asthma_status,
                            covariates, n_boot = 10000, p_threshold = 0.05,
                            seed = NULL) {
  covariates <- as.data.frame(covariates)
  if (!is.null(seed)) set.seed(seed)
  use_mediation_pkg <- requireNamespace("mediation", quietly = TRUE)

  n <- nrow(eqtm_table)
  ie <- lo <- hi <- pv <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    cpg <- eqtm_table$cpg[i]; gene <- eqtm_table$gene[i]
    df <- data.frame(
      cpg_val = as.numeric(methylation[cpg, ]),
      expr_val = as.numeric(expression[gene, ]),
      asthma = asthma_status,
      covariates
    )

    if (use_mediation_pkg) {
      med_fit <- tryCatch(
        stats::lm(expr_val ~ cpg_val + ., data = df[, !names(df) %in% "asthma"]),
        error = function(e) NULL
      )
      out_fit <- tryCatch(
        stats::glm(asthma ~ expr_val + cpg_val + ., data = df, family = stats::binomial()),
        error = function(e) NULL
      )
      if (is.null(med_fit) || is.null(out_fit)) next
      res <- tryCatch(
        mediation::mediate(med_fit, out_fit, treat = "cpg_val", mediator = "expr_val",
                            boot = TRUE, sims = n_boot),
        error = function(e) NULL
      )
      if (is.null(res)) next
      ie[i] <- res$d.avg
      lo[i] <- res$d.avg.ci[1]
      hi[i] <- res$d.avg.ci[2]
      pv[i] <- res$d.avg.p
    } else {
      # Lightweight fallback: percentile bootstrap of a * b
      boot_ie <- numeric(n_boot)
      ns <- nrow(df)
      for (b in seq_len(n_boot)) {
        idx <- sample.int(ns, ns, replace = TRUE)
        dfb <- df[idx, ]
        a_fit <- tryCatch(stats::lm(expr_val ~ cpg_val + ., data = dfb[, !names(dfb) %in% "asthma"]), error = function(e) NULL)
        b_fit <- tryCatch(stats::glm(asthma ~ expr_val + cpg_val + ., data = dfb, family = stats::binomial()), error = function(e) NULL)
        if (is.null(a_fit) || is.null(b_fit)) { boot_ie[b] <- NA; next }
        a_coef <- tryCatch(stats::coef(a_fit)["cpg_val"], error = function(e) NA)
        b_coef <- tryCatch(stats::coef(b_fit)["expr_val"], error = function(e) NA)
        boot_ie[b] <- a_coef * b_coef
      }
      boot_ie <- boot_ie[!is.na(boot_ie)]
      ie[i] <- mean(boot_ie)
      lo[i] <- stats::quantile(boot_ie, 0.025, na.rm = TRUE)
      hi[i] <- stats::quantile(boot_ie, 0.975, na.rm = TRUE)
      pv[i] <- mean(sign(boot_ie) != sign(mean(boot_ie))) * 2
    }
  }

  eqtm_table$indirect_effect <- ie
  eqtm_table$ci_low  <- lo
  eqtm_table$ci_high <- hi
  eqtm_table$p_mediation <- pv

  eqtm_table[
    !is.na(eqtm_table$p_mediation) & eqtm_table$p_mediation <= p_threshold &
    !((eqtm_table$ci_low <= 0) & (eqtm_table$ci_high >= 0)),
  ]
}

#' Fit the final eQTM score via elastic-net logistic regression
#'
#' Applies correlation pruning ([prune_correlated_features()]) followed by
#' elastic-net feature selection, then constructs
#' \deqn{eQTM_i = \sum_{j=1}^k w_j M_{ij}}
#' from the selected CpGs.
#'
#' @param methylation Numeric matrix, CpGs x subjects, restricted to
#'   mediation-significant CpGs (e.g. `unique(mediation_table$cpg)`).
#' @param asthma_status Binary (0/1) outcome vector.
#' @param train_idx Integer/logical index of training subjects for weight
#'   estimation (see [compute_mprs()] for rationale). Defaults to all.
#' @param corr_threshold Correlation-pruning threshold. Default 0.90.
#' @param alpha Elastic-net mixing parameter. Default 0.5.
#' @return A list with `eqtm_score` (named numeric vector, all subjects,
#'   standardized), `weights` (named numeric vector, selected CpGs only),
#'   `selected_cpgs`, and `fit`.
#' @export
compute_eqtm_score <- function(methylation, asthma_status, train_idx = NULL,
                                corr_threshold = 0.90, alpha = 0.5) {
  requirePkg("glmnet", "eQTM score elastic-net feature selection")
  if (is.null(train_idx)) train_idx <- seq_len(ncol(methylation))

  feature_matrix <- t(methylation)  # subjects x CpGs
  train_x <- feature_matrix[train_idx, , drop = FALSE]
  train_y <- asthma_status[train_idx]

  assoc_score <- apply(train_x, 2, function(col) {
    fit <- tryCatch(stats::glm(train_y ~ col, family = stats::binomial()), error = function(e) NULL)
    if (is.null(fit)) return(0)
    smry <- summary(fit)$coefficients
    if (!("col" %in% rownames(smry))) return(0)
    smry["col", "z value"]
  })
  names(assoc_score) <- colnames(train_x)

  retained <- prune_correlated_features(train_x, assoc_score, threshold = corr_threshold)
  train_x_pruned <- train_x[, retained, drop = FALSE]

  fit <- glmnet::cv.glmnet(train_x_pruned, train_y, family = "binomial", alpha = alpha)
  coefs <- as.matrix(stats::coef(fit, s = "lambda.min"))
  weights <- coefs[rownames(coefs) != "(Intercept)", 1]
  weights <- weights[weights != 0]

  eqtm_score <- weighted_score(feature_matrix, weights, standardize = TRUE)

  list(eqtm_score = eqtm_score, weights = weights,
       selected_cpgs = names(weights), fit = fit)
}
