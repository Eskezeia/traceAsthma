#' Build a TF-target regulon set from a curated database, restricted to a
#' target gene universe
#'
#' Thin wrapper around `dorothea`-style regulon tables: filters to
#' high-confidence interactions (confidence classes A-C by default,
#' matching the TRACE-Asthma specification) and, optionally, to a
#' user-supplied target gene universe (e.g. the trans-eQTM target genes
#' identified by [discover_eqtm()]), which is the recommended
#' data-driven TF discovery strategy described in the package vignette.
#'
#' @param regulon_db Data frame with columns `tf`, `target`, `confidence`
#'   (values in `"A","B","C","D","E"`), and `mor` (mode of regulation, +1
#'   for activation / -1 for repression). Typically
#'   `dorothea::dorothea_hs` or equivalent.
#' @param confidence_classes Character vector of confidence classes to
#'   retain. Default `c("A","B","C")`.
#' @param target_universe Optional character vector restricting `target` to
#'   a specific gene set (e.g. trans-eQTM targets).
#' @return Filtered regulon data frame with columns `tf`, `target`, `mor`,
#'   `confidence`.
#' @export
build_regulons <- function(regulon_db, confidence_classes = c("A", "B", "C"),
                            target_universe = NULL) {
  out <- regulon_db[regulon_db$confidence %in% confidence_classes, ]
  if (!is.null(target_universe)) {
    out <- out[out$target %in% target_universe, ]
  }
  out
}

#' Test regulon enrichment among a candidate target gene set
#'
#' Implements the data-driven TF discovery step: given a set of trans-eQTM
#' target genes, tests each candidate TF's regulon for enrichment among
#' those genes via a one-sided hypergeometric test, against a background of
#' all genes tested in the trans-eQTM analysis. TFs surviving this test
#' become candidate mediating TFs for the CpG -> TF -> gene mediation
#' analysis, rather than being chosen a priori from prior literature.
#'
#' @param regulons Output of [build_regulons()].
#' @param target_genes Character vector of genes to test for enrichment
#'   (e.g. significant trans-eQTM target genes).
#' @param background_genes Character vector of all genes in the tested
#'   universe (e.g. all genes included in the trans-eQTM search space).
#' @param fdr_threshold FDR threshold for retaining a candidate TF.
#'   Default 0.05.
#' @return Data frame with columns `tf`, `n_overlap`, `regulon_size`,
#'   `p_value`, `fdr`, sorted by `fdr`, filtered to significant TFs.
#' @export
test_regulon_enrichment <- function(regulons, target_genes, background_genes,
                                     fdr_threshold = 0.05) {
  tfs <- unique(regulons$tf)
  n_bg <- length(background_genes)
  n_target <- length(intersect(target_genes, background_genes))

  res <- do.call(rbind, lapply(tfs, function(tf) {
    regulon_genes <- unique(regulons$target[regulons$tf == tf])
    regulon_genes <- intersect(regulon_genes, background_genes)
    n_overlap <- length(intersect(regulon_genes, target_genes))
    regulon_size <- length(regulon_genes)
    # P(X >= n_overlap), X ~ Hypergeometric(n_target, n_bg - n_target, regulon_size)
    p <- stats::phyper(n_overlap - 1, n_target, n_bg - n_target, regulon_size, lower.tail = FALSE)
    data.frame(tf = tf, n_overlap = n_overlap, regulon_size = regulon_size, p_value = p)
  }))

  res$fdr <- stats::p.adjust(res$p_value, method = "BH")
  res <- res[order(res$fdr), ]
  res[!is.na(res$fdr) & res$fdr <= fdr_threshold, ]
}

#' Infer transcription factor activity from expression data and regulons
#'
#' \deqn{A_t = NES_t}
#'
#' Wraps `viper::viper()` / `decoupleR`'s VIPER-equivalent
#' (`decoupleR::run_viper()`) to estimate per-subject TF activity as the
#' normalized enrichment score (NES) of each TF's regulon among
#' differentially-expressed genes, rather than relying on the TF's own
#' transcript abundance. Prefers `viper` if installed, falls back to
#' `decoupleR`.
#'
#' @param expression Numeric matrix, genes x subjects (normalized;
#'   VST/TPM/log-CPM as used elsewhere in the pipeline).
#' @param regulons Output of [build_regulons()], with columns `tf`,
#'   `target`, `mor`.
#' @param method Either `"viper"` or `"decoupleR"`. Default `"viper"` if
#'   available, otherwise `"decoupleR"`.
#' @param min_targets Minimum regulon size to retain a TF. Default 5.
#' @return Numeric matrix, TFs x subjects, of inferred activity scores
#'   (NES).
#' @export
infer_tf_activity <- function(expression, regulons, method = c("viper", "decoupleR"),
                               min_targets = 5) {
  method <- match.arg(method)

  regulon_sizes <- table(regulons$target[regulons$target %in% rownames(expression)],
                          regulons$tf[regulons$target %in% rownames(expression)])
  # (kept for potential future use / diagnostics; not required for the calls below)

  if (method == "viper" && requireNamespace("viper", quietly = TRUE)) {
    reg_list <- split(regulons, regulons$tf)
    reg_list <- reg_list[vapply(reg_list, function(r) nrow(r) >= min_targets, logical(1))]
    viper_regulon <- lapply(reg_list, function(r) {
      list(tfmode = stats::setNames(r$mor, r$target), likelihood = rep(1, nrow(r)))
    })
    class(viper_regulon) <- "regulon"
    activity <- viper::viper(expression, viper_regulon, verbose = FALSE)
    return(activity)
  }

  requirePkg("decoupleR", "TF activity inference (VIPER unavailable)")
  net <- regulons[, c("tf", "target", "mor")]
  colnames(net) <- c("source", "target", "mor")
  res <- decoupleR::run_viper(mat = expression, network = net, minsize = min_targets)
  # decoupleR returns long format; pivot to TF x sample matrix
  wide <- stats::reshape(
    as.data.frame(res[res$statistic == "viper", c("source", "condition", "score")]),
    idvar = "source", timevar = "condition", direction = "wide"
  )
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$source
  colnames(mat) <- sub("^score\\.", "", colnames(mat))
  mat
}
