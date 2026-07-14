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
#'   (e.g. significant trans-eQTM target genes). If empty, or if no TF's
#'   regulon overlaps `background_genes`, this function returns a
#'   zero-row (but correctly structured) result rather than erroring --
#'   check `nrow()` of the result and see the Details section below for
#'   what an empty result means and how to respond to it.
#' @param background_genes Character vector of all genes in the tested
#'   universe (e.g. all genes included in the trans-eQTM search space).
#' @param fdr_threshold FDR threshold for retaining a candidate TF.
#'   Default 0.05.
#' @return Data frame with columns `tf`, `n_overlap`, `regulon_size`,
#'   `p_value`, `fdr`, sorted by `fdr`, filtered to significant TFs. If
#'   `target_genes` is empty or no candidate TF has any regulon overlap
#'   with `background_genes`, returns a zero-row data frame with these
#'   same columns (not an error) -- this commonly happens when an
#'   upstream discovery step (e.g. [discover_eqtm()] restricted to
#'   trans-only pairs) found zero significant trans associations, which
#'   is a legitimate result at modest sample sizes, not necessarily a bug.
#' @export
test_regulon_enrichment <- function(regulons, target_genes, background_genes,
                                     fdr_threshold = 0.05) {
  empty_result <- data.frame(tf = character(0), n_overlap = integer(0),
                              regulon_size = integer(0), p_value = numeric(0),
                              fdr = numeric(0))

  if (length(target_genes) == 0) {
    warning("test_regulon_enrichment(): `target_genes` is empty; returning a zero-row result. ",
            "This typically means the upstream discovery step found no significant trans (or ",
            "other target-defining) associations to test for TF enrichment -- a legitimate ",
            "outcome, especially at modest sample sizes, not necessarily an error.")
    return(empty_result)
  }

  tfs <- unique(regulons$tf)
  if (length(tfs) == 0) {
    warning("test_regulon_enrichment(): `regulons` has no rows (e.g. after restricting to a ",
            "target_universe with no overlap); returning a zero-row result.")
    return(empty_result)
  }

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

  if (is.null(res) || nrow(res) == 0) return(empty_result)

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
#' @param min_targets Minimum regulon size to retain a TF. Default 5. Applied
#'   to each TF's target set *after* intersecting with `rownames(expression)`
#'   -- if you restricted `regulons` to a small target universe upstream
#'   (e.g. via `build_regulons(..., target_universe = ...)` for a data-driven
#'   TF discovery step), pass the TF's *full* regulon here instead, or this
#'   function will correctly, but perhaps surprisingly, drop TFs that had
#'   plenty of true targets but too few survived the earlier restriction.
#' @return Numeric matrix, TFs x subjects, of inferred activity scores
#'   (NES).
#' @export
infer_tf_activity <- function(expression, regulons, method = c("viper", "decoupleR"),
                               min_targets = 5) {
  method <- match.arg(method)

  # Validate regulon sizes up front, in terms of genes actually present in
  # `expression`, and fail with a specific, actionable message naming which
  # TFs are short and by how much -- rather than letting viper/decoupleR
  # surface their own internal "network is empty" error, which gives no
  # indication of *why* or *which TF* caused it.
  regulons <- regulons[regulons$target %in% rownames(expression), , drop = FALSE]
  tf_sizes <- table(regulons$tf)
  tf_sizes_all <- stats::setNames(as.integer(tf_sizes), names(tf_sizes))
  eligible_tfs <- names(tf_sizes_all)[tf_sizes_all >= min_targets]

  if (length(eligible_tfs) == 0) {
    short_tfs <- if (length(tf_sizes_all) > 0) {
      paste(sprintf("%s (%d target%s in expression)", names(tf_sizes_all), tf_sizes_all,
                     ifelse(tf_sizes_all == 1, "", "s")), collapse = "; ")
    } else {
      "none of the supplied regulon's TFs have any targets present in rownames(expression)"
    }
    stop(
      "infer_tf_activity(): no TF has at least min_targets = ", min_targets,
      " targets overlapping rownames(expression). Per-TF target counts after intersecting ",
      "with expression: ", short_tfs, ".\n",
      "This commonly happens when `regulons` was restricted to a small target_universe ",
      "upstream (e.g. via build_regulons(..., target_universe = <trans-eQTM genes>)) for a ",
      "TF-discovery/enrichment step -- that restriction is appropriate for discovery, but for ",
      "activity inference you should pass each candidate TF's FULL regulon instead: ",
      "e.g. `regulons_full <- build_regulons(regulon_db); ",
      "regulons_full <- regulons_full[regulons_full$tf %in% candidate_tfs, ]`. ",
      "Alternatively, lower `min_targets`.",
      call. = FALSE
    )
  }
  if (length(eligible_tfs) < length(tf_sizes_all)) {
    dropped <- setdiff(names(tf_sizes_all), eligible_tfs)
    warning(sprintf(
      "infer_tf_activity(): dropping %d TF(s) with fewer than min_targets = %d targets in expression: %s.",
      length(dropped), min_targets, paste(dropped, collapse = ", ")
    ))
  }
  regulons <- regulons[regulons$tf %in% eligible_tfs, , drop = FALSE]

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
