#' Run an end-to-end self-diagnostic of the installed package
#'
#' The single command a new user (or a community member verifying an
#' installation, a new release, or a bug report) should run first. Exercises
#' every pipeline stage -- QC, validation, eQTM discovery with ground-truth
#' recovery checking, mediation, eQTM scoring, TF regulon discovery and
#' activity inference, TRACE scoring, the integrated risk model, and
#' validation metrics -- against a realistic simulated cohort with a known
#' causal architecture, and prints a clear PASS/FAIL/SKIPPED report per
#' stage.
#'
#' Stages requiring an optional package that isn't installed are reported
#' as `SKIPPED` (not `FAIL`) with the exact `install.packages()`/
#' `BiocManager::install()` call needed -- see [list_dependencies()]. This
#' function is itself the corrected, hardened version of the
#' community-reported diagnostic pattern that surfaced the
#' \code{infer_tf_activity()}/\code{test_regulon_enrichment()} issues fixed
#' in v0.3.2-0.3.3 (see \code{NEWS.md}): it always uses each candidate TF's
#' \emph{full} regulon for activity inference, never the
#' \code{target_universe}-restricted regulon used only for TF discovery,
#' which is exactly the mistake that produced the
#' \code{"Network is empty after intersecting..."} error in earlier
#' hand-written diagnostic scripts.
#'
#' @param cohort A cohort list in the format returned by
#'   [simulate_trace_asthma_cohort()] (or [example_cohort]/
#'   [simulated_benchmark_cohort]). Defaults to `simulated_benchmark_cohort`
#'   -- the bundled, fixed-seed dataset -- so that every user's diagnostic
#'   run is directly comparable.
#' @param verbose Logical; if `TRUE` (default), prints a running stage-by-
#'   stage report as it goes, not just the final summary.
#' @return Invisibly, a `traceAsthma_diagnostics` object (data frame with
#'   columns `stage`, `status` ("PASS"/"FAIL"/"SKIPPED"), `detail`);
#'   printing it shows a formatted report.
#' @export
#' @examples
#' \donttest{
#' diag <- run_package_diagnostics()
#' print(diag)
#' }
run_package_diagnostics <- function(cohort = NULL, verbose = TRUE) {
  if (is.null(cohort)) {
    utils::data("simulated_benchmark_cohort", package = "traceAsthma", envir = environment())
    cohort <- simulated_benchmark_cohort
  }

  results <- list()
  ctx <- new.env()  # carries intermediate objects forward between stages

  .run_stage <- function(name, expr_fun) {
    out <- tryCatch(
      list(status = "PASS", detail = expr_fun(), error = NULL),
      error = function(e) {
        msg <- conditionMessage(e)
        status <- if (grepl("is required for .* but is not installed", msg) ||
                       grepl("^UPSTREAM_SKIPPED:", msg) ||
                       grepl("there is no package called", msg) ||
                       grepl("could not find function", msg)) "SKIPPED" else "FAIL"
        msg <- sub("^UPSTREAM_SKIPPED: ", "", msg)
        list(status = status, detail = msg, error = e)
      }
    )
    results[[name]] <<- data.frame(stage = name, status = out$status,
                                    detail = if (is.character(out$detail)) out$detail else "ok",
                                    stringsAsFactors = FALSE)
    if (verbose) {
      icon <- switch(out$status, PASS = "[PASS]", FAIL = "[FAIL]", SKIPPED = "[SKIP]")
      cat(sprintf("%-8s %-38s %s\n", icon, name,
                  if (out$status == "PASS" && is.character(out$detail)) out$detail else
                    if (out$status != "PASS") out$detail else ""))
    }
    invisible(out$status == "PASS")
  }

  # Guard helper: if a required upstream object is missing *because* the
  # stage that would have produced it was SKIPPED (missing optional
  # package), propagate SKIPPED rather than misreporting a cascading
  # non-issue as FAIL, which would otherwise wrongly tell the user to file
  # a bug report for what is really just an uninstalled dependency.
  .require_upstream <- function(obj, upstream_stage_name, obj_description) {
    if (!is.null(obj)) return(invisible(TRUE))
    upstream_status <- if (upstream_stage_name %in% names(results)) results[[upstream_stage_name]]$status else NA
    if (identical(upstream_status, "SKIPPED")) {
      stop("UPSTREAM_SKIPPED: skipped because '", upstream_stage_name, "' was skipped (missing optional package)",
           call. = FALSE)
    }
    stop(obj_description, " unavailable from prior stage.", call. = FALSE)
  }

  if (verbose) cat(sprintf("traceAsthma self-diagnostic (package version %s)\n%s\n",
                            as.character(utils::packageVersion("traceAsthma")), strrep("-", 70)))

  # --- Stage: package/dependency overview ---
  .run_stage("dependency_overview", function() {
    deps <- if (verbose) list_dependencies() else utils::capture.output(list_dependencies())
    n_deps <- if (verbose) nrow(deps) else length(grep("^  [a-zA-Z]", deps))
    sprintf("%d optional packages listed; see above for per-stage detail", n_deps)
  })

  # --- Stage: data validation ---
  .run_stage("validate_matrix_methylation", function() {
    validate_matrix(cohort$methylation, type = "methylation"); "ok"
  })
  .run_stage("validate_matrix_expression", function() {
    validate_matrix(cohort$expression, type = "expression"); "ok"
  })
  .run_stage("validate_clinical", function() { validate_clinical(cohort$clinical); "ok" })
  .run_stage("validate_cohort_alignment", function() {
    validate_cohort_alignment(methylation = cohort$methylation, expression = cohort$expression,
                               clinical = cohort$clinical)
    "ok"
  })

  # --- Stage: eQTM discovery + ground-truth recovery ---
  .run_stage("discover_eqtm", function() {
    ctx$eqtm_table <- discover_eqtm(cohort$methylation, cohort$expression, cohort$eqtm_pairs,
                                     cohort$covariates, same_chr = cohort$same_chr, fdr_threshold = 1.0)
    sprintf("%d candidate pairs tested", nrow(ctx$eqtm_table))
  })
  .run_stage("eqtm_ground_truth_recovery", function() {
    if (is.null(ctx$eqtm_table) || nrow(ctx$eqtm_table) == 0) stop("no eQTM pairs available to check")
    causal_keys <- paste(cohort$eqtm_pairs$cpg[cohort$same_chr], cohort$eqtm_pairs$gene[cohort$same_chr])
    is_causal <- paste(ctx$eqtm_table$cpg, ctx$eqtm_table$gene) %in% causal_keys
    if (!any(is_causal) || all(is_causal)) stop("cannot compare causal vs. null pairs (need both present)")
    med_causal <- stats::median(ctx$eqtm_table$p_value[is_causal])
    med_null <- stats::median(ctx$eqtm_table$p_value[!is_causal])
    if (!(med_causal < med_null)) {
      stop(sprintf("causal pairs did NOT rank ahead of null pairs (median P causal=%.4g, null=%.4g)",
                    med_causal, med_null))
    }
    sprintf("median P causal=%.4g < null=%.4g (correctly separated)", med_causal, med_null)
  })

  # --- Stage: asthma-relevant filtering + mediation ---
  .run_stage("filter_asthma_relevant", function() {
    ctx$filtered <- filter_asthma_relevant(ctx$eqtm_table, cohort$methylation, cohort$expression,
                                            cohort$clinical$asthma, cohort$covariates, p_threshold = 1.0)
    sprintf("%d pairs retained", nrow(ctx$filtered))
  })
  .run_stage("test_mediation", function() {
    ctx$mediation_table <- test_mediation(ctx$filtered, cohort$methylation, cohort$expression,
                                           cohort$clinical$asthma, cohort$covariates,
                                           n_boot = 200, p_threshold = 1.0, seed = 1)
    sprintf("%d pairs with mediation evidence tested", nrow(ctx$mediation_table))
  })

  # --- Stage: eQTM score (requires glmnet) ---
  .run_stage("compute_eqtm_score", function() {
    sig_cpgs <- unique(ctx$mediation_table$cpg)
    if (length(sig_cpgs) < 1) stop("no CpGs available for scoring")
    ctx$eqtm_fit <- compute_eqtm_score(cohort$methylation[sig_cpgs, , drop = FALSE], cohort$clinical$asthma)
    sprintf("score built from %d CpGs", length(ctx$eqtm_fit$selected_cpgs))
  })

  # --- Stage: TF discovery (full-regulon-for-inference pattern, see NEWS 0.3.2/0.3.3) ---
  .run_stage("tf_regulon_discovery", function() {
    trans_targets <- unique(ctx$mediation_table$gene[ctx$mediation_table$type == "trans"])
    if (length(trans_targets) == 0) trans_targets <- unique(ctx$eqtm_table$gene[ctx$eqtm_table$type == "trans"])
    ctx$trans_targets <- trans_targets
    regulons_candidate <- build_regulons(cohort$regulon_db, target_universe = trans_targets)
    ctx$enrichment <- test_regulon_enrichment(regulons_candidate, target_genes = trans_targets,
                                               background_genes = rownames(cohort$expression), fdr_threshold = 1.0)
    if (nrow(ctx$enrichment) == 0) stop("no candidate TFs found (zero trans-eQTM targets available at this sample size)")
    sprintf("%d candidate TF(s): %s", nrow(ctx$enrichment), paste(head(ctx$enrichment$tf, 5), collapse = ", "))
  })
  .run_stage("infer_tf_activity", function() {
    .require_upstream(if (!is.null(ctx$enrichment) && nrow(ctx$enrichment) > 0) ctx$enrichment else NULL,
                       "tf_regulon_discovery", "Candidate TF list")
    # IMPORTANT: full regulon for the candidate TFs, NOT the target_universe-restricted
    # one used for discovery above -- this is the fix for the community-reported bug.
    regulons_full <- build_regulons(cohort$regulon_db)
    regulons_final <- regulons_full[regulons_full$tf %in% ctx$enrichment$tf, ]
    ctx$tf_activity <- infer_tf_activity(cohort$expression, regulons_final)
    sprintf("activity inferred for %d TF(s) x %d subjects", nrow(ctx$tf_activity), ncol(ctx$tf_activity))
  })
  .run_stage("compute_trace_score", function() {
    .require_upstream(ctx$tf_activity, "infer_tf_activity", "TF activity matrix")
    ctx$trace_fit <- compute_trace_score(ctx$tf_activity, cohort$clinical$asthma,
                                          fdr_threshold = 1.0, stability_threshold = 0, n_boot = 100)
    sprintf("TRACE score built from TF(s): %s", paste(names(ctx$trace_fit$weights), collapse = ", "))
  })

  # --- Stage: integrated risk model + validation ---
  .run_stage("fit_risk_model", function() {
    .require_upstream(ctx$eqtm_fit, "compute_eqtm_score", "eQTM score")
    .require_upstream(ctx$trace_fit, "compute_trace_score", "TRACE score")
    ctx$model_data <- data.frame(asthma = cohort$clinical$asthma, eqtm = ctx$eqtm_fit$eqtm_score,
                                  trace = ctx$trace_fit$trace_score, age = cohort$clinical$age,
                                  sex = cohort$clinical$sex)
    ctx$risk_model <- fit_risk_model(ctx$model_data, predictors = c("eqtm", "trace", "age", "sex"))
    "model fitted"
  })
  .run_stage("validate_model", function() {
    .require_upstream(ctx$risk_model, "fit_risk_model", "Risk model")
    pred <- stats::predict(ctx$risk_model, type = "response")
    val <- validate_model(ctx$model_data$asthma, pred)
    sprintf("AUROC = %.3f", val$auroc)
  })

  report <- structure(do.call(rbind, results), class = c("traceAsthma_diagnostics", "data.frame"))
  rownames(report) <- NULL

  if (verbose) {
    n_pass <- sum(report$status == "PASS"); n_fail <- sum(report$status == "FAIL")
    n_skip <- sum(report$status == "SKIPPED")
    cat(strrep("-", 70), "\n")
    cat(sprintf("%d PASS, %d FAIL, %d SKIPPED (missing optional packages)\n", n_pass, n_fail, n_skip))
    if (n_fail > 0) {
      cat("\nInstallation or code issue detected -- please report this at\n",
          "https://github.com/Eskezeia/traceAsthma/issues with the full output above.\n", sep = "")
    } else if (n_skip > 0) {
      cat("\nNo failures. Install the SKIPPED stages' packages (see list_dependencies()) to test them too.\n")
    } else {
      cat("\nAll stages passed. Installation verified end to end.\n")
    }
  }

  invisible(report)
}

#' @export
print.traceAsthma_diagnostics <- function(x, ...) {
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}
