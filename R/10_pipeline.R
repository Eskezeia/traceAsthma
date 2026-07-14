#' Run the full TRACE-Asthma discovery pipeline end to end
#'
#' Orchestrates QC through final model fitting on a single (training)
#' cohort, calling each pipeline stage in sequence with sensible defaults.
#' Intended for discovery-cohort model development; for unbiased
#' performance estimation wrap the relevant stages in [nested_cv()], and
#' for applying an already-fitted model to new patients use
#' [predict_asthma_risk()] instead of this function.
#'
#' This is a convenience orchestrator, not a black box: every stage it
#' calls is independently exported and documented (see `?discover_eqtm`,
#' `?compute_trace_score`, etc.), and advanced users are encouraged to
#' call the stages directly for full control over intermediate QC and
#' filtering decisions, which is what a mechanistic, clinically-deployed
#' model warrants.
#'
#' @param genotype_dosage Optional numeric matrix, subjects x variants.
#'   Omit to skip the MPRS layer entirely (see Details).
#' @param mprs_effect_sizes Optional named list of per-trait effect-size
#'   vectors, required if `genotype_dosage` is supplied.
#' @param methylation Numeric matrix, CpGs x subjects.
#' @param expression Numeric matrix, genes x subjects.
#' @param eqtm_pairs Data frame of candidate CpG-gene pairs, as expected by
#'   [discover_eqtm()].
#' @param regulon_db Regulon database, as expected by [build_regulons()].
#' @param clinical Data frame, subjects x clinical covariates (must
#'   include `asthma`; may include `age`, `sex`, `bmi`, `smoking`).
#' @param covariates Numeric matrix/data frame, subjects x eQTM/expression
#'   model covariates (age, sex, ancestry PCs, smoking, cell composition,
#'   surrogate variables).
#' @param train_idx Integer/logical index of training subjects, for use
#'   inside [nested_cv()]. Defaults to all subjects (i.e. a single,
#'   non-cross-validated discovery fit).
#' @param seed Optional integer seed for reproducibility of stochastic
#'   steps (stability selection, mediation bootstrap).
#' @return A list with all intermediate objects (`eqtm_table`,
#'   `mediation_table`, `eqtm_fit`, `tf_activity`, `trace_fit`,
#'   `mprs_fit`, `risk_model`, `nested_comparison`) plus a ready-to-deploy
#'   `deployable_model` (see [build_deployable_model()]).
#' @export
run_trace_asthma_pipeline <- function(genotype_dosage = NULL, mprs_effect_sizes = NULL,
                                       methylation, expression, eqtm_pairs, regulon_db,
                                       clinical, covariates, train_idx = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(train_idx)) train_idx <- seq_len(ncol(methylation))
  asthma_status <- clinical$asthma

  message("[1/7] Discovering cis/trans eQTMs...")
  eqtm_table <- discover_eqtm(methylation, expression, eqtm_pairs, covariates)

  message("[2/7] Filtering for independent CpG- and gene-level asthma association...")
  eqtm_table <- filter_asthma_relevant(eqtm_table, methylation, expression,
                                        asthma_status, covariates)

  message("[3/7] Testing CpG -> expression -> asthma mediation...")
  mediation_table <- test_mediation(eqtm_table, methylation, expression, asthma_status,
                                     covariates, seed = seed)

  message("[4/7] Fitting eQTM score (correlation pruning + elastic net)...")
  eqtm_meth <- methylation[unique(mediation_table$cpg), , drop = FALSE]
  eqtm_fit <- compute_eqtm_score(eqtm_meth, asthma_status, train_idx = train_idx)

  message("[5/7] Discovering candidate TFs from trans-eQTM targets and inferring activity...")
  trans_targets <- unique(mediation_table$gene[mediation_table$type == "trans"])
  regulons_candidate <- build_regulons(regulon_db, target_universe = trans_targets)
  enrichment <- test_regulon_enrichment(
    regulons_candidate, target_genes = trans_targets, background_genes = rownames(expression)
  )
  regulons_final <- regulons_candidate[regulons_candidate$tf %in% enrichment$tf, ]
  tf_activity <- infer_tf_activity(expression, regulons_final)

  message("[6/7] Fitting TRACE score (screen + prune + stability select + elastic net)...")
  trace_fit <- compute_trace_score(tf_activity, asthma_status, train_idx = train_idx, seed = seed)

  mprs_fit <- NULL
  if (!is.null(genotype_dosage)) {
    message("      Fitting MPRS...")
    trait_prs <- sapply(names(mprs_effect_sizes), function(trait) {
      compute_prs(genotype_dosage, mprs_effect_sizes[[trait]])
    })
    mprs_fit <- compute_mprs(trait_prs, asthma_status, train_idx = train_idx)
  }

  message("[7/7] Assembling analysis dataset and fitting integrated risk model...")
  model_data <- data.frame(
    asthma = asthma_status,
    eqtm = eqtm_fit$eqtm_score,
    trace = trace_fit$trace_score
  )
  if (!is.null(mprs_fit)) model_data$mprs <- mprs_fit$mprs
  for (cl in intersect(c("age", "sex", "bmi", "smoking"), colnames(clinical))) {
    model_data[[cl]] <- clinical[[cl]]
  }

  risk_model <- fit_risk_model(model_data[train_idx, , drop = FALSE])
  nested_comparison <- tryCatch(compare_nested_models(model_data[train_idx, , drop = FALSE]),
                                 error = function(e) NULL)

  deployable <- build_deployable_model(
    eqtm_weights = eqtm_fit$weights,
    tf_weights = trace_fit$weights,
    regulons = regulons_final,
    risk_model = risk_model,
    mprs_weights = if (!is.null(mprs_fit)) mprs_fit$weights else NULL,
    score_center_scale = list(
      eqtm  = c(mean = mean(eqtm_fit$eqtm_score[train_idx]), sd = stats::sd(eqtm_fit$eqtm_score[train_idx])),
      trace = c(mean = mean(trace_fit$trace_score[train_idx]), sd = stats::sd(trace_fit$trace_score[train_idx]))
    )
  )

  list(
    eqtm_table = eqtm_table,
    mediation_table = mediation_table,
    eqtm_fit = eqtm_fit,
    tf_activity = tf_activity,
    trace_fit = trace_fit,
    mprs_fit = mprs_fit,
    model_data = model_data,
    risk_model = risk_model,
    nested_comparison = nested_comparison,
    deployable_model = deployable
  )
}
