#' Detect which TRACE-Asthma layers a given set of inputs supports
#'
#' Real users rarely have every data type described in the TRACE-Asthma
#' design (genotype, methylation, expression, and a full clinical panel
#' including BMI and smoking). This function inspects whatever you *do*
#' have and reports which layers of the pipeline can be built, so
#' [run_trace_asthma_auto()] can construct the best model your data
#' supports rather than requiring the full specification.
#'
#' The only hard requirement is methylation, expression, and an asthma
#' outcome for the same subjects -- everything else (genotype/MPRS,
#' BMI, smoking, ancestry) is optional and used opportunistically.
#'
#' @param methylation Numeric matrix, CpGs x subjects. Required.
#' @param expression Numeric matrix, genes x subjects. Required.
#' @param clinical Data frame with an `asthma` (0/1) column and, if
#'   available, any of `age`, `sex`, `bmi`, `smoking`,
#'   `ancestry_pc1`, `ancestry_pc2`, .... Required (at minimum, `asthma`).
#' @param genotype_dosage Optional numeric matrix, subjects x variants.
#' @param mprs_effect_sizes Optional named list of per-trait effect-size
#'   vectors (see [compute_mprs()]).
#' @param regulon_db Optional TF-target regulon database (see
#'   [build_regulons()]); if omitted, the TRACE layer cannot be built.
#' @param eqtm_pairs Optional data frame of candidate CpG-gene pairs (see
#'   [discover_eqtm()]); if omitted, all CpG-gene combinations present in
#'   both `methylation` and any `eqtm_pairs`-free discovery mode are
#'   considered too computationally expensive to enumerate automatically,
#'   so eQTM discovery is skipped and a warning is issued (supply
#'   `eqtm_pairs` for a targeted candidate list, e.g. from prior
#'   differential methylation / expression analysis).
#' @return An object of class `traceAsthma_plan`: a list describing which
#'   layers are available (`mprs`, `eqtm`, `trace`) and which clinical
#'   covariates were found (`clinical_covariates`), plus the assembled
#'   `covariates` data frame to use in eQTM/TF association models
#'   (built automatically from whichever of `age`, `sex`, and
#'   `ancestry_pc*` columns are present in `clinical`, if `covariates`
#'   is not supplied directly).
#' @export
detect_available_data <- function(methylation, expression, clinical,
                                   genotype_dosage = NULL, mprs_effect_sizes = NULL,
                                   regulon_db = NULL, eqtm_pairs = NULL) {
  if (missing(methylation) || missing(expression) || missing(clinical)) {
    stop("methylation, expression, and clinical (with an `asthma` column) are the minimum required inputs.",
         call. = FALSE)
  }
  if (!("asthma" %in% colnames(clinical))) {
    stop("`clinical` must contain an `asthma` (0/1) outcome column.", call. = FALSE)
  }

  n_subjects <- ncol(methylation)
  if (ncol(expression) != n_subjects || nrow(clinical) != n_subjects) {
    stop("methylation, expression, and clinical must all describe the same subjects ",
         "in the same order (ncol(methylation) == ncol(expression) == nrow(clinical)).",
         call. = FALSE)
  }

  mprs_available <- !is.null(genotype_dosage) && !is.null(mprs_effect_sizes)
  trace_available <- !is.null(regulon_db)
  eqtm_available <- !is.null(eqtm_pairs)

  ancestry_cols <- grep("^ancestry_pc", colnames(clinical), value = TRUE)
  clinical_risk_covariates <- intersect(c("age", "sex", "bmi", "smoking"), colnames(clinical))
  model_covariates <- c(intersect(c("age", "sex"), colnames(clinical)), ancestry_cols)

  covariates <- if (length(model_covariates) > 0) {
    df <- clinical[, model_covariates, drop = FALSE]
    if ("sex" %in% colnames(df) && !is.numeric(df$sex)) df$sex <- as.numeric(as.factor(df$sex)) - 1
    rownames(df) <- colnames(methylation)
    df
  } else {
    NULL
  }

  plan <- list(
    n_subjects = n_subjects,
    mprs_available = mprs_available,
    eqtm_available = eqtm_available,
    trace_available = trace_available,
    clinical_risk_covariates = clinical_risk_covariates,
    model_covariates = model_covariates,
    covariates = covariates
  )
  structure(plan, class = "traceAsthma_plan")
}

#' @export
print.traceAsthma_plan <- function(x, ...) {
  yn <- function(b) if (b) "yes" else "no"
  cat("TRACE-Asthma data availability plan\n")
  cat(sprintf("  Subjects                         : %d\n", x$n_subjects))
  cat(sprintf("  eQTM candidate pairs supplied     : %s\n", yn(x$eqtm_available)))
  cat(sprintf("  MPRS layer (genotype + effects)   : %s\n", yn(x$mprs_available)))
  cat(sprintf("  TRACE layer (regulon database)    : %s\n", yn(x$trace_available)))
  cat(sprintf("  Clinical covariates for final model: %s\n",
              if (length(x$clinical_risk_covariates)) paste(x$clinical_risk_covariates, collapse = ", ") else "none"))
  cat(sprintf("  eQTM/TF association-model covariates: %s\n",
              if (length(x$model_covariates)) paste(x$model_covariates, collapse = ", ") else "none (intercept-only)"))
  if (!x$eqtm_available) {
    cat("\n  NOTE: no `eqtm_pairs` supplied -> eQTM discovery will be skipped.\n",
        "        Supply a candidate CpG-gene pair list (e.g. from a prior\n",
        "        differential methylation/expression shortlist) to enable it.\n", sep = "")
  }
  if (!x$trace_available) {
    cat("\n  NOTE: no `regulon_db` supplied -> TRACE layer will be skipped.\n",
        "        A public regulon set (e.g. DoRothEA) is sufficient; see\n",
        "        vignette(\"trace-asthma-pipeline\").\n", sep = "")
  }
  if (!x$mprs_available) {
    cat("\n  NOTE: no genotype/MPRS supplied -> genetic layer will be skipped.\n",
        "        The model will be built from eQTM + TRACE + clinical data only.\n", sep = "")
  }
  invisible(x)
}

#' Automated end-to-end pipeline that adapts to whatever data you have
#'
#' The recommended entry point for most users. Wraps
#' [detect_available_data()] and [run_trace_asthma_pipeline()]: builds the
#' maximal model your supplied data supports, skipping layers you don't
#' have data for (MPRS, TRACE) rather than requiring the full
#' specification, and prints a clear report of what was included. At
#' minimum this needs `methylation`, `expression`, and `clinical$asthma`
#' for the same set of subjects -- e.g. a user with only CpG methylation,
#' gene expression, asthma status, age, sex, and ancestry can run this
#' directly and get an eQTM(+ TRACE, if a regulon database is supplied)
#' model with age/sex/ancestry-adjusted association testing, with no
#' genotype or BMI/smoking data required.
#'
#' @inheritParams detect_available_data
#' @param train_idx,seed As in [run_trace_asthma_pipeline()].
#' @return A list with `plan` (the `traceAsthma_plan` used) and everything
#'   [run_trace_asthma_pipeline()] returns for the layers that were built.
#'   If neither `eqtm_pairs` nor `regulon_db` were supplied, only QC-level
#'   outputs and the availability report are returned (there is nothing
#'   to model).
#' @export
run_trace_asthma_auto <- function(methylation, expression, clinical,
                                   genotype_dosage = NULL, mprs_effect_sizes = NULL,
                                   regulon_db = NULL, eqtm_pairs = NULL,
                                   train_idx = NULL, seed = NULL) {
  plan <- detect_available_data(methylation, expression, clinical, genotype_dosage,
                                 mprs_effect_sizes, regulon_db, eqtm_pairs)
  print(plan)

  if (!plan$eqtm_available) {
    message("\nNo eQTM candidate pairs supplied; returning the availability plan only. ",
            "Supply `eqtm_pairs` and re-run to build a model.")
    return(list(plan = plan))
  }

  covariates <- plan$covariates
  if (is.null(covariates)) {
    covariates <- data.frame(intercept_only = rep(1, plan$n_subjects))
    rownames(covariates) <- colnames(methylation)
  }

  if (plan$trace_available) {
    message("\nRunning full pipeline (eQTM + TRACE", if (plan$mprs_available) " + MPRS" else "", ")...\n")
    result <- run_trace_asthma_pipeline(
      genotype_dosage = if (plan$mprs_available) genotype_dosage else NULL,
      mprs_effect_sizes = if (plan$mprs_available) mprs_effect_sizes else NULL,
      methylation = methylation, expression = expression,
      eqtm_pairs = eqtm_pairs, regulon_db = regulon_db,
      clinical = clinical, covariates = covariates,
      train_idx = train_idx, seed = seed
    )
    result$plan <- plan
    return(result)
  }

  # eQTM-only path (no regulon database -> no TRACE layer)
  message("\nRunning eQTM-only pipeline (no regulon database supplied for the TRACE layer)...\n")
  asthma_status <- clinical$asthma
  eqtm_table <- discover_eqtm(methylation, expression, eqtm_pairs, covariates)
  eqtm_table <- filter_asthma_relevant(eqtm_table, methylation, expression, asthma_status, covariates)
  mediation_table <- test_mediation(eqtm_table, methylation, expression, asthma_status, covariates, seed = seed)
  eqtm_meth <- methylation[unique(mediation_table$cpg), , drop = FALSE]
  eqtm_fit <- compute_eqtm_score(eqtm_meth, asthma_status, train_idx = train_idx)

  model_data <- data.frame(asthma = asthma_status, eqtm = eqtm_fit$eqtm_score)
  for (cl in plan$clinical_risk_covariates) model_data[[cl]] <- clinical[[cl]]
  if (!is.null(train_idx)) {
    risk_model <- fit_risk_model(model_data[train_idx, , drop = FALSE])
  } else {
    risk_model <- fit_risk_model(model_data)
  }

  list(
    plan = plan, eqtm_table = eqtm_table, mediation_table = mediation_table,
    eqtm_fit = eqtm_fit, model_data = model_data, risk_model = risk_model
  )
}

#' List traceAsthma's dependencies and what each is used for
#'
#' Prints (and invisibly returns) a table of every non-base package
#' `traceAsthma` uses, whether it is required (`Imports`) or
#' optional-per-stage (`Suggests`), and which pipeline stage(s) need it --
#' so a user can install only what their planned analysis requires rather
#' than the full stack.
#'
#' @return Invisibly, a data frame with columns `package`, `type`, `stage`.
#' @export
list_dependencies <- function() {
  deps <- data.frame(
    package = c("glmnet", "pROC", "PRROC", "mediation", "lavaan", "viper",
                "decoupleR", "minfi", "sva", "dorothea", "xgboost", "randomForest",
                "rms", "dcurves", "shiny", "writexl", "testthat", "knitr", "rmarkdown"),
    type = c("Suggests", "Suggests", "Suggests", "Suggests", "Suggests", "Suggests",
             "Suggests", "Suggests", "Suggests", "Suggests", "Suggests", "Suggests",
             "Suggests", "Suggests", "Suggests", "Suggests", "Suggests", "Suggests", "Suggests"),
    stage = c(
      "eQTM/TRACE/MPRS elastic-net feature selection (compute_eqtm_score, compute_trace_score, compute_mprs, stability_select)",
      "Discrimination metrics: AUC, ROC curves, nested-model comparison, decision curve x-axis (compare_nested_models, validate_model, plot_roc_curves)",
      "AUPRC computation (validate_model)",
      "CpG -> expression -> asthma mediation testing (test_mediation); optional, has a base-R bootstrap fallback",
      "Structural equation modeling of the MPRS -> eQTM -> TRACE -> Asthma pathway (fit_mediation_sem, plot_mediation_path)",
      "Transcription factor activity inference, preferred backend (infer_tf_activity)",
      "Transcription factor activity inference, fallback backend if viper is unavailable (infer_tf_activity)",
      "Reference-based cell-type deconvolution helper input (estimate_cell_proportions); not a hard dependency, user-supplied function",
      "Surrogate variable analysis for RNA-seq batch correction (recommended upstream step, not called directly by this package)",
      "Source of curated TF-target regulons (build_regulons input; download from Bioconductor)",
      "Benchmark nonlinear risk model alternative to logistic regression (optional, user-driven)",
      "Benchmark nonlinear risk model alternative to logistic regression (optional, user-driven)",
      "Alternative calibration/validation utilities (optional; base-R calibration is implemented in validate_model)",
      "Decision curve analysis (decision_curve)",
      "Point-of-care prediction app (inst/shiny/app.R)",
      "Native .xlsx table export (export_table); falls back to CSV if absent",
      "Running the package's unit test suite (development only)",
      "Building vignettes (development/CRAN submission only)",
      "Building vignettes (development/CRAN submission only)"
    ),
    stringsAsFactors = FALSE
  )
  cat("traceAsthma dependencies (all optional / Suggests, installed on demand):\n\n")
  for (i in seq_len(nrow(deps))) {
    cat(sprintf("  %-14s %s\n", deps$package[i], deps$stage[i]))
  }
  cat("\nBase-R-only functionality (no dependencies beyond stats/utils/tools) includes:\n",
      "  QC (qc_genotype, qc_methylation, qc_rnaseq), discover_eqtm, filter_asthma_relevant,\n",
      "  build_regulons, test_regulon_enrichment, fit_risk_model, all format_*_table()\n",
      "  functions, and all plot_*() figure functions (base graphics only).\n", sep = "")
  invisible(deps)
}
