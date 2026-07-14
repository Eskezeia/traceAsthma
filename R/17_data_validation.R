#' Validate a methylation, expression, or genotype dosage matrix
#'
#' Checks structural and value-range assumptions before any analysis
#' function runs, so a malformed real-world dataset fails immediately with
#' a specific, actionable message rather than producing a silent
#' downstream error or, worse, a silently wrong result deep inside
#' [discover_eqtm()] or [infer_tf_activity()].
#'
#' @param x The matrix to validate.
#' @param type One of `"methylation"` (values expected in \[0,1\], or
#'   already M-value transformed -- see `allow_mvalues`), `"expression"`
#'   (non-negative, typically counts or normalized values), or
#'   `"dosage"` (values expected in {0,1,2}, allowing non-integer
#'   imputed dosages when `allow_imputed = TRUE`).
#' @param allow_mvalues Logical, only used when `type = "methylation"`;
#'   if `TRUE`, values outside \[0,1\] are permitted (treated as
#'   logit-transformed M-values) with a message rather than an error.
#'   Default `FALSE`.
#' @param allow_imputed Logical, only used when `type = "dosage"`; if
#'   `TRUE`, non-integer dosages in \[0,2\] (typical of imputed genotypes)
#'   are permitted. Default `TRUE`.
#' @param max_missing_frac Numeric in \[0,1\], maximum allowed fraction of
#'   missing values (overall). Default 0.2. Rows/columns exceeding this
#'   individually are reported by name in the error, not just counted.
#' @param min_samples Minimum number of columns (subjects) required.
#'   Default 3 (the minimum for any variance-based downstream step to be
#'   meaningful; realistic analyses need far more, but this catches
#'   obviously malformed input like a single-column upload).
#' @return Invisibly, `TRUE` if all checks pass. Throws an informative
#'   error naming the specific problem (and offending row/column IDs,
#'   truncated to the first 10) otherwise.
#' @export
#' @examples
#' m <- matrix(runif(20, 0, 1), nrow = 4, dimnames = list(paste0("cg", 1:4), paste0("P", 1:5)))
#' validate_matrix(m, type = "methylation")
validate_matrix <- function(x, type = c("methylation", "expression", "dosage"),
                             allow_mvalues = FALSE, allow_imputed = TRUE,
                             max_missing_frac = 0.2, min_samples = 3) {
  type <- match.arg(type)
  nm <- deparse(substitute(x))

  if (!is.matrix(x)) {
    stop(sprintf("`%s` must be a matrix (features x subjects); got class '%s'. ",
                  nm, paste(class(x), collapse = "/")),
         "If you have a data.frame, convert with as.matrix() after moving the ID column to rownames ",
         "(see import_data(..., as_matrix = TRUE)).", call. = FALSE)
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(sprintf("`%s` must have both row names (feature IDs) and column names (subject IDs).", nm),
         call. = FALSE)
  }
  if (ncol(x) < min_samples) {
    stop(sprintf("`%s` has only %d subject(s); at least %d are required for any meaningful analysis.",
                  nm, ncol(x), min_samples), call. = FALSE)
  }
  if (anyDuplicated(rownames(x))) {
    dup <- unique(rownames(x)[duplicated(rownames(x))])
    stop(sprintf("`%s` has duplicated feature IDs (row names): %s%s",
                  nm, paste(utils::head(dup, 10), collapse = ", "),
                  if (length(dup) > 10) ", ..." else ""), call. = FALSE)
  }
  if (anyDuplicated(colnames(x))) {
    dup <- unique(colnames(x)[duplicated(colnames(x))])
    stop(sprintf("`%s` has duplicated subject IDs (column names): %s%s",
                  nm, paste(utils::head(dup, 10), collapse = ", "),
                  if (length(dup) > 10) ", ..." else ""), call. = FALSE)
  }

  miss_frac <- mean(is.na(x))
  if (miss_frac > max_missing_frac) {
    row_miss <- rowMeans(is.na(x))
    worst <- names(sort(row_miss, decreasing = TRUE))[seq_len(min(10, sum(row_miss > 0)))]
    stop(sprintf(
      "`%s` is %.1f%% missing, exceeding max_missing_frac = %.1f%%. Features with the most missing values: %s. ",
      nm, 100 * miss_frac, 100 * max_missing_frac, paste(worst, collapse = ", ")),
      "Impute or filter missing data before proceeding (e.g. exclude features/subjects above a ",
      "missingness threshold, or impute -- imputation strategy is a scientific decision left to the user, ",
      "not silently applied by traceAsthma).", call. = FALSE)
  }

  vals <- x[!is.na(x)]
  if (type == "methylation" && !allow_mvalues) {
    if (any(vals < 0 | vals > 1)) {
      n_bad <- sum(vals < 0 | vals > 1)
      stop(sprintf(
        "`%s` contains %d value(s) outside [0,1], but type = \"methylation\" expects beta values in [0,1]. ",
        nm, n_bad),
        "If these are logit-transformed M-values, call with allow_mvalues = TRUE.", call. = FALSE)
    }
  }
  if (type == "expression") {
    if (any(vals < 0)) {
      stop(sprintf("`%s` contains negative values, but type = \"expression\" expects non-negative counts or normalized expression values.", nm),
           call. = FALSE)
    }
  }
  if (type == "dosage") {
    if (allow_imputed) {
      if (any(vals < 0 | vals > 2)) {
        stop(sprintf("`%s` contains values outside [0,2], but type = \"dosage\" expects genotype dosages in [0,2].", nm),
             call. = FALSE)
      }
    } else if (any(!(vals %in% c(0, 1, 2)))) {
      stop(sprintf(
        "`%s` contains non-integer or out-of-range values, but type = \"dosage\" with allow_imputed = FALSE ",
        "expects hard-called genotypes in {0,1,2}. Set allow_imputed = TRUE for imputed/dosage data.", nm),
        call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Validate a clinical/covariate data frame
#'
#' @param clinical Data frame to validate.
#' @param required_cols Character vector of columns that must be present.
#'   Default `"asthma"` (the minimum required by every pipeline entry point).
#' @param id_col Character, name of the subject-ID column, if present, used
#'   to check for duplicates. Default `"patient_id"`; set `NULL` to skip
#'   this check if your data frame has no explicit ID column (row names
#'   are then assumed to be subject IDs).
#' @return Invisibly, `TRUE` if all checks pass; throws an informative
#'   error otherwise.
#' @export
#' @examples
#' cl <- data.frame(patient_id = paste0("P", 1:5), asthma = c(0,1,0,1,1),
#'                   age = c(30,40,25,55,33), sex = c(0,1,0,1,0))
#' validate_clinical(cl)
validate_clinical <- function(clinical, required_cols = "asthma", id_col = "patient_id") {
  nm <- deparse(substitute(clinical))
  if (!is.data.frame(clinical)) {
    stop(sprintf("`%s` must be a data frame; got class '%s'.", nm, paste(class(clinical), collapse = "/")),
         call. = FALSE)
  }
  missing_cols <- setdiff(required_cols, colnames(clinical))
  if (length(missing_cols) > 0) {
    stop(sprintf("`%s` is missing required column(s): %s.", nm, paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }
  if ("asthma" %in% colnames(clinical)) {
    bad <- !(clinical$asthma %in% c(0, 1, NA))
    if (any(bad)) {
      stop(sprintf("`%s$asthma` must be coded 0/1 (control/case); found other value(s): %s.",
                    nm, paste(unique(clinical$asthma[bad]), collapse = ", ")), call. = FALSE)
    }
    if (any(is.na(clinical$asthma))) {
      warning(sprintf("`%s$asthma` has %d missing value(s); affected subjects will be dropped by most modeling functions.",
                       nm, sum(is.na(clinical$asthma))))
    }
  }
  if (!is.null(id_col) && id_col %in% colnames(clinical)) {
    if (anyDuplicated(clinical[[id_col]])) {
      dup <- unique(clinical[[id_col]][duplicated(clinical[[id_col]])])
      stop(sprintf("`%s$%s` has duplicated subject IDs: %s.", nm, id_col, paste(utils::head(dup, 10), collapse = ", ")),
           call. = FALSE)
    }
  }
  for (numeric_field in intersect(c("age", "bmi"), colnames(clinical))) {
    v <- clinical[[numeric_field]]
    if (!is.numeric(v)) {
      stop(sprintf("`%s$%s` must be numeric; got class '%s'.", nm, numeric_field, paste(class(v), collapse = "/")),
           call. = FALSE)
    }
    bounds <- if (numeric_field == "age") c(0, 120) else c(10, 80)
    if (any(v < bounds[1] | v > bounds[2], na.rm = TRUE)) {
      warning(sprintf("`%s$%s` has value(s) outside the plausible range [%g, %g]; check for data entry or unit errors.",
                       nm, numeric_field, bounds[1], bounds[2]))
    }
  }
  if ("sex" %in% colnames(clinical) && !all(clinical$sex %in% c(0, 1, NA))) {
    warning(sprintf(
      "`%s$sex` is not coded 0/1. traceAsthma's covariate-building functions expect numeric 0/1 ",
      "(see the data dictionary, vignette(\"data-requirements\")); non-numeric sex codings ",
      "(e.g. \"F\"/\"M\") will be auto-converted where possible but explicit recoding is safer.", nm))
  }
  invisible(TRUE)
}

#' Validate that a set of omics matrices and a clinical table describe the
#' same subjects in the same order
#'
#' The single most common real-world data bug: methylation, expression,
#' and clinical tables assembled by different upstream pipelines/people,
#' with subjects in different orders or a partial ID mismatch. Every
#' function in this package assumes column-order alignment (see the data
#' dictionary vignette); this function checks that assumption explicitly
#' rather than letting a silent misalignment produce a plausible-looking
#' but wrong result.
#'
#' @param ... Named matrices/data frames to check for consistent subject
#'   sets (matrices are checked by `colnames()`, the `clinical` argument,
#'   if present by name, is checked by `rownames()` or its `id_col`).
#' @param id_col Character, the clinical data frame's ID column name, if
#'   its subject identifiers are a column rather than row names. Default
#'   `"patient_id"`.
#' @return Invisibly, `TRUE` if all objects describe an identical,
#'   identically-ordered subject set; throws an informative error
#'   otherwise (including exactly which subjects are missing from which
#'   object, truncated to the first 10 per object).
#' @export
#' @examples
#' m <- matrix(1, 2, 3, dimnames = list(c("a","b"), c("P1","P2","P3")))
#' e <- matrix(1, 2, 3, dimnames = list(c("g1","g2"), c("P1","P2","P3")))
#' cl <- data.frame(patient_id = c("P1","P2","P3"), asthma = c(0,1,0))
#' validate_cohort_alignment(methylation = m, expression = e, clinical = cl)
validate_cohort_alignment <- function(..., id_col = "patient_id") {
  objs <- list(...)
  if (length(objs) < 2) {
    stop("Supply at least two named objects to check for alignment (e.g. methylation, expression, clinical).",
         call. = FALSE)
  }
  if (is.null(names(objs)) || any(names(objs) == "")) {
    stop("All arguments to validate_cohort_alignment() must be named (e.g. `methylation = m, clinical = cl`).",
         call. = FALSE)
  }

  ids <- lapply(names(objs), function(nm) {
    obj <- objs[[nm]]
    if (is.matrix(obj)) return(colnames(obj))
    if (is.data.frame(obj) && id_col %in% colnames(obj)) return(as.character(obj[[id_col]]))
    if (is.data.frame(obj)) return(rownames(obj))
    stop(sprintf("Don't know how to extract subject IDs from `%s` (class '%s').",
                  nm, paste(class(obj), collapse = "/")), call. = FALSE)
  })
  names(ids) <- names(objs)

  reference_name <- names(ids)[1]
  reference <- ids[[1]]
  problems <- character(0)
  for (nm in names(ids)[-1]) {
    if (!identical(ids[[nm]], reference)) {
      missing_here <- setdiff(reference, ids[[nm]])
      extra_here <- setdiff(ids[[nm]], reference)
      order_differs <- length(missing_here) == 0 && length(extra_here) == 0
      detail <- if (order_differs) {
        sprintf("`%s` has the same subjects as `%s` but in a different order.", nm, reference_name)
      } else {
        sprintf(
          "`%s` vs `%s`: %d subject(s) in `%s` missing from `%s` (%s%s); %d subject(s) in `%s` not in `%s` (%s%s).",
          nm, reference_name, length(missing_here), reference_name, nm,
          paste(utils::head(missing_here, 10), collapse = ", "), if (length(missing_here) > 10) ", ..." else "",
          length(extra_here), nm, reference_name,
          paste(utils::head(extra_here, 10), collapse = ", "), if (length(extra_here) > 10) ", ..." else ""
        )
      }
      problems <- c(problems, detail)
    }
  }
  if (length(problems) > 0) {
    stop("Subject misalignment detected across inputs:\n  ", paste(problems, collapse = "\n  "),
         "\nAll matrices/tables passed to traceAsthma pipeline functions must describe the same ",
         "subjects in the same order. Reindex/reorder before proceeding.", call. = FALSE)
  }
  invisible(TRUE)
}
