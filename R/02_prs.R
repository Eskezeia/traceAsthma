#' Compute a trait-specific polygenic risk score
#'
#' \deqn{PRS_{i,t} = \sum_{j=1}^{m_t} \beta_{j,t} G_{ij}}
#'
#' where `G` is genotype dosage and `beta` is a posterior SNP effect
#' estimate for trait `t` (e.g. from PRS-CSx or SBayesRC). This function
#' performs the dosage-weighted sum; obtaining posterior effect sizes is
#' expected to be done upstream with PRS-CSx / SBayesRC command-line tools
#' (see vignette), since these require external LD reference panels well
#' beyond the scope of an R package.
#'
#' @param dosage_matrix Numeric matrix, subjects x variants, of imputed
#'   genotype dosages (0-2).
#' @param effect_sizes Named numeric vector of posterior SNP effect sizes;
#'   names must match `colnames(dosage_matrix)`.
#' @param ancestry Optional factor/character vector (length = nrow) giving
#'   each subject's ancestry group, used to standardize the score within
#'   stratum (recommended for cross-ancestry cohorts).
#' @return Named numeric vector of per-subject trait-specific PRS
#'   (standardized).
#' @export
#' @examples
#' dosage <- matrix(c(0,1,2,1, 2,0,1,1), nrow = 4,
#'                   dimnames = list(paste0("S", 1:4), c("rs1","rs2")))
#' eff <- c(rs1 = 0.12, rs2 = -0.08)
#' compute_prs(dosage, eff)
compute_prs <- function(dosage_matrix, effect_sizes, ancestry = NULL) {
  raw <- weighted_score(dosage_matrix, effect_sizes, standardize = FALSE)
  z_score(raw, strata = ancestry)
}

#' Construct the multi-trait polygenic risk score (MPRS)
#'
#' \deqn{MPRS_i = \sum_t w_t PRS_{i,t}}
#'
#' Combines trait-specific PRS (asthma, IgE, eosinophil count, FEV1,
#' FEV1/FVC, or any user-supplied set) into a single composite score via
#' elastic-net regression against asthma status. Weights must be estimated
#' on a training split only -- see `train_idx` -- to avoid information
#' leakage into later model evaluation.
#'
#' @param prs_matrix Numeric matrix or data frame, subjects x trait-specific
#'   PRS (e.g. columns `PRS_Asthma`, `PRS_IgE`, `PRS_EOS`, `PRS_FEV1`,
#'   `PRS_FEV1_FVC`), typically the output of repeated [compute_prs()]
#'   calls bound with `cbind()`.
#' @param asthma_status Binary (0/1) outcome vector, same length as
#'   `nrow(prs_matrix)`.
#' @param train_idx Integer or logical vector indexing which subjects to
#'   use for weight estimation (e.g. the training fold in a
#'   cross-validation loop). Defaults to all subjects, but you should pass
#'   a training-only index when calling this inside [nested_cv()].
#' @param alpha Elastic-net mixing parameter. Default 0.5.
#' @return A list with `mprs` (named numeric vector, all subjects,
#'   standardized), `weights` (named numeric vector of trait weights), and
#'   `fit` (the underlying `glmnet::cv.glmnet` fit object, for inspection).
#' @export
compute_mprs <- function(prs_matrix, asthma_status, train_idx = NULL, alpha = 0.5) {
  requirePkg("glmnet", "multi-trait PRS weight estimation")
  prs_matrix <- as.matrix(prs_matrix)
  if (is.null(train_idx)) train_idx <- seq_len(nrow(prs_matrix))

  fit <- glmnet::cv.glmnet(
    prs_matrix[train_idx, , drop = FALSE],
    asthma_status[train_idx],
    family = "binomial", alpha = alpha
  )
  coefs <- as.matrix(stats::coef(fit, s = "lambda.min"))
  weights <- coefs[rownames(coefs) != "(Intercept)", 1]
  weights <- weights[colnames(prs_matrix)]  # keep zero-weighted traits, order-matched
  names(weights) <- colnames(prs_matrix)

  mprs <- weighted_score(prs_matrix, weights, standardize = TRUE)
  list(mprs = mprs, weights = weights, fit = fit)
}
