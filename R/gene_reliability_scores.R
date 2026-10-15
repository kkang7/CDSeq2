# CDSeq2 gene-cell-type reliability annotations, experimental v0.2.0.
# Base R only. This file never runs CDSeq2 or differential-expression tests.
# Inputs use the saved composition-review MEAN feature schema, not v0.1 medians.

.cgs2_features <- c("contribution_mean", "expected_log_mean_count", "ref_specificity",
                    "bulk_logCPM_sd", "theta_mean", "bulk_max_donor_share", "composition_loo_r2", "bulk_theta_r")
.cgs2_log_features <- c("contribution_mean", "ref_specificity", "bulk_logCPM_sd",
                        "theta_mean", "bulk_max_donor_share")

.cgs2_design <- function(features) {
  if (!is.data.frame(features) || !all(.cgs2_features %in% names(features)))
    stop("Supply a data.frame with the eight composition-review mean features.", call. = FALSE)
  if (!all(vapply(features[, .cgs2_features, drop = FALSE], is.numeric, logical(1))))
    stop("All predictors must be numeric.", call. = FALSE)
  x <- as.matrix(features[, .cgs2_features, drop = FALSE])
  if (any(!is.finite(x))) stop("Predictors must be finite; missing values are not silently imputed.", call. = FALSE)
  for (name in .cgs2_log_features) {
    if (any(x[, name] < 0)) stop(paste(name, "must be nonnegative"), call. = FALSE)
    x[, name] <- log2(x[, name] + 1e-6)
  }
  for (name in c("contribution_mean", "ref_specificity", "theta_mean", "bulk_max_donor_share"))
    if (any(features[[name]] < 0 | features[[name]] > 1 + 1e-10))
      stop(paste(name, "must lie in [0,1]"), call. = FALSE)
  if (any(abs(features$bulk_theta_r) > 1 + 1e-10))
    stop("bulk_theta_r must lie in [-1,1]", call. = FALSE)
  x
}

#' Read the bundled development calibration (plain R expression; base R only).
#' Read only a trusted calibration file, because dget evaluates R expressions.
#'
#' @param path Character scalar giving the path to a trusted calibration file
#'   saved as an R expression, such as \code{development_model.dput}.
#' @return A calibration list containing the model version, context, gene panel,
#'   endpoint definitions, and fitted reliability models.
load_gene_score_model <- function(path) {
  model <- dget(path)
  if (!is.list(model) || !identical(model$version, "0.2.0") ||
      is.null(model$models) || is.null(model$context))
    stop("Unsupported calibration object", call. = FALSE)
  model
}

.cgs2_predict <- function(features, model) {
  if (!identical(model$features, .cgs2_features)) stop("Feature schema mismatch", call. = FALSE)
  x <- .cgs2_design(features)
  z <- sweep(sweep(x, 2, model$center, "-"), 2, model$scale, "/")
  value <- as.vector(model$intercept + z %*% model$coef)
  outside <- rowSums(sweep(x, 2, model$lower - 1e-9, "<") |
                       sweep(x, 2, model$upper + 1e-9, ">")) > 0
  value <- switch(model$transform,
                  fisher = tanh(value), log1p = expm1(pmin(30, pmax(0, value))),
                  bounded = pmin(1, pmax(0, value)), stop("Unknown target transform", call. = FALSE))
  list(value = value, outside = outside)
}

#' Annotate expression-variation and DE-effect reliability for gene-type pairs.
#'
#' @param features data.frame, one row per gene-type-scenario; contains gene,
#' cell_type, cohort_n and the eight exact saved-report columns listed above.
#' Duplicate gene-type keys are allowed across caller-defined scenarios; input
#' row order is preserved. Scenario labels are never score predictors.
#' @param model Object returned by load_gene_score_model().
#' @param context Named list specifying the actual experiment context. Every
#' calibration context field must match. See model$context and README.md.
#' @param allow_extrapolation Default FALSE returns NA when context, gene panel,
#' or observed feature ranges are unsupported. TRUE permits explicitly flagged
#' exploratory predictions; it cannot manufacture an unsupported endpoint.
#'
#' @return data.frame of predictions, separate 0-100 score indices, and status
#' columns. No row is removed. Scores are neither probabilities nor p/q-values.
#' Group-effect score predicts absolute log2CPM contrast error; it does not
#' estimate the chance that a gene is differentially expressed. The injection
#' head is trained only on moderate realized Astro/Oligo injections. Recipient
#' leakage estimates are conditional on the named experimentally injected source.
#' All bundled calibration models remain DEVELOPMENT-ONLY.
score_gene_celltype_v2 <- function(features, model, context, allow_extrapolation = FALSE) {
  if (!identical(model$version, "0.2.0")) stop("Unsupported model version", call. = FALSE)
  if (!all(c("gene", "cell_type", "cohort_n") %in% names(features)))
    stop("gene, cell_type and cohort_n are required", call. = FALSE)
  if (anyNA(features$gene) || anyNA(features$cell_type) || anyNA(features$cohort_n))
    stop("Missing identifiers or cohort_n", call. = FALSE)
  if (!is.list(context)) stop("context must be an explicit named list", call. = FALSE)
  if (length(allow_extrapolation) != 1 || is.na(allow_extrapolation) || !is.logical(allow_extrapolation))
    stop("allow_extrapolation must be TRUE or FALSE", call. = FALSE)
  .cgs2_design(features)
  mismatch <- any(vapply(names(model$context), function(k)
    !isTRUE(all.equal(context[[k]], model$context[[k]], check.attributes = FALSE)), logical(1)))
  out <- features[, c("gene", "cell_type"), drop = FALSE]
  for (head in names(model$head_definitions)) {
    out[[head]] <- rep(NA_real_, nrow(out))
    out[[paste0(head, "_status")]] <- rep("unsupported_cell_type", nrow(out))
  }
  for (key in names(model$models)) {
    bits <- strsplit(key, "::", fixed = TRUE)[[1]]
    head <- bits[1]; ct <- bits[2]
    idx <- which(as.character(features$cell_type) == ct)
    if (!length(idx)) next
    pred <- .cgs2_predict(features[idx, , drop = FALSE], model$models[[key]])
    bad <- pred$outside | mismatch | features$cohort_n[idx] != model$context$cohort_n |
      !(features$gene[idx] %in% model$gene_panel)
    if (!allow_extrapolation) pred$value[bad] <- NA_real_
    out[[head]][idx] <- pred$value
    out[[paste0(head, "_status")]][idx] <- ifelse(bad,
                                                  if (allow_extrapolation) "extrapolation" else "unsupported_context_or_features", "development_only")
  }
  out$variation_score <- 100 * pmin(1, pmax(0, out$predicted_r_normalized))
  out$residual_variation_score <- 100 * pmin(1, pmax(0, out$predicted_r_residual))
  out$deg_effect_score <- 100 / (1 + out$predicted_effect_absolute_error / 0.5)
  out$injection_fidelity_score <- 100 * out$predicted_injection_fidelity
  out$deg_detection_probability <- NA_real_
  out$deg_detection_status <- "not_estimated_sparse_positive_calls"
  out
}

#' Score system 1: cross-sample variation and residual-scale distortion.
#'
#' @param features A data frame containing \code{gene}, \code{cell_type},
#'   \code{cohort_n}, and the eight numeric composition-review mean features
#'   documented in \code{score_gene_celltype_v2}.
#' @param model A calibration list returned by \code{load_gene_score_model}.
#' @param context A named list describing the actual application context.
#'   Its fields are checked against \code{model$context}.
#' @param allow_extrapolation Logical scalar. If \code{FALSE}, unsupported
#'   contexts, genes, or feature ranges receive missing scores. If \code{TRUE},
#'   exploratory predictions are returned with an extrapolation status.
#' @return A data frame preserving input row order, with gene and cell-type
#'   identifiers, variation scores, predicted correlations, residual-scale
#'   distortion, and endpoint-specific status columns.
score_variation_reliability <- function(features, model, context, allow_extrapolation = FALSE) {
  out <- score_gene_celltype_v2(features, model, context, allow_extrapolation)
  keep <- names(out) %in% c("gene", "cell_type") |
    grepl("variation|predicted_r_|sd_distortion", names(out))
  out[, keep, drop = FALSE]
}

#' Score system 2: group-effect error, injection fidelity and recipient leakage.
#'
#' @param features A data frame containing \code{gene}, \code{cell_type},
#'   \code{cohort_n}, and the eight numeric composition-review mean features
#'   documented in \code{score_gene_celltype_v2}.
#' @param model A calibration list returned by \code{load_gene_score_model}.
#' @param context A named list describing the actual application context.
#'   Its fields are checked against \code{model$context}.
#' @param allow_extrapolation Logical scalar. If \code{FALSE}, unsupported
#'   contexts, genes, or feature ranges receive missing scores. If \code{TRUE},
#'   exploratory predictions are returned with an extrapolation status.
#' @return A data frame preserving input row order, with gene and cell-type
#'   identifiers, DE-effect reliability scores, injection fidelity, recipient
#'   leakage predictions, and status columns. The detection probability is
#'   missing because that endpoint has not been calibrated.
score_deg_reliability <- function(features, model, context, allow_extrapolation = FALSE) {
  out <- score_gene_celltype_v2(features, model, context, allow_extrapolation)
  keep <- names(out) %in% c("gene", "cell_type") |
    grepl("deg_|effect_|injection_|leakage_", names(out))
  out[, keep, drop = FALSE]
}

#' Fixed-penalty continuous calibration building block for future development.
#' This fits a reliability regression only; it never fits CDSeq2.
#' Training and validation groups must be split by the caller before use.
#'
#' @param features A training data frame containing the eight numeric
#'   composition-review mean predictors documented in
#'   \code{score_gene_celltype_v2}. Only training rows should be supplied.
#' @param transformed_target A finite numeric vector with one entry per row of
#'   \code{features}, already transformed for the chosen endpoint. For example,
#'   use Fisher-transformed correlations or \code{log1p} absolute errors.
#' @param penalty A positive finite numeric scalar giving the ridge penalty
#'   on the sum of squared slopes in the mean-squared-error objective.
#'   The intercept is not penalized. Defaults to \code{0.1}.
#' @return A list containing feature names, training centers and scales,
#'   coefficients, intercept, penalty, feature-range bounds, and training row
#'   count. This is a calibration building block; endpoint transformation and
#'   context metadata must be supplied when assembling a full scoring model.
fit_reliability_ridge <- function(features, transformed_target, penalty = 0.1) {
  x <- .cgs2_design(features); y <- as.numeric(transformed_target)
  if (length(y) != nrow(x) || any(!is.finite(y)) || !length(y) ||
      length(penalty) != 1 || !is.finite(penalty) || penalty <= 0)
    stop("Invalid calibration target or penalty", call. = FALSE)
  center <- colMeans(x)
  scale <- sqrt(colMeans(sweep(x, 2, center, "-")^2))
  scale[scale < 1e-12] <- 1
  z <- sweep(sweep(x, 2, center, "-"), 2, scale, "/")
  intercept <- mean(y)
  coef <- solve(crossprod(z) / nrow(x) + diag(penalty, ncol(x)),
                crossprod(z, y - intercept) / nrow(x))
  list(features = .cgs2_features, center = unname(center), scale = unname(scale),
       coef = as.vector(coef), intercept = intercept, penalty = penalty,
       lower = apply(x, 2, min), upper = apply(x, 2, max), n_rows = nrow(x))
}
