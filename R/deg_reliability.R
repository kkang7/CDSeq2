# Package interface 0.23.1; frozen statistical models are unchanged from v0.23.
# No fitting, DE tests, package attachment, workers, or filesystem writes here.

.cdsdeg_assert <- function(ok, ...) {
  if (!isTRUE(ok)) stop(..., call. = FALSE)
}
.cdsdeg_key <- function(d) paste(d$contrast, d$gene, d$cell_type, sep = "\r")
.cdsdeg_ids <- function(d, label) {
  .cdsdeg_assert(is.data.frame(d) && !anyDuplicated(names(d)), label, " must be a data.frame with unique column names.")
  .cdsdeg_assert(all(c("contrast", "gene", "cell_type") %in% names(d)), label, " needs contrast, gene, cell_type.")
  for (nm in c("contrast", "gene", "cell_type")) {
    d[[nm]] <- as.character(d[[nm]])
    .cdsdeg_assert(!anyNA(d[[nm]]) && all(nzchar(d[[nm]])) &&
      !any(grepl("\r", d[[nm]], fixed = TRUE)), "Invalid ", label, "$", nm)
  }
  .cdsdeg_assert(!anyDuplicated(.cdsdeg_key(d)), label, " has duplicated contrast/gene/cell_type rows.")
  d
}
.cdsdeg_de <- function(de) {
  .cdsdeg_assert(is.data.frame(de), "de must be a data.frame.")
  if (!"contrast" %in% names(de)) de$contrast <- if ("dataset" %in% names(de)) as.character(de$dataset) else rep("contrast1", nrow(de))
  de <- .cdsdeg_ids(de, "de")
  required <- c("estimated_effect", "p_value", "q_global")
  .cdsdeg_assert(nrow(de) > 0L && all(required %in% names(de)), "de needs estimated_effect, p_value, q_global and at least one tested row.")
  for (nm in required) .cdsdeg_assert(is.numeric(de[[nm]]) && all(is.finite(de[[nm]])),
    "de$", nm, " must be finite numeric values. Supply the eligible tested universe; exclude untested rows.")
  .cdsdeg_assert(all(de$p_value >= 0 & de$p_value <= 1 & de$q_global >= 0 & de$q_global <= 1), "p_value and q_global must lie in [0, 1].")
  .cdsdeg_assert(all(de$q_global + 1e-12 >= de$p_value), "q_global cannot be smaller than the original p_value.")
  de[, c("contrast", "gene", "cell_type", required), drop = FALSE]
}
.cdsdeg_cross <- function(de) {
  n <- nrow(de)
  out <- data.frame(estimated_abs_effect = log1p(abs(de$estimated_effect)),
    estimated_sign = sign(de$estimated_effect),
    estimated_log_p_support = pmin(50, -log10(pmax(de$p_value, 1e-50))),
    other_max_abs_effect = rep(NA_real_, n), other_mean_abs_effect = NA_real_,
    self_abs_share = NA_real_, other_same_sign_fraction = NA_real_, other_effect_count = 0L)
  for (ii in split(seq_len(n), paste(de$contrast, de$gene, sep = "\r"))) for (i in ii) {
    other <- de$estimated_effect[setdiff(ii, i)]
    out$other_effect_count[i] <- length(other)
    if (!length(other)) next
    out$other_max_abs_effect[i] <- log1p(max(abs(other)))
    out$other_mean_abs_effect[i] <- log1p(mean(abs(other)))
    out$self_abs_share[i] <- abs(de$estimated_effect[i]) / (abs(de$estimated_effect[i]) + sum(abs(other)) + 1e-8)
    out$other_same_sign_fraction[i] <- mean(sign(other) == sign(de$estimated_effect[i]))
  }
  out$cross_type_available <- out$other_effect_count > 0L
  out
}
.cdsdeg_equal <- function(a, b, tol = 1e-8) {
  length(a) == length(b) && identical(is.na(a), is.na(b)) && all(abs(a[!is.na(a)] - b[!is.na(b)]) <= tol)
}
.cdsdeg_features <- function(de, features, model) {
  .cdsdeg_assert(is.data.frame(features), "features must be a data.frame from prepare_deg_features().")
  if (!"contrast" %in% names(features)) {
    if ("dataset" %in% names(features)) features$contrast <- as.character(features$dataset)
    else {
      .cdsdeg_assert(length(unique(de$contrast)) == 1L, "Multiple contrasts require explicit features$contrast; positional alignment is not accepted.")
      features$contrast <- rep(de$contrast[1L], nrow(features))
    }
  }
  features <- .cdsdeg_ids(features, "features")
  .cdsdeg_assert(setequal(.cdsdeg_key(de), .cdsdeg_key(features)), "de and features keys must match exactly.")
  features <- features[match(.cdsdeg_key(de), .cdsdeg_key(features)), , drop = FALSE]
  rownames(features) <- NULL
  need <- unique(c(model$core$own, model$core$cross))
  .cdsdeg_assert(all(c(need, "feature_schema", "cross_type_available") %in% names(features)), "Missing prepared predictor columns.")
  .cdsdeg_assert(!anyNA(features$feature_schema) && all(features$feature_schema == model$feature_schema), "Feature schema mismatch; raw observables are not prepared predictors.")
  .cdsdeg_assert(all(vapply(features[need], is.numeric, TRUE)) && !any(is.infinite(as.matrix(features[need]))), "Prepared features must be numeric; use NA, not infinity, for missing support.")
  .cdsdeg_assert(is.logical(features$cross_type_available) && !anyNA(features$cross_type_available), "cross_type_available must be logical with no NA.")
  derived <- .cdsdeg_cross(de)
  for (nm in names(derived)) .cdsdeg_assert(.cdsdeg_equal(as.numeric(features[[nm]]), as.numeric(derived[[nm]])),
    "Predictor ", nm, " disagrees with de. Include every eligible tested cell type for each gene and regenerate features.")
  features
}

#' Prepare predictors for the frozen CDSeq2 DEG scorer
#'
#' Computes the unchanged v0.4 feature definitions used by the v0.23 validation.
#' Supplies transformed predictors, not a reliability score. No DE test is run.
#' @param de Data frame of all eligible tests for ONE contrast, containing gene,
#'   cell_type, estimated_effect, p_value, q_global. Optional contrast identifier.
#'   Preserve original global BH values; do not filter to significant calls.
#' @param bulk Numeric genes-by-samples observed count matrix, with dimnames.
#' @param reference_mean Numeric genes-by-cell-types shared reference profile;
#'   columns sum to one over the SAME full modeled gene panel as bulk.
#' @param rna_proportions Numeric samples-by-cell-types estimated RNA fractions;
#'   rows sum to one. These are RNA fractions, not cell-number fractions.
#' @return Data frame keyed by contrast/gene/cell_type in de order, with the
#'   transformed frozen predictors and feature_schema='deg_features_v0.4.0'.
#' @details At least K+2 samples and a full-rank composition design are required.
#'   Matrix axes are aligned by names; profiles are never silently renormalized.
#'   Features pool the samples in the contrast, matching the frozen model.
#' @export
#' @examples
#' data(cdseq2_toy)
#' x <- cdseq2_toy
#' de <- expand.grid(gene = rownames(x$data$reference_mean),
#'     cell_type = colnames(x$data$reference_mean), stringsAsFactors = FALSE)
#' de$estimated_effect <- 0.1
#' de$p_value <- 0.5
#' de$q_global <- 0.5
#' features <- prepare_deg_features(de, x$data$validation$bulk,
#'     x$data$reference_mean, x$validation_fractions)
#' head(features)
prepare_deg_features <- function(de, bulk, reference_mean, rna_proportions) {
  de <- .cdsdeg_de(de)
  .cdsdeg_assert(length(unique(de$contrast)) == 1L, "Prepare each contrast separately using its own sample subset, then rbind the feature tables.")
  o <- .cdsdeg_observables(bulk, reference_mean, rna_proportions)
  d <- de; d$scenario <- de$contrast; d$called <- de$q_global <= .05
  o$scenario <- de$contrast[1L]
  f <- .cdsdeg_prepare_canonical(o, d, list(feature_schema = "pilot_logCPM_r",
    effect_scale = "log2CPM_contrast", call_definition = "supplied_original_global_BH"))
  f$contrast <- f$scenario; f$scenario <- NULL; f$called <- NULL
  f[, c("contrast", "gene", "cell_type", setdiff(names(f), c("contrast", "gene", "cell_type"))), drop = FALSE]
}

.cdsdeg_core <- function(de, x, core) {
  n <- nrow(de)
  out <- data.frame(deg_score = rep(NA_real_, n), score_status = ifelse(de$q_global <= .1, "missing_own_predictors", "outside_scoring_universe"),
    model_variant = rep(NA_character_, n), model_family = NA_character_, model_id = NA_character_, clipped_features = rep("", n))
  own <- rowSums(!is.finite(as.matrix(x[, core$own, drop = FALSE]))) == 0L
  cross <- x$cross_type_available & rowSums(!is.finite(as.matrix(x[, core$cross, drop = FALSE]))) == 0L
  types <- unique(sub("::.*", "", names(core$reference)))
  out$score_status[de$q_global <= .1 & !de$cell_type %in% types] <- "unsupported_cell_type"
  for (ct in intersect(unique(de$cell_type), types)) for (cr in c(FALSE, TRUE)) {
    ii <- which(de$q_global <= .1 & de$cell_type == ct & own & cross == cr)
    if (!length(ii)) next
    variants <- if (cr) c("cross_type", "own_type") else "own_type"
    id <- paste(ct, variants, sep = "::"); id <- id[id %in% names(core$reference)][1L]
    .cdsdeg_assert(!is.na(id), "Frozen model route is missing for ", ct)
    for (fam in c("reference", "expanded")) {
      jj <- ii[if (fam == "reference") de$q_global[ii] <= .05 else de$q_global[ii] > .05]
      if (!length(jj)) next
      m <- core[[fam]][[id]]; .cdsdeg_assert(!is.null(m), "Frozen model missing: ", fam, " / ", id)
      pp <- m$preprocessing; raw <- as.matrix(x[jj, pp$features, drop = FALSE])
      outside <- sweep(raw, 2L, pp$lower - 1e-9, "<") | sweep(raw, 2L, pp$upper + 1e-9, ">")
      flags <- apply(outside, 1L, function(z) paste(pp$features[z], collapse = ";"))
      for (j in seq_len(ncol(raw))) raw[, j] <- pmin(pp$upper[j], pmax(pp$lower[j], raw[, j]))
      nd <- as.data.frame(sweep(sweep(raw, 2L, pp$center, "-"), 2L, pp$scale, "/"))
      names(nd) <- unname(m$variable_map)
      p <- as.numeric(stats::predict(m$fit, newdata = nd, type = "response"))
      .cdsdeg_assert(all(is.finite(p) & p >= 0 & p <= 1), "Invalid frozen model prediction.")
      out$deg_score[jj] <- p; out$clipped_features[jj] <- flags
      out$score_status[jj] <- ifelse(nzchar(flags), "scored_clipped_to_training_ranges", "scored_within_training_ranges")
      out$model_variant[jj] <- sub(".*::", "", id); out$model_family[jj] <- fam; out$model_id[jj] <- paste(fam, id, sep = "::")
    }
  }
  out
}
.cdsdeg_annotation_design <- function(de, x, types, own) {
  key <- paste(de$contrast, de$gene, sep = "\r"); keys <- unique(key)
  gi <- match(key, keys); ix <- split(seq_len(nrow(de)), factor(gi, levels = seq_along(keys)))
  rx <- x[, own, drop = FALSE]; rx$signed_effect <- sign(de$estimated_effect) * log1p(abs(de$estimated_effect))
  gx <- data.frame(max_abs_effect = vapply(ix, function(ii) max(log1p(abs(de$estimated_effect[ii]))), 0),
    max_p_support = vapply(ix, function(ii) max(pmin(50, -log10(pmax(de$p_value[ii], 1e-50)))), 0),
    n_q10 = vapply(ix, function(ii) sum(de$q_global[ii] <= .1), 0L), n_eligible = lengths(ix))
  for (ct in types) {
    ii <- which(de$cell_type == ct); jj <- match(keys, key[ii]); valid <- !is.na(jj); rr <- ii[jj[valid]]
    values <- list(effect = rx$signed_effect, psupport = rx$estimated_log_p_support,
      contribution = rx$contribution_mean, specificity = rx$ref_specificity, count = rx$expected_log_mean_count)
    for (f in names(values)) {
      nm <- paste0("ct_", make.names(ct), "_", f); v <- rep(NA_real_, length(keys)); v[valid] <- values[[f]][rr]
      gx[[nm]] <- v; rx[[nm]] <- v[gi]
    }
    nm <- paste0("ct_", make.names(ct), "_eligible"); gx[[nm]] <- as.integer(valid); rx[[nm]] <- as.integer(valid)[gi]
  }
  list(x = rx, gx = gx, gi = gi, ix = ix, trigger = gx$n_q10 > 0L)
}
.cdsdeg_annotation_predict <- function(m, x, ii) {
  pp <- m$pp; raw <- as.matrix(x[ii, pp$input, drop = FALSE]); miss <- !is.finite(raw)
  for (j in seq_along(pp$input)) raw[miss[, j], j] <- pp$median[j]
  if (length(pp$missing)) {
    v <- 1L * miss[, pp$missing, drop = FALSE]; colnames(v) <- paste0(pp$missing, "__missing"); raw <- cbind(raw, v)
  }
  raw <- raw[, pp$features, drop = FALSE]
  clipped <- rep(FALSE, length(ii))
  if (ncol(raw)) {
    clipped <- rowSums(sweep(raw, 2L, pp$lower - 1e-9, "<") | sweep(raw, 2L, pp$upper + 1e-9, ">")) > 0L
    for (j in seq_len(ncol(raw))) raw[, j] <- pmin(pp$upper[j], pmax(pp$lower[j], raw[, j]))
    raw <- sweep(sweep(raw, 2L, pp$center, "-"), 2L, pp$scale, "/")
  }
  nd <- as.data.frame(raw); names(nd) <- m$variables
  if (!ncol(raw)) nd <- data.frame(row = seq_along(ii))
  p <- as.numeric(stats::predict(m$fit, newdata = nd, type = "response"))
  .cdsdeg_assert(length(p) == length(ii) && all(is.finite(p) & p >= 0 & p <= 1), "Invalid annotation prediction.")
  list(score = p, clipped = clipped, imputed = rowSums(miss) > 0L)
}
.cdsdeg_annotations <- function(de, x, ann, enabled) {
  n <- nrow(de); empty <- rep(NA_real_, n)
  out <- data.frame(gene_score = empty, gene_threshold = empty, gene_supported = FALSE, gene_pass = FALSE,
    attribution_risk = empty, risk_threshold = empty, risk_threshold_supported = FALSE, flag_high_risk = FALSE,
    source_score = empty, source_threshold = empty, source_threshold_supported = FALSE,
    source_selected = FALSE, candidate_types = rep("", n), alternative_types = rep("", n),
    n_candidate_types = integer(n), annotation_status = rep(if (enabled) "outside_triggered_gene" else "not_requested", n),
    gene_clipped = rep(NA, n), gene_imputed = NA, risk_clipped = NA, risk_imputed = NA,
    source_clipped = NA, source_imputed = NA, stringsAsFactors = FALSE)
  if (!enabled) return(out)
  b <- .cdsdeg_annotation_design(de, x, ann$types, ann$own)
  # Named-type models were trained on these identities only. Do not reinterpret new labels.
  bad_contrasts <- unique(de$contrast[!de$cell_type %in% ann$types])
  supported_context <- !de$contrast %in% bad_contrasts
  gs <- rep(NA_real_, nrow(b$gx)); gclip <- gimp <- rep(NA, length(gs))
  first <- vapply(b$ix, `[`, 1L, 1L)
  for (id in names(ann$models)) {
    q <- strsplit(id, "::", fixed = TRUE)[[1L]]; method <- q[1L]; ct <- q[2L]
    isgene <- method == "Gene_named"
    ii <- if (isgene) which(b$trigger & supported_context[first]) else which(de$cell_type == ct & supported_context &
      if (method == "Risk_own") de$q_global <= .1 else b$trigger[b$gi])
    if (!length(ii)) next
    p <- .cdsdeg_annotation_predict(ann$models[[id]], if (isgene) b$gx else b$x, ii)
    if (isgene) { gs[ii] <- p$score; gclip[ii] <- p$clipped; gimp[ii] <- p$imputed }
    else if (method == "Risk_own") {
      out$attribution_risk[ii] <- p$score; out$risk_clipped[ii] <- p$clipped; out$risk_imputed[ii] <- p$imputed
    } else { out$source_score[ii] <- p$score; out$source_clipped[ii] <- p$clipped; out$source_imputed[ii] <- p$imputed }
  }
  cut <- function(method, ct) {
    t <- ann$thresholds; v <- t$threshold[match(paste(method, ct), paste(t$method, t$cell_type))]
    v[!is.finite(v)] <- NA_real_; v
  }
  out$gene_score <- gs[b$gi]; out$gene_clipped <- gclip[b$gi]; out$gene_imputed <- gimp[b$gi]
  out$gene_threshold <- rep(cut("Gene_named", "GENE"), n)
  out$gene_supported <- is.finite(out$gene_score) & is.finite(out$gene_threshold)
  out$gene_pass <- out$gene_supported & out$gene_score >= out$gene_threshold
  out$risk_threshold <- cut("Risk_own", de$cell_type)
  out$risk_threshold_supported <- is.finite(out$attribution_risk) & is.finite(out$risk_threshold)
  out$flag_high_risk <- out$risk_threshold_supported & out$attribution_risk >= out$risk_threshold
  out$source_threshold <- cut("Source_own", de$cell_type)
  out$source_threshold_supported <- is.finite(out$source_score) & is.finite(out$source_threshold) & out$gene_supported
  eligible <- out$source_threshold_supported & out$gene_pass & out$source_score >= out$source_threshold
  for (ii in b$ix) {
    jj <- ii[eligible[ii]]
    jj <- utils::head(jj[order(-out$source_score[jj], de$cell_type[jj], method = "radix")], 2L)
    out$source_selected[jj] <- TRUE
    out$candidate_types[ii] <- paste(de$cell_type[jj], collapse = ";")
    out$n_candidate_types[ii] <- length(jj)
    for (i in ii) out$alternative_types[i] <- paste(de$cell_type[setdiff(jj, i)], collapse = ";")
  }
  out$annotation_status[b$trigger[b$gi] & supported_context] <- "exploratory_annotations"
  out$annotation_status[!supported_context] <- "unsupported_type_set"
  out
}

.cdsdeg_context <- function(context, expected) {
  if (is.null(context)) return("context_not_declared")
  .cdsdeg_assert(is.list(context) && !is.null(names(context)) && !anyDuplicated(names(context)), "context must be NULL or a named list.")
  ok <- vapply(names(expected), function(nm) length(context[[nm]]) == 1L && !is.na(context[[nm]]) &&
    as.character(context[[nm]]) == as.character(expected[[nm]]), TRUE)
  if (all(ok)) "declared_benchmark_context" else "context_differs_or_incomplete"
}

#' Score the reliability of cell-type-specific DEG calls
#'
#' Uses frozen GAM models validated in v0.23; does not refit CDSeq2 or call DEGs.
#' The main target is a nonzero true effect in the called type with correct sign,
#' among original BH candidates. No effect-size tolerance/cutoff is imposed.
#' @param de All eligible DE tests, including nonsignificant cell types. Required:
#'   gene, cell_type, estimated_effect, p_value, q_global. Optional contrast;
#'   dataset is accepted as its alias, otherwise 'contrast1'. q_global contains
#'   original BH values adjusted over all eligible gene/type tests per contrast.
#' @param features Prepared, transformed predictors from prepare_deg_features().
#'   Keys are matched by names, never by row position. For multiple contrasts,
#'   provide features$contrast. Supply either features or all three matrices.
#' @param bulk,reference_mean,rna_proportions Matrices for one contrast, as
#'   documented in prepare_deg_features(). Not needed when features is supplied.
#' @param model Frozen model from load_deg_scorer(); NULL loads the package model.
#' @param annotations Include exploratory gene, risk, and source annotations.
#' @param bh_threshold Reporting threshold from 0 to 0.10 inclusive; default 0.10. Models
#'   always route at original q=0.05 and abstain above q=0.10.
#' @param score_threshold Reporting cutoff from 0 to 1 inclusive; default 0.75.
#' @param context Optional named provenance list. Compare model$context for the
#'   benchmark conditions. Unspecified/different settings are marked exploratory.
#' @param warn Warn about context differences and unknown cell-type labels.
#' @return A data.frame in de order. deg_score and all annotation scores use 0--1.
#'   retained is the original-call decision; retained_plus_gene additionally
#'   applies the frozen gene gate. score_status distinguishes unsupported/NA
#'   scores from low scores. Attribution flags never override retained. Source
#'   candidates are exploratory and never reassign the original DEG.
#' @details Training used paired limma-trend coefficients on log2(CPM+0.5), with
#'   20 pairs and seven ROSMAP PFC types. Scores are not p-values or FDR guarantees.
#'   Threshold changes affect reporting only, not models or annotations. Annotation
#'   triggers always use q<=0.10, including when bh_threshold is lower.
#'   Vascular/unknown types have no core model. Inputs contain no ground truth.
#' @export
#' @examples
#' model <- load_deg_scorer()
#' # Bundled prediction controls demonstrate the input schema, not new validation.
#' scores <- score_deg(model$self_test$de, features = model$self_test$features,
#'     model = model, context = model$context)
#' head(scores)
score_deg <- function(de, features = NULL, bulk = NULL, reference_mean = NULL,
                      rna_proportions = NULL, model = NULL, annotations = TRUE,
                      bh_threshold = .10, score_threshold = .75,
                      context = NULL, warn = TRUE) {
  de <- .cdsdeg_de(de)
  for (nm in c("annotations", "warn")) .cdsdeg_assert(is.logical(get(nm)) && length(get(nm)) == 1L && !is.na(get(nm)), nm, " must be TRUE or FALSE.")
  .cdsdeg_assert(is.numeric(bh_threshold) && length(bh_threshold) == 1L && is.finite(bh_threshold) && bh_threshold >= 0 && bh_threshold <= .1, "bh_threshold must be between 0 and 0.10.")
  .cdsdeg_assert(is.numeric(score_threshold) && length(score_threshold) == 1L && is.finite(score_threshold) && score_threshold >= 0 && score_threshold <= 1, "score_threshold must be between 0 and 1.")
  if (is.null(model)) model <- load_deg_scorer()
  .cdsdeg_model_valid(model)
  .cdsdeg_assert(requireNamespace("mgcv", quietly = TRUE), "Install mgcv to predict with the frozen GAMs.")
  matrices <- list(bulk, reference_mean, rna_proportions)
  if (is.null(features)) {
    .cdsdeg_assert(all(!vapply(matrices, is.null, TRUE)), "Supply features, or bulk + reference_mean + rna_proportions.")
    features <- prepare_deg_features(de, bulk, reference_mean, rna_proportions)
  } else .cdsdeg_assert(all(vapply(matrices, is.null, TRUE)), "Supply features OR matrices, not both.")
  x <- .cdsdeg_features(de, features, model)
  cs <- .cdsdeg_context(context, model$context)
  if (warn && cs != "declared_benchmark_context") warning("Frozen ROSMAP PFC calibration: supplied context differs or is incomplete. Treat scores as exploratory; see model$context.", call. = FALSE)
  if (warn && any(!de$cell_type %in% model$types)) warning("Unknown cell-type labels abstain; named-type annotations are disabled for affected contrasts. Do not automatically map unrelated types.", call. = FALSE)
  core <- .cdsdeg_core(de, x, model$core)
  ann <- .cdsdeg_annotations(de, x, model$annotations, annotations)
  out <- cbind(de, core, ann)
  out$bh_pass <- de$q_global <= bh_threshold
  out$retained <- out$bh_pass & is.finite(out$deg_score) & out$deg_score >= score_threshold
  out$retained_plus_gene <- out$retained & out$gene_pass
  out$report_status <- ifelse(!out$bh_pass, "outside_reporting_BH_threshold", ifelse(out$retained,
    "original_call_supported", ifelse(out$n_candidate_types > 0L, "candidate_locations",
      ifelse(out$gene_pass, "location_unresolved", ifelse(!annotations, "original_call_not_retained", "insufficient_gene_evidence_or_support")))))
  out$risk_state <- ifelse(!out$risk_threshold_supported, "unavailable", ifelse(out$flag_high_risk, "high_attribution_risk", "below_warning_threshold"))
  out$interpretation <- ifelse(out$report_status == "original_call_supported",
    "Original call passes the score cutoff; attribution warnings and locations are supplementary.",
    ifelse(out$report_status == "candidate_locations", paste0("Exploratory change location: ", gsub(";", " and/or ", out$candidate_types), ". Alternative-type significance and direction are not established."),
      ifelse(out$report_status == "location_unresolved", "Gene-level evidence passes; location is unresolved.",
        ifelse(!out$bh_pass, "Original BH q exceeds the reporting threshold.", "Original call not retained; inspect score and support status."))))
  out$model_version <- model$model_version; out$interface_version <- "0.23.1"; out$context_status <- cs
  out$bh_threshold <- bh_threshold; out$score_threshold <- score_threshold
  rownames(out) <- NULL
  attr(out, "scoring_provenance") <- list(model_version = model$model_version, context = context,
    training_context = model$context, source_checksums = model$source_checksums,
    score_scale = "0-1", annotation_trigger_q = .1, gene_gate_optional = TRUE)
  out
}

.cdsdeg_model_valid <- function(model) {
  .cdsdeg_assert(is.list(model) && identical(model$format, "CDSeq2_DEG_frozen_1") &&
    identical(model$feature_schema, "deg_features_v0.4.0") && identical(model$model_version, "v0.23_core_v0.18"),
    "Incompatible model bundle; use the supplied frozen model RDS.")
  .cdsdeg_assert(length(model$core$reference) == 12L && length(model$core$expanded) == 12L && length(model$annotations$models) == 13L, "Incomplete model bundle.")
  invisible(model)
}

#' Load the frozen DEG model bundle
#' @param path Optional RDS path to a complete frozen model. NULL loads the bundled
#' model from its metadata, core and annotation files in extdata.
#'   When sourcing files before package installation, supply this path explicitly.
#' @param check Run bundled numerical prediction checks on the current mgcv build.
#' @return Frozen model list for score_deg(). No training or file changes occur.
#' @export
#' @examples
#' model <- load_deg_scorer()
#' model$context
load_deg_scorer <- function(path = NULL, check = TRUE) {
  .cdsdeg_assert(is.logical(check) && length(check) == 1L && !is.na(check), "check must be TRUE or FALSE.")
  # Assisted-by: OpenAI Codex. Split storage preserves the frozen model.
  if (is.null(path)) {
    root <- system.file("extdata", package = "CDSeq2", mustWork = TRUE)
    model <- readRDS(file.path(root, "deg_scorer_metadata.rds"))
    model$core <- readRDS(file.path(root, "deg_scorer_core.rds"))
    model$annotations <- readRDS(file.path(root, "deg_scorer_annotations.rds"))
  } else {
    .cdsdeg_assert(is.character(path) && length(path) == 1L &&
      !is.na(path) && nzchar(path) && file.exists(path),
      "Supply a readable RDS path to a complete frozen DEG model.")
    model <- readRDS(path)
  }
  .cdsdeg_model_valid(model)
  if (check) check_deg_scorer(model)
  model
}

#' Check frozen DEG predictions against bundled numerical controls
#' @param model Model list returned by load_deg_scorer(check=FALSE).
#' @param tolerance Maximum absolute difference on the 0--1 scale (default 1e-7).
#' @return Invisible data.frame of checks; stops on mismatched scores/decisions.
#' @details This checks current-runtime prediction compatibility, not biological
#'   generalization. Controls were recorded from the frozen validation workflow.
#' @export
#' @examples
#' model <- load_deg_scorer(check = FALSE)
#' check_deg_scorer(model)
check_deg_scorer <- function(model, tolerance = 1e-7) {
  .cdsdeg_model_valid(model)
  .cdsdeg_assert(is.numeric(tolerance) && length(tolerance) == 1L && is.finite(tolerance) && tolerance > 0, "Invalid tolerance.")
  s <- model$self_test
  .cdsdeg_assert(is.list(s) && all(c("de", "features", "expected") %in% names(s)), "Model bundle lacks prediction controls.")
  r <- score_deg(s$de, features = s$features, model = model, context = model$context, warn = FALSE)
  rows <- lapply(names(s$expected), function(nm) {
    a <- r[[nm]]; b <- s$expected[[nm]]
    ok <- if (is.numeric(b)) .cdsdeg_equal(a, b, tolerance) else identical(a, b)
    delta <- if (is.numeric(b) && any(is.finite(b))) max(abs(a[is.finite(b)] - b[is.finite(b)])) else NA_real_
    data.frame(check = nm, passed = ok, max_abs_difference = delta)
  })
  out <- do.call(rbind, rows)
  .cdsdeg_assert(all(out$passed), "Frozen prediction compatibility check failed: ", paste(out$check[!out$passed], collapse = ", "), ". Keep the supplied model/code together; inspect mgcv compatibility.")
  invisible(out)
}
