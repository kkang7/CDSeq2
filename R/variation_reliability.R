# Application-specific GAM6/GAM7 training and prediction interface, revision 2.4.0.
# No automatic pretrained model; no CDSeq2 fitting or p-value recalculation.
.vscore_features <- c('contribution_mean', 'expected_log_mean_count',
  'ref_specificity', 'bulk_logCPM_sd', 'bulk_max_donor_share',
  'composition_loo_r2', 'bulk_theta_r')
.vscore_log_features <- c('contribution_mean', 'ref_specificity',
  'bulk_logCPM_sd', 'bulk_max_donor_share')
.vscore_assert <- function(ok, ...) {
  if (!isTRUE(ok)) stop(..., call. = FALSE)
}
.vscore_ids <- function(x, what) {
  .vscore_assert(!is.null(x) && length(x) > 0L && !anyNA(x) &&
    all(nzchar(x)) && !anyDuplicated(x), what, ' must have unique nonempty names.')
}
.vscore_matrix <- function(x, what) {
  x <- as.matrix(x)
  .vscore_assert(is.numeric(x) && length(dim(x)) == 2L && all(dim(x) > 0L),
    what, ' must be a nonempty numeric matrix.')
  .vscore_ids(rownames(x), paste(what, 'rows'))
  .vscore_ids(colnames(x), paste(what, 'columns'))
  .vscore_assert(all(is.finite(x)) && all(x >= 0),
    what, ' must contain finite nonnegative values.')
  x
}
.vscore_cor <- function(a, b) {
  if (length(a) < 3L || any(!is.finite(c(a, b))) ||
      stats::sd(a) < 1e-12 || stats::sd(b) < 1e-12) return(NA_real_)
  unname(stats::cor(a, b))
}
.vscore_x <- function(d, features = .vscore_features) {
  x <- as.matrix(d[, features, drop = FALSE]); storage.mode(x) <- 'double'
  for (nm in .vscore_log_features) x[, nm] <- log2(x[, nm] + 1e-6)
  x
}
.vscore_key <- function(d) paste(d$gene, d$cell_type, sep = '\r')

# The historical seven-feature definition remains unchanged. Model selection
# is explicit and local; no caller or package-global feature vector is mutated.
.vscore_spec <- function(model = c('GAM6', 'GAM7')) {
  model <- match.arg(model)
  list(name = model,
    features = if (model == 'GAM6') setdiff(.vscore_features, 'composition_loo_r2') else .vscore_features,
    logs = .vscore_log_features,
    schema = paste0('CDSeq2_variation_', model, '_v1'),
    method = if (model == 'GAM6') 'GAM6_no_composition_loo_r2' else 'GAM7_no_theta')
}
.vscore_model_spec <- function(scorer) {
  .vscore_assert(is.list(scorer), 'Scorer must be a list.')
  candidates <- lapply(c('GAM6', 'GAM7'), .vscore_spec)
  ok <- vapply(candidates, function(z) identical(scorer$schema, z$schema) &&
    identical(scorer$method, z$method) && identical(scorer$features, z$features), logical(1))
  .vscore_assert(sum(ok) == 1L,
    'Unsupported scorer schema. Use train_variation_scorer(); for the saved GAM6 ablation use import_variation_scorer().')
  z <- candidates[[which(ok)]]
  .vscore_assert(identical(scorer$log_features, z$logs) ||
    (z$name == 'GAM7' && is.null(scorer$log_features)), 'Scorer transformation metadata differs.')
  z
}
.vscore_join <- function(data, predictions) {
  cbind(data, predictions[setdiff(names(predictions), names(data))])
}
.vscore_sentinel <- function(features, scorer) {
  p <- .vscore_predict_table(features, scorer)
  list(features = features, predictions = p[c('gene', 'cell_type', 'predicted_r',
    'variation_score', 'outside_training_range', 'range_features')])
}


#' Prepare phenotype-independent variation reliability features
#'
#' Computes the unchanged seven observable features over the complete modeled gene
#' panel. Names are aligned explicitly; no genes are silently intersected.
#' Features use the supplied reference mean and estimated RNA proportions,
#' not CDSeq2 posterior allocation variability. No phenotype is required.
#' GAM6 uses six features; composition_loo_r2 is returned as a diagnostic only.
#' bulk_theta_r remains the correlation of bulk CPM (not logCPM) with RNA fractions.
#'
#' @param bulk Numeric genes-by-samples matrix of observed bulk counts, with
#'   row and column names. Every sample must have a positive panel library.
#' @param reference_mean Numeric genes-by-cell-types matrix, normalized to sum
#'   to one over genes in each type. Use the same support-reference definition as training. For the ROSMAP
#'   example this is the unsmoothed pooled-reference mean, distinct from the
#'   smoothed CDSeq2 fitting prior; do not substitute estimated GEPs.
#' @param rna_fractions Numeric samples-by-cell-types matrix of estimated RNA
#'   fractions, with rows summing to one. These are not cell-number fractions.
#' @param scorer Optional object from [load_variation_scorer()]. If supplied,
#'   require the recorded gene panel and cell-type set. Matching these axes
#'   does not establish transportability to another tissue or experiment.
#' @return A data frame with one row per gene and cell type, seven predictors,
#'   `cohort_n`, and `theta_mean` (a diagnostic, not a GAM predictor).
#' @export
prepare_variation_features <- function(bulk, reference_mean, rna_fractions,
                                      scorer = NULL) {
  bulk <- .vscore_matrix(bulk, 'bulk')
  mu <- .vscore_matrix(reference_mean, 'reference_mean')
  theta <- .vscore_matrix(rna_fractions, 'rna_fractions')
  genes <- rownames(bulk); samples <- colnames(bulk); types <- colnames(mu)
  .vscore_assert(ncol(bulk) >= 3L, 'At least three samples are required.')
  .vscore_assert(setequal(rownames(mu), genes), 'Reference and bulk gene panels differ.')
  .vscore_assert(setequal(rownames(theta), samples), 'Fraction and bulk samples differ.')
  .vscore_assert(setequal(colnames(theta), types), 'Reference and fraction cell types differ.')
  mu <- mu[genes, types, drop = FALSE]; theta <- theta[samples, types, drop = FALSE]
  .vscore_assert(all(colSums(bulk) > 0), 'Zero bulk panel library; features undefined.')
  .vscore_assert(max(abs(colSums(mu) - 1)) < 1e-6, 'Reference columns must sum to one.')
  .vscore_assert(max(abs(rowSums(theta) - 1)) < 1e-6, 'RNA fraction rows must sum to one.')
  if (!is.null(scorer)) {
    check_variation_scorer(scorer, replay = FALSE)
    .vscore_assert(setequal(genes, scorer$gene_panel),
      'Modeled gene panel differs from the recorded benchmark. Revalidate the new context.')
    .vscore_assert(setequal(types, names(scorer$models)),
      'Cell-type set differs from the recorded benchmark.')
  }
  B <- sweep(bulk, 2L, colSums(bulk), '/') * 1e6
  logB <- log2(B + .5)
  X <- cbind(1, theta[, -ncol(theta), drop = FALSE]); qrX <- qr(X)
  residuals <- t(as.matrix(stats::lm.fit(X, t(logB))$residuals))
  tss <- rowSums((logB - rowMeans(logB))^2)
  H <- qr.Q(qrX)[, seq_len(qrX$rank), drop = FALSE]; h <- rowSums(H^2)
  loo <- if (all(h < 1 - 1e-8))
    1 - rowSums(sweep(residuals, 2L, 1 - h, '/')^2) / pmax(tss, 1e-12)
    else rep(NA_real_, length(genes))
  # Accumulate sample means without a genes x types x samples allocation array.
  contribution <- allocated <- matrix(0, nrow(mu), ncol(mu), dimnames = dimnames(mu))
  for (i in seq_along(samples)) {
    x <- sweep(mu, 2L, theta[i, ], '*')
    w <- x / pmax(rowSums(x), 1e-300)
    contribution <- contribution + w
    allocated <- allocated + w * bulk[, i]
  }
  contribution <- contribution / length(samples); allocated <- allocated / length(samples)
  bulk_sd <- apply(logB, 1L, stats::sd)
  max_share <- apply(bulk, 1L, max) / pmax(rowSums(bulk), 1)
  rows <- lapply(seq_along(types), function(k) data.frame(
    gene = genes, cell_type = types[k], cohort_n = length(samples),
    contribution_mean = contribution[, k],
    expected_log_mean_count = log2(allocated[, k] + .5),
    ref_specificity = mu[, k] / pmax(rowSums(mu), 1e-300),
    bulk_logCPM_sd = bulk_sd, bulk_max_donor_share = max_share,
    composition_loo_r2 = loo,
    bulk_theta_r = apply(B, 1L, function(x) .vscore_cor(x, theta[, k])),
    theta_mean = mean(theta[, k]), stringsAsFactors = FALSE, row.names = NULL))
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, 'modeled_gene_panel') <- genes
  attr(out, 'feature_schema') <- 'CDSeq2_GAM7_features_v1'
  out
}

#' Load a user-trained variation scorer
#'
#' @param file Explicit RDS path created by [save_variation_scorer()]. No default
#'   brain or other pretrained model is automatically loaded.
#' @param verify Replay the saved feature/prediction sentinel after loading.
#' @return The trained scorer object.
#' @export
load_variation_scorer <- function(file, verify = TRUE) {
  .vscore_assert(!missing(file) && length(file) == 1L && nzchar(file) && file.exists(file),
    'Supply the path to your trained scorer RDS file.')
  model <- readRDS(file)
  check_variation_scorer(model, replay = verify)
  model
}

#' Check the frozen variation scorer and saved prediction replay
#'
#' @param scorer A scorer object from [load_variation_scorer()].
#' @param replay Whether to replay its saved sentinel predictions.
#' @param tolerance Absolute numerical tolerance for replay.
#' @return Invisibly, TRUE; otherwise stops with a diagnostic.
#' @export
check_variation_scorer <- function(scorer, replay = TRUE, tolerance = 1e-8) {
  .vscore_assert(requireNamespace('mgcv', quietly = TRUE), 'Install mgcv to use the variation scorer.')
  spec <- .vscore_model_spec(scorer)
  .vscore_features <- spec$features
  .vscore_assert(is.list(scorer$models), 'Invalid cell-type models.')
  .vscore_ids(names(scorer$models), 'Model cell types')
  .vscore_ids(scorer$gene_panel, 'Model gene panel')
  for (ct in names(scorer$models)) {
    m <- scorer$models[[ct]]
    .vscore_assert(identical(m$kind, 'gam') && inherits(m$fit, 'gam') &&
      identical(names(m$lower), .vscore_features) && identical(names(m$upper), .vscore_features) &&
      all(is.finite(c(m$lower, m$upper))) && all(m$lower <= m$upper),
      'Invalid GAM or training bounds: ', ct)
  }
  if (replay) {
    .vscore_assert(is.list(scorer$sentinel) && is.data.frame(scorer$sentinel$features) &&
      is.data.frame(scorer$sentinel$predictions), 'Missing saved replay sentinel.')
    .vscore_replay(scorer$sentinel$features, scorer$sentinel$predictions, scorer, tolerance)
  }
  invisible(TRUE)
}

#' Score reliability of across-sample expression trends
#'
#' Predicts normalized-profile Pearson correlation on the Fisher-z scale with
#' the frozen cell-type GAM, then returns `predicted_r = tanh(z)` and
#' `variation_score = 100 * max(0, predicted_r)`.
#' The index is not a probability of a correct association, a p-value, or an
#' estimate of effect-size accuracy. It does not establish cell-type attribution.
#' Rank within cell type; inspect composition, support, and range flags.
#'
#' @param features Data frame returned by [prepare_variation_features()], or a
#'   saved table with the exact feature definitions required by the scorer and gene/type IDs.
#' @param scorer Explicit object from [train_variation_scorer()] or
#'   [load_variation_scorer()]. No fitted model is supplied automatically.
#' @param range_action Either `"flag"` (the policy used in the reported benchmark: retain finite
#'   predictions and flag extrapolation), or `"abstain"` (explicitly set flagged
#'   predictions to missing). Abstention changes the evaluated policy.
#' @return Data frame preserving row order, with identifiers, `predicted_r`,
#'   `variation_score`, `outside_training_range`, `range_features`, `status`,
#'   `panel_status`, and the two observable support-baseline columns. Missing
#'   predictors and unsupported types retain missing learned scores.
#' @export
score_variation <- function(features, scorer,
                            range_action = c('flag', 'abstain')) {
  range_action <- match.arg(range_action)
  check_variation_scorer(scorer, replay = FALSE)
  .vscore_predict_table(features, scorer, range_action)
}
.vscore_predict_table <- function(features, scorer, range_action = 'flag') {
  .vscore_features <- scorer$features
  .vscore_assert(identical(.vscore_features, .vscore_spec('GAM6')$features) ||
    identical(.vscore_features, .vscore_spec('GAM7')$features), 'Unsupported predictor list.')
  .vscore_assert(is.data.frame(features) && !anyDuplicated(names(features)) &&
    all(c('gene', 'cell_type', .vscore_features) %in% names(features)), 'Feature columns missing or duplicated.')
  .vscore_assert(all(vapply(features[.vscore_features], is.numeric, logical(1))),
    'Predictor columns must be numeric.')
  .vscore_assert(!anyNA(features$gene) && !anyNA(features$cell_type) &&
    all(nzchar(as.character(features$gene))) && all(nzchar(as.character(features$cell_type))) &&
    !anyDuplicated(.vscore_key(features)), 'Invalid or duplicate gene/type identifiers.')
  n <- nrow(features)
  out <- data.frame(gene = as.character(features$gene), cell_type = as.character(features$cell_type),
    predicted_r = rep(NA_real_, n), variation_score = rep(NA_real_, n),
    outside_training_range = rep(NA, n), range_features = rep('', n),
    status = rep('undefined_feature', n), stringsAsFactors = FALSE)
  # Domain checks apply only to defined values. Undefined values cause abstention below.
  for (nm in c('contribution_mean', 'ref_specificity', 'bulk_max_donor_share')) {
    v <- features[[nm]]
    .vscore_assert(all(!is.finite(v) | (v >= 0 & v <= 1 + 1e-8)), 'Invalid feature domain: ', nm)
  }
  v <- features$bulk_logCPM_sd
  .vscore_assert(all(!is.finite(v) | v >= 0), 'bulk_logCPM_sd cannot be negative.')
  v <- features$bulk_theta_r
  .vscore_assert(all(!is.finite(v) | abs(v) <= 1 + 1e-8), 'bulk_theta_r outside correlation range.')
  for (ct in unique(out$cell_type)) {
    ii <- which(out$cell_type == ct); m <- scorer$models[[ct]]
    if (is.null(m)) { out$status[ii] <- 'unsupported_cell_type'; next }
    ok <- rowSums(!is.finite(as.matrix(features[ii, .vscore_features, drop = FALSE]))) == 0L
    out$range_features[ii[!ok]] <- 'undefined_feature'
    jj <- ii[ok]; if (!length(jj)) next
    x <- .vscore_x(features[jj, , drop = FALSE], .vscore_features)
    p <- tanh(as.numeric(stats::predict(m$fit, newdata = as.data.frame(x), type = 'response')))
    .vscore_assert(all(is.finite(p)), 'GAM returned nonfinite predictions for ', ct)
    b <- sweep(x, 2L, m$lower - 1e-9, '<') | sweep(x, 2L, m$upper + 1e-9, '>')
    outside <- rowSums(b) > 0L
    out$predicted_r[jj] <- p; out$variation_score[jj] <- 100 * pmax(0, p)
    out$outside_training_range[jj] <- outside
    out$range_features[jj] <- apply(b, 1L, function(v) paste(colnames(x)[v], collapse = ';'))
    out$status[jj] <- ifelse(outside, 'outside_training_range', 'ok')
    if (range_action == 'abstain') {
      out$predicted_r[jj[outside]] <- NA_real_; out$variation_score[jj[outside]] <- NA_real_
      out$status[jj[outside]] <- 'abstained_outside_training_range'
    }
  }
  panel <- attr(features, 'modeled_gene_panel')
  panel_status <- if (is.null(panel)) 'not_recorded' else if
    (setequal(panel, scorer$gene_panel)) 'matches_recorded_panel' else 'different_panel'
  out$panel_status <- rep(panel_status, n)
  out$contribution_mean <- features$contribution_mean
  out$expected_log_mean_count <- features$expected_log_mean_count
  attr(out, 'model_id') <- scorer$model_id
  attr(out, 'context') <- scorer$context
  out
}
.vscore_replay <- function(features, saved, scorer, tolerance = 1e-8) {
  .vscore_assert(is.numeric(tolerance) && length(tolerance) == 1L && is.finite(tolerance) && tolerance >= 0,
    'Invalid replay tolerance.')
  now <- .vscore_predict_table(features, scorer)
  .vscore_assert(is.data.frame(saved) && !anyDuplicated(names(saved)) &&
    all(c('gene','cell_type') %in% names(saved)) && !anyNA(saved$gene) &&
    !anyNA(saved$cell_type), 'Invalid saved replay table.')
  .vscore_assert(!anyDuplicated(.vscore_key(saved)) &&
    setequal(.vscore_key(saved), .vscore_key(now)), 'Replay gene/type keys differ.')
  saved <- saved[match(.vscore_key(now), .vscore_key(saved)), , drop = FALSE]
  columns <- c(predicted_r = 'predicted_r', variation_score = 'variation_score',
    outside_training_range = 'outside_training_range', range_features = 'range_features')
  if (!all(columns %in% names(saved))) columns <- c(
    predicted_r = paste0(scorer$method, '_predicted_r'),
    variation_score = paste0(scorer$method, '_score'),
    outside_training_range = paste0(scorer$method, '_outside_range'),
    range_features = paste0(scorer$method, '_range_features'))
  .vscore_assert(all(unname(columns) %in% names(saved)), 'Saved replay columns missing.')
  max_error <- 0
  for (nm in names(columns)) {
    a <- now[[nm]]; b <- saved[[columns[[nm]]]]
    if (is.numeric(a)) {
      .vscore_assert(identical(is.na(a), is.na(b)) && all(is.finite(a[!is.na(a)])) &&
        all(is.finite(b[!is.na(b)])), 'Replay missingness differs: ', nm)
      e <- max(c(0, abs(a - b)), na.rm = TRUE)
      .vscore_assert(e <= tolerance, 'Saved GAM replay failed: ', nm, '; max error=', e)
      max_error <- max(max_error, e)
    } else .vscore_assert(identical(a, b), 'Saved replay flags/reasons differ: ', nm)
  }
  invisible(max_error)
}

#' Select a descriptive fraction of profiles within each cell type
#'
#' Uses the chosen finite ranking, descending, with lexical gene tie breaks.
#' For comparisons, supply the same rows with all compared rankings finite.
#' This is a ranking utility, not a validated clinical or FDR threshold.
#' It does not modify association p-values or their multiple-testing family.
#'
#' @param scores Output from [score_variation()].
#' @param fraction Fraction to retain in each type, between zero and one.
#' @param ranking Ranking column: variation_score, contribution_mean, or
#'   expected_log_mean_count. Baseline support columns are returned by
#'   [score_variation()]. Higher values rank first for all three.
#' @return Input table with `rank_within_type`, `n_eligible_in_type`, and
#'   `retained` columns. Undefined scores are not eligible.
#' @export
select_variation_profiles <- function(scores, fraction = .2,
    ranking = c('variation_score', 'contribution_mean', 'expected_log_mean_count')) {
  ranking <- match.arg(ranking)
  .vscore_assert(ranking %in% names(scores) && is.numeric(scores[[ranking]]),
    'Requested ranking column is missing or not numeric: ', ranking)
  .vscore_assert(length(fraction) == 1L && is.finite(fraction) &&
    fraction >= 0 && fraction <= 1, 'fraction must be between zero and one.')
  .vscore_assert(is.data.frame(scores) && all(c('gene', 'cell_type', 'variation_score') %in% names(scores)) &&
    !anyNA(scores$gene) && !anyNA(scores$cell_type) && !anyDuplicated(.vscore_key(scores)), 'Invalid scores.')
  out <- scores; out$rank_within_type <- rep(NA_integer_, nrow(out))
  out$n_eligible_in_type <- rep(0L, nrow(out)); out$retained <- rep(FALSE, nrow(out))
  for (ct in unique(out$cell_type)) {
    all_i <- which(out$cell_type == ct)
    ii <- all_i[is.finite(out[[ranking]][all_i])]
    ii <- ii[order(-out[[ranking]][ii], out$gene[ii], method = 'radix')]
    out$rank_within_type[ii] <- seq_along(ii)
    out$n_eligible_in_type[all_i] <- length(ii)
    if (length(ii) && fraction > 0) out$retained[head(ii, floor(fraction * length(ii) + 1e-9))] <- TRUE
  }
  out
}

.vscore_hash <- function(x) {
  f <- tempfile('variation_hash_'); on.exit(unlink(f), add = TRUE)
  writeBin(serialize(x, NULL, version = 2), f)
  unname(tools::md5sum(f))
}
.vscore_with_seed <- function(seed, fun) {
  .vscore_assert(length(seed) == 1L && is.finite(seed) && seed >= 0 &&
    seed <= .Machine$integer.max && seed == floor(seed), 'Invalid seed.')
  had <- exists('.Random.seed', .GlobalEnv, inherits = FALSE)
  if (had) old <- get('.Random.seed', .GlobalEnv)
  kinds <- RNGkind()
  on.exit({do.call(RNGkind, as.list(kinds))
    if (had) assign('.Random.seed', old, .GlobalEnv)
    else if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) rm('.Random.seed', envir=.GlobalEnv)})
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed);fun()
}
.vscore_array <- function(x, name) {
  .vscore_assert(is.array(x) && is.numeric(x) && length(dim(x)) == 3L &&
    all(dim(x) > 0L) && all(is.finite(x)) && all(x >= 0),
    name, ' must be a finite nonnegative genes x types x samples array.')
  for (j in seq_len(3L)) .vscore_ids(dimnames(x)[[j]], paste(name,'axis',j))
  x
}
.vscore_donors <- function(donor_ids, samples, reference_donor_ids) {
  .vscore_assert(is.character(donor_ids) && length(donor_ids) == length(samples) &&
    !is.null(names(donor_ids)) && !anyDuplicated(names(donor_ids)) &&
    setequal(names(donor_ids), samples) && !anyNA(donor_ids) && all(nzchar(donor_ids)),
    'donor_ids must map each named sample to its original biological donor ID.')
  .vscore_ids(reference_donor_ids, 'reference_donor_ids')
  .vscore_assert(!any(donor_ids %in% reference_donor_ids),
    'Benchmark target donors overlap reference donors. Split biological donors before simulation.')
  donor_ids[samples]
}
.vscore_profile_cor <- function(a, b, method = 'pearson') {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L || stats::sd(a[ok]) < 1e-12 || stats::sd(b[ok]) < 1e-12) return(NA_real_)
  unname(stats::cor(a[ok], b[ok], method = method))
}

#' Prepare training or validation targets from a known-component benchmark
#'
#' Features come only from observable inputs; responses come from correlations
#' between estimated and known within-type normalized profiles. The known
#' components must come from a simulation or independent measurement, never
#' from treating the CDSeq2 estimate itself as truth.
#'
#' @param bulk,reference_mean,rna_fractions Inputs to [prepare_variation_features()].
#' @param estimated_counts Posterior mean allocated reads, a named genes-by-types-
#'   samples array, summing over types to bulk counts (relative tolerance 1e-6).
#' @param truth_counts Known component array on the same axes. May be expected
#'   or realized counts, or directly aggregated observed counts. Declare the
#'   oracle explicitly and use the same definition across training/validation.
#' @param donor_ids Named character vector mapping samples to biological donor
#'   IDs. Repeated simulations must retain original donor IDs for leakage checks.
#' @param reference_donor_ids Unique biological donor IDs used for the reference.
#' @param cohort_id Unique identifier for this benchmark realization.
#' @param oracle Description of the target. For direct aggregation specify
#'   `"direct_observed_component_counts"`; the historical default is
#'   `"expected_components"`. No simulation is required by this function.
#' @param min_donors Minimum distinct biological donors per cohort (default 10).
#'   This is a software safeguard, not a validated minimum sample size.
#' @return A benchmark object accepted by [train_variation_scorer()] and
#'   [validate_variation_scorer()]. Zero-library profiles remain unavailable;
#'   their target correlations are not computed on a selected donor subset.
#' @export
prepare_variation_benchmark <- function(bulk, reference_mean, rna_fractions,
    estimated_counts, truth_counts, donor_ids, reference_donor_ids, cohort_id,
    oracle = 'expected_components', min_donors = 10L) {
  f <- prepare_variation_features(bulk, reference_mean, rna_fractions)
  bulk <- as.matrix(bulk); mu <- as.matrix(reference_mean)
  genes <- rownames(bulk); types <- colnames(mu); samples <- colnames(bulk)
  donor_ids <- .vscore_donors(donor_ids, samples, reference_donor_ids)
  .vscore_assert(length(min_donors)==1L && is.finite(min_donors) && min_donors>=3 &&
    min_donors==floor(min_donors), 'min_donors must be an integer of at least three.')
  .vscore_assert(length(unique(donor_ids)) >= min_donors, 'Too few distinct benchmark donors.')
  .vscore_assert(is.character(cohort_id) && length(cohort_id)==1L && !is.na(cohort_id) && nzchar(cohort_id), 'Invalid cohort_id.')
  .vscore_assert(is.character(oracle) && length(oracle)==1L && !is.na(oracle) && nzchar(oracle), 'Specify the oracle.')
  aa <- .vscore_array(estimated_counts,'estimated_counts'); tt <- .vscore_array(truth_counts,'truth_counts')
  for (a in list(aa,tt)) for (j in seq_len(3L)) .vscore_assert(
    setequal(dimnames(a)[[j]], list(genes,types,samples)[[j]]), 'Component array axes differ from bulk/reference.')
  aa <- aa[genes,types,samples,drop=FALSE]; tt <- tt[genes,types,samples,drop=FALSE]
  .vscore_assert(max(abs(apply(aa,c(1,3),sum)-bulk)/pmax(bulk,1)) < 1e-6,
    'estimated_counts do not conserve bulk counts. Supply allocated reads, not normalized expression.')
  r <- sp <- rep(NA_real_, nrow(f)); diagnostics <- list()
  for (k in seq_along(types)) {
    a <- matrix(aa[,k,], nrow=length(genes)); t <- matrix(tt[,k,],nrow=length(genes))
    la <- colSums(a); lt <- colSums(t); valid_library <- all(la>0 & lt>0)
    ii <- which(f$cell_type==types[k])
    if (valid_library) {
      a <- sweep(a,2L,la,'/');t <- sweep(t,2L,lt,'/')
      r[ii] <- vapply(seq_along(genes),function(g) .vscore_profile_cor(a[g,],t[g,]),numeric(1))
      sp[ii] <- vapply(seq_along(genes),function(g) .vscore_profile_cor(a[g,],t[g,],'spearman'),numeric(1))
    }
    diagnostics[[k]] <- data.frame(cell_type=types[k],samples=length(samples),
      distinct_donors=length(unique(donor_ids)),zero_estimated_libraries=sum(la==0),
      zero_truth_libraries=sum(lt==0),defined_targets=sum(is.finite(r[ii])))
  }
  f$cohort_id <- cohort_id;f$observed_r <- r;f$observed_spearman <- sp
  out <- list(data=f, gene_panel=genes,cell_types=types,donor_ids=unique(unname(donor_ids)),
    reference_donor_ids=reference_donor_ids,cohort_id=cohort_id,oracle=oracle,
    reference_hash=.vscore_hash(mu[genes,types,drop=FALSE]),diagnostics=do.call(rbind,diagnostics))
  out$fingerprint <- .vscore_hash(out)
  class(out) <- c('cdseq2_variation_benchmark','list');out
}

.vscore_benchmarks <- function(benchmarks) {
  if (inherits(benchmarks,'cdseq2_variation_benchmark')) benchmarks <- list(benchmarks)
  .vscore_assert(is.list(benchmarks) && length(benchmarks)>0L && all(vapply(benchmarks,
    inherits, logical(1), what='cdseq2_variation_benchmark')), 'Supply prepared variation benchmark objects.')
  ids <- vapply(benchmarks,`[[`,character(1),'cohort_id')
  .vscore_assert(!anyDuplicated(ids),'Benchmark cohort_id values must be unique.')
  base <- benchmarks[[1L]]
  for (b in benchmarks) {
    .vscore_assert(setequal(b$gene_panel,base$gene_panel) && setequal(b$cell_types,base$cell_types),
      'Training/validation benchmarks must use the same complete panel and cell-type set.')
    .vscore_assert(identical(b$oracle,base$oracle),'Do not mix oracle definitions.')
    z <- unclass(b);z$fingerprint <- NULL
    .vscore_assert(identical(.vscore_hash(z),b$fingerprint),'Benchmark was modified after preparation.')
  }
  targets <- unique(unlist(lapply(benchmarks,`[[`,'donor_ids'),use.names=FALSE))
  refs <- unique(unlist(lapply(benchmarks,`[[`,'reference_donor_ids'),use.names=FALSE))
  .vscore_assert(!any(targets %in% refs),'A target donor in one benchmark is a reference donor in another.')
  benchmarks
}
.vscore_fit <- function(d, min_profiles, features = .vscore_features) {
  active_features <- features
  ok <- is.finite(d$observed_r) & rowSums(!is.finite(as.matrix(d[active_features])))==0L
  dat <- d[ok,,drop=FALSE]
  .vscore_assert(nrow(dat)>=min_profiles && length(unique(dat$gene))>=min_profiles,
    'Insufficient distinct genes with defined targets/features in a cell-type fit. ',
    'Add benchmark support or revise the design; no fallback model was fitted.')
  x <- .vscore_x(dat, active_features); reg <- as.data.frame(x)
  reg$.y <- atanh(pmax(-.995,pmin(.995,dat$observed_r)))
  nu <- vapply(reg[colnames(x)],function(z)length(unique(z)),integer(1)); active <- names(nu)[nu>1L]
  terms <- vapply(active,function(nm) if(nu[[nm]]<5L) nm else paste0('s(',nm,', k=5, bs="cs")'),character(1))
  s <- mgcv::s
  formula <- stats::as.formula(paste('.y ~',if(length(terms))paste(terms,collapse=' + ')else'1'),env=environment())
  fit <- mgcv::gam(formula,data=reg,family=stats::gaussian(),method='REML',
    control=mgcv::gam.control(nthreads=1L))
  .vscore_assert(isTRUE(fit$converged),'GAM did not converge.')
  if (!is.null(fit$outer.info$conv)) .vscore_assert(identical(fit$outer.info$conv,'full convergence'),
    'GAM smoothing optimization did not fully converge.')
  for (nm in c('formula','terms','pterms','pred.formula')) if(!is.null(fit[[nm]]))
    environment(fit[[nm]]) <- asNamespace('mgcv')
  list(kind='gam',fit=fit,terms=terms,dropped=names(nu)[nu<=1L],
    lower=apply(x,2L,min),upper=apply(x,2L,max),n_training=nrow(dat),n_distinct_genes=length(unique(dat$gene)))
}
.vscore_evaluate <- function(d,retentions,model_name='GAM7') {
  .vscore_assert(is.numeric(retentions) && length(retentions)>0L && all(is.finite(retentions)) &&
    all(retentions>=0 & retentions<=1),'Invalid retention fractions.')
  summaries <- curves <- list(); ns <- nc <- 0L
  keys <- unique(d[c('cohort_id','cell_type')])
  for (u in seq_len(nrow(keys))) {
    z <- d[d$cohort_id==keys$cohort_id[u] & d$cell_type==keys$cell_type[u],,drop=FALSE]
    both <- is.finite(z$observed_r)&is.finite(z$predicted_r)
    ns<-ns+1L;summaries[[ns]]<-data.frame(cohort_id=keys$cohort_id[u],cell_type=keys$cell_type[u],
      profiles=nrow(z),defined_truth=sum(is.finite(z$observed_r)),finite_score=sum(is.finite(z$predicted_r)),
      evaluated=sum(both),MAE=if(any(both))mean(abs(z$predicted_r[both]-z$observed_r[both]))else NA_real_,
      RMSE=if(any(both))sqrt(mean((z$predicted_r[both]-z$observed_r[both])^2))else NA_real_,
      rank_Spearman=.vscore_profile_cor(z$predicted_r,z$observed_r,'spearman'),
      range_flagged=sum(z$outside_training_range %in% TRUE))
    # Common rankability depends on observed predictors, never on truth availability.
    common <- is.finite(z$variation_score)&is.finite(z$contribution_mean)&is.finite(z$expected_log_mean_count)
    rankers <- stats::setNames(c('variation_score','contribution_mean','expected_log_mean_count'),
      c(model_name,'contribution','expected_count'))
    for (method in names(rankers)) {
      ix <- which(common);ix <- ix[order(-z[[rankers[[method]]]][ix],z$gene[ix],method='radix')]
      for (fraction in sort(unique(retentions))) {
        selected <- head(ix,floor(length(ix)*fraction+1e-9)); r<-z$observed_r[selected]; sp<-z$observed_spearman[selected]
        good<-is.finite(r); gs<-is.finite(sp);nc<-nc+1L
        curves[[nc]]<-data.frame(cohort_id=keys$cohort_id[u],cell_type=keys$cell_type[u],method=method,
          retention=fraction,all_profiles=nrow(z),rankable=sum(common),retained=length(selected),
          defined_Pearson=sum(good),mean_Pearson=if(any(good))mean(r[good])else NA_real_,
          median_Pearson=if(any(good))stats::median(r[good])else NA_real_,
          fraction_Pearson_ge_0_5=if(any(good))mean(r[good]>=.5)else NA_real_,
          defined_Spearman=sum(gs),median_Spearman=if(any(gs))stats::median(sp[gs])else NA_real_)
      }
    }
  }
  list(by_cohort_and_type=do.call(rbind,summaries),retention=do.call(rbind,curves))
}

#' Train an application-specific variation scorer from known-component benchmarks
#'
#' Fits separate GAM6 (default) or historical GAM7 models by cell type. The response
#' is Fisher-transformed normalized Pearson recovery, not association correctness.
#' Grouped gene cross-validation is a development diagnostic; it does not hold
#' out donors. Use [validate_variation_scorer()] with untouched biological donors.
#'
#' @param benchmarks A prepared benchmark or a list of such objects.
#' @param context Named list describing `application`, `tissue`,
#'   `deconvolution_method`, and `fit_settings`. Include disease/condition and
#'   simulation design as additional entries when relevant.
#' @param gene_groups Optional named character vector mapping every gene to
#'   its linked-gene group. Default groups by gene identity. Supply groups for
#'   linked injections, correlated construction units, or other known dependencies.
#' @param folds Number of grouped development folds, at least two (default five).
#' @param seed Seed for group-to-fold allocation; caller RNG state is restored.
#' @param min_profiles Minimum distinct eligible genes in each type/fold fit,
#'   default 100, at least 50. This is a safeguard, not a sample-size guarantee.
#' @param retentions Fractions for per-cohort/type ranking diagnostics.
#' @param model GAM6 excludes composition_loo_r2; GAM7 reproduces the historical
#'   seven-feature specification. New fits default to GAM6.
#' @return A user-trained scorer, including grouped OOF predictions, diagnostic
#'   metrics, donor provenance and replay sentinel. It is not automatically
#'   declared independently validated, and no held-out validation data are refitted.
#' @export
train_variation_scorer <- function(benchmarks, context, gene_groups=NULL, folds=5L,
    seed=1L,min_profiles=100L,retentions=seq(0,1,.05), model=c('GAM6','GAM7')) {
  spec <- .vscore_spec(match.arg(model)); .vscore_features <- spec$features
  .vscore_assert(requireNamespace('mgcv',quietly=TRUE),'Install mgcv.')
  bs <- .vscore_benchmarks(benchmarks);base<-bs[[1L]]
  required<-c('application','tissue','deconvolution_method','fit_settings')
  .vscore_assert(is.list(context)&&all(required%in%names(context))&&
    all(vapply(context[required],function(x)length(x)>0L,logical(1))),
    'Provide context with application, tissue, deconvolution_method, and fit_settings.')
  .vscore_assert(length(folds)==1L&&is.finite(folds)&&folds>=2L&&folds==floor(folds),'folds must be an integer >=2.')
  .vscore_assert(length(min_profiles)==1L&&is.finite(min_profiles)&&min_profiles>=50L&&min_profiles==floor(min_profiles),
    'min_profiles must be an integer >=50.')
  genes<-base$gene_panel
  if(is.null(gene_groups))gene_groups<-stats::setNames(genes,genes)
  .vscore_assert(is.character(gene_groups)&&!anyNA(gene_groups)&&all(nzchar(gene_groups))&&
    !is.null(names(gene_groups))&&!anyDuplicated(names(gene_groups))&&setequal(names(gene_groups),genes),
    'gene_groups must name every gene exactly once.')
  group_names<-sort(unique(gene_groups),method='radix')
  .vscore_assert(length(group_names)>=folds,'Too few gene groups for requested folds.')
  gf<-.vscore_with_seed(seed,function()stats::setNames(sample(rep(seq_len(folds),length.out=length(group_names))),group_names))
  d<-do.call(rbind,lapply(bs,`[[`,'data'));rownames(d)<-NULL
  d$gene_group<-unname(gene_groups[d$gene]);d$fold<-unname(gf[d$gene_group])
  types<-base$cell_types; models<-list();oof<-list();nr<-0L
  for(ct in types) {
    dd<-d[d$cell_type==ct,,drop=FALSE]
    message('Training variation GAM | ',ct,' | ',folds,' grouped folds + full fit')
    for(fold in seq_len(folds)) {
      m<-.vscore_fit(dd[dd$fold!=fold,,drop=FALSE],min_profiles,.vscore_features)
      test<-dd[dd$fold==fold,,drop=FALSE]
      for(cid in unique(test$cohort_id)) {
        test_c<-test[test$cohort_id==cid,,drop=FALSE]
        pp<-.vscore_predict_table(test_c,list(features=.vscore_features,models=stats::setNames(list(m),ct),gene_panel=genes,model_id='OOF',context=context))
        nr<-nr+1L;oof[[nr]]<-.vscore_join(test_c,pp)
      }
    }
    models[[ct]]<-.vscore_fit(dd,min_profiles,.vscore_features)
  }
  oof<-do.call(rbind,oof);rownames(oof)<-NULL
  scorer<-list(schema=spec$schema,method=spec$method,features=.vscore_features,log_features=spec$logs,
    models=models,gene_panel=genes,context=context,fit_settings=context$fit_settings,
    target='Across-donor normalized Pearson correlation',oracle=base$oracle,
    donor_ids=unique(unlist(lapply(bs,`[[`,'donor_ids'),use.names=FALSE)),
    reference_donor_ids=unique(unlist(lapply(bs,`[[`,'reference_donor_ids'),use.names=FALSE)),
    reference_hashes=unique(vapply(bs,`[[`,character(1),'reference_hash')),
    training_benchmarks=vapply(bs,`[[`,character(1),'fingerprint'),
    validation_status='grouped_gene_development_CV_only',
    development=list(predictions=oof,metrics=.vscore_evaluate(oof,retentions,spec$name),
      fold_by_gene=stats::setNames(unname(gf[gene_groups[genes]]),genes)),
    provenance=list(interface_version='2.4.0',created=as.character(Sys.time()),
      R=R.version.string,mgcv=as.character(utils::packageVersion('mgcv')),seed=seed,folds=folds))
  scorer$model_id<-paste0(spec$name,'_',.vscore_hash(list(models=models,training=scorer$training_benchmarks,context=context)))
  sentinel<-bs[[1L]]$data
  sentinel<-sentinel[unlist(lapply(types,function(ct)head(which(sentinel$cell_type==ct),10L))),,drop=FALSE]
  scorer$sentinel<-.vscore_sentinel(sentinel,scorer)
  class(scorer)<-c('cdseq2_variation_scorer','list');check_variation_scorer(scorer);scorer
}

#' Evaluate a trained variation scorer on untouched donor benchmarks
#'
#' Verifies that biological target donors were absent from both training and
#' reference sets. Sharing the reference donors between training and validation
#' is allowed; validation target donors cannot contribute to either reference.
#' No model fitting or threshold optimization is performed. Once examined,
#' validation data must not be reused as untouched evidence after model tuning.
#'
#' @param scorer User-trained model from [train_variation_scorer()].
#' @param benchmarks New donor benchmark object(s).
#' @param retentions Fractions used for descriptive retention curves.
#' @return Predictions, per-cohort/type metrics, ranking curves and donor
#'   separation checks. This evaluates profile recovery, not downstream
#'   association power, FDR or cell-type attribution.
#' @export
validate_variation_scorer <- function(scorer,benchmarks,retentions=seq(0,1,.05)) {
  check_variation_scorer(scorer)
  .vscore_assert(length(scorer$donor_ids)>0L&&length(scorer$reference_donor_ids)>0L,
    'Training donor provenance missing; cannot certify donor separation.')
  bs<-.vscore_benchmarks(benchmarks);rows<-checks<-list()
  for(i in seq_along(bs)) {
    b<-bs[[i]]
    .vscore_assert(setequal(b$gene_panel,scorer$gene_panel)&&setequal(b$cell_types,names(scorer$models)),
      'Validation panel or cell types differ from training.')
    .vscore_assert(identical(b$oracle,scorer$oracle),'Validation oracle differs from training.')
    .vscore_assert(!any(b$donor_ids%in%c(scorer$donor_ids,scorer$reference_donor_ids))&&
      !any(b$reference_donor_ids%in%scorer$donor_ids),
      'Validation and training/reference donors overlap. Use new biological donors, not new sample labels.')
    p<-score_variation(b$data,scorer)
    rows[[i]]<-.vscore_join(b$data,p)
    checks[[i]]<-data.frame(cohort_id=b$cohort_id,n_donors=length(b$donor_ids),
      training_target_overlap=0L,training_reference_overlap=0L,
      reference_previously_seen=b$reference_hash%in%scorer$reference_hashes)
  }
  d<-do.call(rbind,rows)
  list(model_id=scorer$model_id,status='completed_independent_donor_profile_evaluation',
    predictions=d,metrics=.vscore_evaluate(d,retentions,.vscore_model_spec(scorer)$name),donor_checks=do.call(rbind,checks),
    scope='Profile fidelity; downstream association must be evaluated separately. No automatic acceptance decision.')
}

#' Save an explicitly trained variation scorer
#' @param scorer Model returned by [train_variation_scorer()].
#' @param file Destination RDS path.
#' @return Invisibly, the saved path. Existing files are never overwritten.
#' @export
save_variation_scorer <- function(scorer,file) {
  check_variation_scorer(scorer)
  .vscore_assert(length(file)==1L&&!is.na(file)&&nzchar(file)&&!file.exists(file),'Choose a new model file path.')
  dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile('variation_scorer_',tmpdir=dirname(file));on.exit(unlink(tmp),add=TRUE)
  saveRDS(scorer,tmp,version=2);load_variation_scorer(tmp)
  .vscore_assert(file.rename(tmp,file),'Cannot save scorer.');invisible(normalizePath(file))
}

#' Simulate bulk mixtures from donor-resolved component counts
#'
#' Preserves the observed within-type donor expression profiles and simulates
#' multinomial reads under chosen library sizes and RNA fractions. This simple
#' matched-context generator does not simulate every biological or platform
#' mismatch. Use distinct biological donor sets for reference, development and
#' validation. Shared reference means alone are insufficient input: they do not
#' supply known donor-to-donor expression variation.
#'
#' @param component_counts Named genes-by-types-by-samples array of observed
#'   donor-resolved single-cell/nucleus component counts. Every type/sample
#'   library must be positive; absent types are not imputed.
#' @param donor_ids,reference_donor_ids Biological donor identifiers as in
#'   [prepare_variation_benchmark()].
#' @param library_sizes Optional named sample library totals (positive integers).
#'   Default uses the total component reads for each sample.
#' @param rna_fractions Optional named samples-by-types RNA fractions. Default
#'   uses the observed component library proportions.
#' @param seed Simulation seed; caller RNG state is restored.
#' @return A list containing bulk counts, expected `truth_counts`,
#'   `realized_counts`, donor provenance and simulation settings. Run CDSeq2 on
#'   the simulated bulk separately, then call [prepare_variation_benchmark()].
#' @export
simulate_variation_mixtures <- function(component_counts, donor_ids,
    reference_donor_ids, library_sizes=NULL, rna_fractions=NULL, seed=1L) {
  a<-.vscore_array(component_counts,'component_counts')
  genes<-dimnames(a)[[1]];types<-dimnames(a)[[2]];samples<-dimnames(a)[[3]]
  donor_ids<-.vscore_donors(donor_ids,samples,reference_donor_ids)
  lib<-apply(a,c(2,3),sum)
  .vscore_assert(all(lib>0),'Every source type/donor must have a positive component library; absent types need a separate design.')
  if(is.null(library_sizes))library_sizes<-colSums(lib)
  .vscore_assert(is.numeric(library_sizes)&&!anyDuplicated(names(library_sizes))&&
    setequal(names(library_sizes),samples)&&all(is.finite(library_sizes))&&
    all(library_sizes>0 & library_sizes<=.Machine$integer.max & library_sizes==floor(library_sizes)),
    'library_sizes must be named positive integer totals within R multinomial limits.')
  library_sizes<-library_sizes[samples]
  if(is.null(rna_fractions))rna_fractions<-t(sweep(lib,2,colSums(lib),'/'))
  theta<-.vscore_matrix(rna_fractions,'rna_fractions')
  .vscore_assert(setequal(rownames(theta),samples)&&setequal(colnames(theta),types),
    'Simulation RNA fraction axes differ.')
  theta<-theta[samples,types,drop=FALSE]
  .vscore_assert(max(abs(rowSums(theta)-1))<1e-6,'Simulation RNA fraction rows must sum to one.')
  q<-sweep(a,c(2,3),lib,'/')
  expected<-sweep(q,c(2,3),t(theta*library_sizes),'*')
  realized<-.vscore_with_seed(seed,function(){
    b<-array(0,dim(a),dimnames(a))
    for(i in seq_along(samples)) {
      l<-as.numeric(stats::rmultinom(1L,library_sizes[i],theta[i,]))
      for(k in seq_along(types))if(l[k]>0)b[,k,i]<-as.numeric(stats::rmultinom(1L,l[k],q[,k,i]))
    };b
  })
  bulk<-apply(realized,c(1,3),sum)
  .vscore_assert(identical(as.numeric(colSums(bulk)),as.numeric(library_sizes)),'Simulated count conservation failed.')
  list(bulk=bulk,truth_counts=expected,realized_counts=realized,
    donor_ids=donor_ids,reference_donor_ids=reference_donor_ids,
    settings=list(generator='donor_component_multinomial_v1',seed=seed,
      library_sizes=library_sizes,rna_fractions=theta,
      oracle='expected_components',component_hash=.vscore_hash(a)))
}

#' Train variation GAMs from saved observable features and recovery targets
#'
#' This scoring-only entry point preserves saved gene-group folds and eligibility
#' rules. It uses the same numerical fitter as [train_variation_scorer()].
#' No deconvolution is run. Gene-group CV does not establish new-donor validity.
#'
#' @param data Data frame containing gene, cell_type, cohort_id, observed_r and
#'   the selected model's columns from [prepare_variation_features()]. Optional logical
#'   training_eligible controls training inclusion; omitted means all rows.
#'   Optional observed_spearman is an evaluation diagnostic only. Supply known
#'   component recovery correlations; never use estimates as their own truth.
#' @param context Named list containing application, tissue,
#'   deconvolution_method and fit_settings, plus relevant simulation context.
#' @param gene_groups Named vector mapping every gene to a dependency group.
#' @param fold_by_gene Named integer vector giving every gene its saved fold.
#'   All occurrences of a gene/group must share a fold. Labels may start at zero.
#' @param donor_ids,reference_donor_ids Optional biological donor provenance.
#'   Supply both, or neither. Without these, the result cannot pass the donor
#'   separation checks in [validate_variation_scorer()].
#' @param oracle Explicit description of the benchmark recovery target.
#' @param workers Number of fresh PSOCK workers; default one. Each GAM uses one
#'   mgcv thread. Worker libraries are matched to the current R library paths.
#' @param cache_dir Optional directory for individually verified model caches.
#'   An active/stale lock stops reuse; it is never removed automatically.
#' @param min_profiles Minimum distinct eligible genes in each fit; default 100.
#' @param retentions Fractions used for profile diagnostics.
#' @param model GAM6 (default) or historical GAM7. For reproducing the ablation,
#'   supply its original training_eligible mask and saved folds; the ablation
#'   retained seven-feature eligibility even though GAM6 uses only six predictors.
#' @return An explicit user-trained model with OOF predictions, fit audit,
#'   input fingerprint, saved folds, context and prediction replay examples.
#' @export
train_variation_scorer_from_features <- function(data, context, gene_groups,
    fold_by_gene, donor_ids=NULL, reference_donor_ids=NULL,
    oracle='observed_component_normalized_Pearson', workers=1L,
    cache_dir=NULL, min_profiles=100L, retentions=seq(0,1,.05), model=c('GAM6','GAM7')) {
  spec <- .vscore_spec(match.arg(model)); .vscore_features <- spec$features
  .vscore_assert(requireNamespace('mgcv',quietly=TRUE),'Install mgcv.')
  required<-c('gene','cell_type','cohort_id','observed_r',.vscore_features)
  .vscore_assert(is.data.frame(data)&&nrow(data)>0L&&all(required%in%names(data)),
    'Missing training identifiers, features or observed_r.')
  d<-data;rownames(d)<-NULL
  for(nm in c('gene','cell_type','cohort_id')) {
    d[[nm]]<-as.character(d[[nm]])
    .vscore_assert(!anyNA(d[[nm]])&&all(nzchar(d[[nm]])),'Invalid training identifiers.')
  }
  .vscore_assert(!anyDuplicated(paste(d$cohort_id,d$cell_type,d$gene,sep='\r')),
    'Duplicated cohort/gene/type rows.')
  .vscore_assert(is.numeric(d$observed_r)&&all(is.na(d$observed_r)|
    (is.finite(d$observed_r)&abs(d$observed_r)<=1+1e-8)), 'Invalid observed_r target.')
  .vscore_assert(all(vapply(d[.vscore_features],is.numeric,logical(1))), 'Features must be numeric.')
  for(nm in .vscore_log_features) .vscore_assert(all(!is.finite(d[[nm]])|d[[nm]]>=0),
    'Negative nonnegative predictor: ',nm)
  if(!'training_eligible'%in%names(d))d$training_eligible<-TRUE
  .vscore_assert(is.logical(d$training_eligible)&&!anyNA(d$training_eligible),'training_eligible must be logical.')
  if(!'observed_spearman'%in%names(d))d$observed_spearman<-NA_real_
  genes<-unique(d$gene);types<-sort(unique(d$cell_type)); groups<-as.character(gene_groups)
  names(groups)<-names(gene_groups)
  for(v in list(groups,fold_by_gene)) .vscore_assert(!is.null(names(v))&&
    !anyDuplicated(names(v))&&setequal(names(v),genes)&&!anyNA(v),'Gene maps must name every gene once.')
  .vscore_assert(all(nzchar(groups))&&is.numeric(fold_by_gene)&&all(is.finite(fold_by_gene))&&
    all(fold_by_gene==floor(fold_by_gene))&&length(unique(fold_by_gene))>=2L,'Invalid group/fold map.')
  .vscore_assert(all(vapply(split(fold_by_gene[names(groups)],groups),function(x)
    length(unique(x))==1L,logical(1))),'A linked gene group crosses folds.')
  for(ix in split(seq_len(nrow(d)),paste(d$cohort_id,d$cell_type,sep='\r')))
    .vscore_assert(setequal(d$gene[ix],genes),'Each cohort/type must contain the complete modeled gene panel.')
  ctx<-c('application','tissue','deconvolution_method','fit_settings')
  .vscore_assert(is.list(context)&&all(ctx%in%names(context))&&
    all(vapply(context[ctx],function(x)length(x)>0L,logical(1))),'Incomplete application context.')
  .vscore_assert(length(oracle)==1L&&is.character(oracle)&&!is.na(oracle)&&nzchar(oracle),'Specify oracle.')
  if(length(donor_ids)||length(reference_donor_ids)) {
    .vscore_ids(donor_ids,'donor_ids');.vscore_ids(reference_donor_ids,'reference_donor_ids')
    .vscore_assert(!any(donor_ids%in%reference_donor_ids),'Target/reference donors overlap.')
  }
  .vscore_assert(length(workers)==1L&&is.finite(workers)&&workers>=1&&workers==floor(workers),'Invalid workers.')
  .vscore_assert(length(min_profiles)==1L&&is.finite(min_profiles)&&min_profiles>=50&&
    min_profiles==floor(min_profiles),'min_profiles must be an integer >=50.')
  d$gene_group<-unname(groups[d$gene]);d$fold<-unname(fold_by_gene[d$gene])
  runtime<-list(R=R.version.string,mgcv=as.character(utils::packageVersion('mgcv')),
    platform=R.version$platform,BLAS=unname(extSoftVersion()['BLAS']))
  specification<-list(features=.vscore_features,logs=.vscore_log_features,
    fitter=paste(deparse(.vscore_fit),collapse='\n'))
  fingerprint<-.vscore_hash(list(d,context,groups,fold_by_gene,oracle,runtime,specification,min_profiles,
    interface_version='2.4.0',method=spec$method,
    training_task=paste(deparse(.vscore_training_task),collapse='\n')))
  if(!is.null(cache_dir)) {
    dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
    lock<-file.path(cache_dir,'TRAINING.lock')
    .vscore_assert(dir.create(lock,showWarnings=FALSE),'Another training call is active or a stale lock remains: ',lock)
    writeLines(c(paste('PID',Sys.getpid()),paste('HOST',Sys.info()[['nodename']])),file.path(lock,'owner.txt'))
    on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  }
  folds<-sort(unique(d$fold));tasks<-unlist(lapply(seq_along(types),function(k)
    lapply(seq_len(length(folds)+1L),function(j) list(cell_type=types[k],
      fold=if(j<=length(folds))folds[j]else NA_real_,file=sprintf('type%02d_fit%02d.rds',k,j)))),recursive=FALSE)
  # Clone only pure numerical helpers into a base-parent environment. This avoids
  # serializing the CDSeq2 namespace or a load_all() environment into workers.
  module<-new.env(parent=baseenv())
  helpers<-c('.vscore_assert','.vscore_hash','.vscore_x','.vscore_fit','.vscore_training_task')
  for(nm in helpers){f<-get(nm,envir=environment());environment(f)<-module;assign(nm,f,module)}
  module$.vscore_features<-.vscore_features;module$.vscore_log_features<-.vscore_log_features
  module$d<-d;module$min_profiles<-min_profiles;module$fingerprint<-fingerprint;module$cache_dir<-cache_dir
  worker<-function(task) .vscore_training_task(task,d,min_profiles,fingerprint,cache_dir)
  environment(worker)<-module
  message(length(tasks),' scoring fits; no CDSeq2 deconvolution. Workers: ',min(workers,length(tasks)))
  if(workers>1L) {
    cl<-parallel::makePSOCKcluster(min(workers,length(tasks)),rscript_args='--vanilla')
    on.exit(parallel::stopCluster(cl),add=TRUE)
    info<-parallel::clusterCall(cl,function(paths){.libPaths(paths);list(R=R.version.string,
      mgcv=as.character(utils::packageVersion('mgcv')))},.libPaths())
    .vscore_assert(all(vapply(info,function(z)identical(z$R,runtime$R)&&identical(z$mgcv,runtime$mgcv),logical(1))),
      'Worker R/mgcv versions differ.')
    ans<-parallel::parLapplyLB(cl,tasks,worker)
  } else ans<-lapply(tasks,function(t){message(t$cell_type,' | ',if(is.na(t$fold))'full'else paste('fold',t$fold));worker(t)})
  models<-list();oof<-d;oof$predicted_r<-NA_real_;oof$outside_training_range<-NA
  audit<-list()
  for(i in seq_along(ans)) {
    z<-ans[[i]]
    if(is.na(z$fold))models[[z$cell_type]]<-z$model else {
      oof$predicted_r[z$rows]<-z$prediction;oof$outside_training_range[z$rows]<-z$outside
    }
    audit[[i]]<-data.frame(cell_type=z$cell_type,fold=if(is.na(z$fold))'full'else as.character(z$fold),
      training_rows=z$model$n_training,training_genes=z$model$n_distinct_genes,
      training_groups=z$training_groups,cached=z$cached)
  }
  oof$variation_score<-100*pmax(0,oof$predicted_r)
  model<-list(schema=spec$schema,method=spec$method,features=.vscore_features,log_features=spec$logs,
    models=models,gene_panel=genes,context=context,fit_settings=context$fit_settings,
    target='Across-donor normalized Pearson correlation',oracle=oracle,
    donor_ids=donor_ids,reference_donor_ids=reference_donor_ids,reference_hashes=character(),
    training_benchmarks=fingerprint,validation_status='grouped_gene_development_CV_only',
    donor_provenance=if(length(donor_ids))'provided'else'not_in_saved_feature_table',
    development=list(predictions=oof,metrics=.vscore_evaluate(oof,retentions,spec$name),
      fold_by_gene=fold_by_gene,fit_audit=do.call(rbind,audit)),
    provenance=list(interface_version='2.4.0',runtime=runtime,created=as.character(Sys.time()),
      input_fingerprint=fingerprint,specification=specification))
  model$model_id<-paste0(spec$name,'_',fingerprint)
  sent<-d[d$cohort_id==d$cohort_id[1],,drop=FALSE]
  sent<-sent[unlist(lapply(types,function(ct)head(which(sent$cell_type==ct),10L))),,drop=FALSE]
  model$sentinel<-.vscore_sentinel(sent,model)
  class(model)<-c('cdseq2_variation_scorer','list');check_variation_scorer(model);model
}
.vscore_training_task <- function(task,d,min_profiles,fingerprint,cache_dir) {
  path<-if(!is.null(cache_dir))file.path(cache_dir,task$file)else NULL
  if(!is.null(path)&&file.exists(path)) {
    .vscore_assert(file.exists(paste0(path,'.md5'))&&identical(unname(tools::md5sum(path)),
      readLines(paste0(path,'.md5'),warn=FALSE)),'Changed/incomplete scoring checkpoint: ',path)
    z<-readRDS(path)
    .vscore_assert(identical(z$fingerprint,fingerprint)&&identical(z$cell_type,task$cell_type)&&
      identical(z$fold,task$fold),'Checkpoint context/runtime differs. Use a new cache_dir.')
    z$cached<-TRUE;return(z)
  }
  train<-which(d$cell_type==task$cell_type&d$training_eligible&
    (is.na(task$fold)|d$fold!=task$fold))
  test<-if(is.na(task$fold))integer()else which(d$cell_type==task$cell_type&d$fold==task$fold)
  .vscore_assert(!length(intersect(d$gene_group[train],d$gene_group[test])),'Grouped CV leakage.')
  m<-.vscore_fit(d[train,,drop=FALSE],min_profiles)
  x<-.vscore_x(d[test,,drop=FALSE]);valid<-rowSums(!is.finite(x))==0L
  pred<-rep(NA_real_,length(test));outside<-rep(NA,length(test))
  if(any(valid)) {
    pred[valid]<-tanh(as.numeric(stats::predict(m$fit,newdata=as.data.frame(x[valid,,drop=FALSE]),type='response')))
    .vscore_assert(all(is.finite(pred[valid])), 'GAM returned nonfinite OOF predictions.')
    outside[valid]<-rowSums(sweep(x[valid,,drop=FALSE],2,m$lower-1e-9,'<')|
      sweep(x[valid,,drop=FALSE],2,m$upper+1e-9,'>'))>0L
  }
  z<-list(fingerprint=fingerprint,cell_type=task$cell_type,fold=task$fold,model=m,
    rows=test,prediction=pred,outside=outside,training_groups=length(unique(d$gene_group[train])),cached=FALSE)
  if(!is.null(path)) {
    tmp<-tempfile(tmpdir=dirname(path));on.exit(unlink(tmp),add=TRUE);saveRDS(z,tmp,version=2)
    .vscore_assert(file.rename(tmp,path),'Cannot write checkpoint.')
    writeLines(unname(tools::md5sum(path)),paste0(path,'.md5'))
  }
  z
}


#' Import the explicitly fitted GAM6 ablation without refitting
#'
#' Adapts the saved direct-aggregation GAM6 schema to the package interface.
#' Fitted GAMs, bounds, features, model ID and original provenance are preserved.
#' The original sentinel must replay before an object is returned. Optional full
#' replay tables check all saved predictions in addition to the original sentinel.
#' This import does not turn post-validation model development into independent
#' validation and does not fit or supply a pretrained model.
#'
#' @param scorer An explicit GAM6 ablation object or its RDS file path.
#' @param features,predictions Optional full saved feature and prediction tables,
#'   supplied together, each with unique gene/cell_type keys for one cohort.
#'   Prediction columns use the names returned by score_variation().
#' @param tolerance Maximum absolute replay error; default 1e-8.
#' @return A package-compatible scorer; the supplied object/file is not changed.
#' @export
import_variation_scorer <- function(scorer, features=NULL, predictions=NULL,
                                   tolerance=1e-8) {
  if (is.character(scorer)) {
    .vscore_assert(length(scorer)==1L && !is.na(scorer) && file.exists(scorer),
      'Supply an existing GAM6 ablation RDS file.')
    scorer <- readRDS(scorer)
  }
  .vscore_assert(is.list(scorer) &&
    identical(scorer$schema,'CDSeq2_GAM6_no_composition_loo_ablation_v1') &&
    identical(scorer$method,'GAM6_no_composition_loo_r2'),
    'Import supports only the saved direct GAM6 ablation. Use load_variation_scorer() for package models.')
  spec <- .vscore_spec('GAM6')
  .vscore_assert(identical(scorer$features,spec$features) &&
    identical(scorer$log_features,spec$logs), 'Ablation feature definitions differ.')
  .vscore_assert(is.list(scorer$context) && is.character(scorer$model_id) &&
    length(scorer$model_id)==1L && !is.na(scorer$model_id) && nzchar(scorer$model_id),
    'Ablation context/model ID missing.')
  .vscore_ids(scorer$donor_ids, 'Ablation training donors')
  .vscore_ids(scorer$reference_donor_ids, 'Ablation reference donors')
  .vscore_assert(!any(scorer$donor_ids %in% scorer$reference_donor_ids),
    'Ablation target/reference donors overlap.')
  .vscore_assert(identical(is.null(features),is.null(predictions)),
    'Supply both full replay features and predictions, or neither.')
  source_schema <- scorer$schema
  scorer$schema <- spec$schema
  scorer$validation_status <- 'exploratory_after_validation_review'
  scorer$fit_settings <- scorer$context$fit_settings
  if (is.null(scorer$reference_hashes)) scorer$reference_hashes <- character()
  scorer$import_provenance <- list(interface_version='2.4.0',
    source_schema=source_schema, source_model_id=scorer$model_id,
    imported=as.character(Sys.time()), refitted=FALSE,
    source_status=scorer$status, full_replay_rows=if(is.null(features))0L else nrow(features))
  class(scorer) <- unique(c('cdseq2_variation_scorer',class(scorer)))
  check_variation_scorer(scorer, replay=TRUE, tolerance=tolerance)
  if (!is.null(features)) .vscore_replay(features,predictions,scorer,tolerance)
  scorer
}
