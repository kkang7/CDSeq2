#' CDSeq4: donor-calibrated cell-type-specific reference precision
#'
#' `cdseq4()` is an empirical-Bayes extension of `CDSeq2::cdseq2()` for
#' individual-level complete deconvolution. It retains the CDSeq2 likelihood,
#' Gibbs sampler, Rao--Blackwellized estimator, and output structure, while
#' replacing one global expression-prior precision with a named precision for
#' each cell type.
#'
#' Reference cells are aggregated into donor-by-cell-type pseudobulks. For each
#' cell type, a Dirichlet--multinomial marginal likelihood estimates the
#' relative amount of donor-to-donor profile dispersion. Because absolute
#' dispersion estimates can be unstable in a compact reference, CDSeq4
#' conservatively calibrates the cell-type estimates around a prespecified
#' global precision. Thus, the global regularization scale is retained while
#' cell types with greater reference-donor variability receive relatively lower
#' precision and more stable cell types receive relatively higher precision.
#'
#' The calibration uses only the single-cell reference and its donor labels. The resulting named precision vector is passed through the native `reference_precision` interface of the revised `cdseq2()` implementation.
#' Bulk-sample truth, target condition labels, and held-out donor expression are
#' never used.
#'
#' @param bulk_data Gene-by-sample raw bulk count matrix.
#' @param mcmc_iterations Number of MCMC iterations passed to
#'   `CDSeq2::cdseq2()`.
#' @param cpu_number Number of workers passed to `CDSeq2::cdseq2()`.
#' @param reference_scRNA_count Gene-by-cell reference count matrix.
#' @param reference_scRNA_annotation Optional cell-type label for every
#'   reference cell. When `NULL`, labels are read from `reference_metadata`.
#' @param reference_metadata Optional data frame containing reference cell IDs,
#'   donors, and cell types. Rows are aligned to the columns of
#'   `reference_scRNA_count` by `reference_cell_id_col` or row names.
#' @param reference_scRNA_donor Optional donor label for every reference cell.
#'   When both this vector and `reference_metadata` are supplied, they must
#'   agree after cell-ID alignment.
#' @param reference_cell_id_col Column in `reference_metadata` containing cell
#'   IDs.
#' @param reference_donor_col Column in `reference_metadata` containing donor
#'   IDs.
#' @param reference_cell_type_col Column in `reference_metadata` containing
#'   cell-type labels.
#' @param reference_aggregation Reference profile aggregation passed to
#'   `CDSeq2::cdseq2()`. The default is `"pooled_counts"` so the first CDSeq4
#'   trial differs from the CDSeq2 baseline only through prior precision.
#' @param reference_min_cells_per_donor_type Minimum donor-by-cell-type cell
#'   count used by the underlying CDSeq2 reference aggregation.
#' @param reference_proportion Cell-type proportion-prior setting passed to
#'   `CDSeq2::cdseq2()`.
#' @param reference_precision Positive scalar defining the global precision
#'   around which cell-type-specific values are calibrated.
#' @param precision_min_cells_per_donor_type Minimum number of reference cells
#'   required for a donor-by-cell-type pseudobulk to inform precision.
#' @param precision_min_donors Minimum number of eligible reference donors
#'   required to estimate a cell-type precision. Cell types below this threshold
#'   use `reference_precision` unchanged.
#' @param precision_raw_bounds Positive lower and upper bounds used when
#'   optimizing the raw Dirichlet--multinomial concentration.
#' @param precision_multiplier_bounds Allowed lower and upper multipliers of
#'   `reference_precision` after relative calibration. The defaults constrain
#'   the first trial to `[0.3, 3]` times the CDSeq2 baseline.
#' @param precision_shrinkage_target_donors Number of eligible donors at which
#'   the nominal reliability weight reaches one before applying
#'   `precision_shrinkage_weight_cap`.
#' @param precision_shrinkage_weight_cap Maximum weight assigned to the raw
#'   reference-donor dispersion estimate. Values below one retain shrinkage
#'   toward the global CDSeq2 precision.
#' @param precision_profile_pseudocount Nonnegative pseudocount used only when
#'   estimating the cell-type mean profile from donor pseudobulks.
#' @param precision_profile_floor Positive floor applied to the estimated mean
#'   profile before evaluating the Dirichlet--multinomial likelihood.
#' @param precision_mean_aggregation Either `"equal_donor"` or
#'   `"pooled_counts"` for estimating the mean profile used by the precision
#'   likelihood. The default prevents high-depth donors from dominating the
#'   dispersion calibration.
#' @param precision_fallback Either `"global"`, which retains the global
#'   precision when estimation is unsupported, or `"error"`.
#' @param smooth_ref Whether to smooth the CDSeq2 reference profile.
#' @param smoothing_alpha Reference-smoothing pseudocount passed to CDSeq2.
#' @param smoothing_scope Reference-smoothing scope passed to CDSeq2.
#' @param genes_use Optional genes passed to CDSeq2 and used for precision
#'   calibration.
#' @param exclude_genes Optional genes excluded from CDSeq2 and precision
#'   calibration.
#' @param drop_all_zero_genes Whether to omit genes with no counts from the
#'   calibration gene universe, matching the CDSeq2 setting.
#' @param aggregate_to Optional annotation aggregation passed to CDSeq2. The
#'   current donor-precision implementation requires `NULL` so that precision
#'   labels correspond directly to reference cell-type labels.
#' @param rng RNG implementation selector passed to CDSeq2.
#' @param save_last_n_mc Number of retained MCMC draws.
#' @param n_chains Number of chains passed to CDSeq2 when supported.
#' @param seed MCMC seed passed to CDSeq2 when supported.
#' @param primary_gep_estimator CDSeq2 expression estimator. The default is
#'   `"rao_blackwell"`.
#' @param return_posterior_sd Whether CDSeq2 should return posterior SDs.
#' @param return_mcmc_arrays Whether CDSeq2 should retain MCMC arrays.
#' @param return_phi_draw_means Whether CDSeq2 should return normalized-profile
#'   draw means.
#' @param return_assignments Whether CDSeq2 should retain assignment output.
#' @param return_aggregated_arrays Whether CDSeq2 should retain aggregated
#'   arrays.
#' @param low_information_depth Optional low-information depth threshold passed
#'   to CDSeq2.
#' @param low_information_prior_weight Optional low-information prior weight
#'   passed to CDSeq2.
#' @param reference_correlation_warning Correlation threshold for warning about
#'   highly similar reference profiles.
#' @param adjust_prior Whether to apply CDSeq2's optional prior adjustment.
#' @param adjust_prior_weight Weight used by CDSeq2's optional prior adjustment.
#' @param verbose Print CDSeq4 and CDSeq2 progress messages.
#' @param ... Additional named arguments passed to `CDSeq2::cdseq2()`.
#'
#' @return A CDSeq2-compatible fit with class `cdseq4_fit` prepended. Existing
#'   components such as `cdseq_result` are preserved. Added components include
#'   `global_reference_precision`, `cell_type_reference_precision`,
#'   `precision_diagnostics`, `precision_calibration_settings`,
#'   `precision_strategy`, and `reference_metadata_audit`.
#'
#' @examples
#' \dontrun{
#' data(cdseq_example_data, package = "CDSeq2")
#'
#' fit4 <- cdseq4(
#'   bulk_data = cdseq_example_data$bulk_data,
#'   reference_scRNA_count = cdseq_example_data$reference_scRNA_count,
#'   reference_scRNA_annotation =
#'     cdseq_example_data$reference_scRNA_annotation,
#'   reference_metadata = cdseq_example_data$reference_metadata,
#'   reference_precision = 1e6,
#'   mcmc_iterations = 1000,
#'   save_last_n_mc = 500,
#'   cpu_number = 6
#' )
#'
#' fit4$cell_type_reference_precision
#' fit4$precision_diagnostics
#' }
#'
#' @export
cdseq4 <- function(
    bulk_data,
    mcmc_iterations = 700,
    cpu_number = 1,
    reference_scRNA_count = NULL,
    reference_scRNA_annotation = NULL,
    reference_metadata = NULL,
    reference_scRNA_donor = NULL,
    reference_cell_id_col = "cell_id",
    reference_donor_col = "reference_donor",
    reference_cell_type_col = "cell_type",
    reference_aggregation = c("pooled_counts", "equal_donor"),
    reference_min_cells_per_donor_type = 1L,
    reference_proportion = 1,
    reference_precision = 1e6,
    precision_min_cells_per_donor_type = 2L,
    precision_min_donors = 5L,
    precision_raw_bounds = c(1, 1e9),
    precision_multiplier_bounds = c(0.3, 3),
    precision_shrinkage_target_donors = 12,
    precision_shrinkage_weight_cap = 0.75,
    precision_profile_pseudocount = 0.5,
    precision_profile_floor = 1e-12,
    precision_mean_aggregation = c("equal_donor", "pooled_counts"),
    precision_fallback = c("global", "error"),
    smooth_ref = TRUE,
    smoothing_alpha = 1,
    smoothing_scope = c("aggregate", "per_cell"),
    genes_use = NULL,
    exclude_genes = NULL,
    drop_all_zero_genes = TRUE,
    aggregate_to = NULL,
    rng = 1,
    save_last_n_mc = 100,
    n_chains = 1L,
    seed = 123L,
    primary_gep_estimator = c("rao_blackwell", "legacy"),
    return_posterior_sd = TRUE,
    return_mcmc_arrays = TRUE,
    return_phi_draw_means = FALSE,
    return_assignments = TRUE,
    return_aggregated_arrays = FALSE,
    low_information_depth = NULL,
    low_information_prior_weight = NULL,
    reference_correlation_warning = 0.98,
    adjust_prior = FALSE,
    adjust_prior_weight = 0.2,
    verbose = FALSE,
    ...
) {
  reference_aggregation <- match.arg(reference_aggregation)
  precision_mean_aggregation <- match.arg(precision_mean_aggregation)
  precision_fallback <- match.arg(precision_fallback)
  smoothing_scope <- match.arg(smoothing_scope)
  primary_gep_estimator <- match.arg(primary_gep_estimator)

  .cdseq4_validate_scalar_settings(
    reference_precision = reference_precision,
    precision_min_cells_per_donor_type =
      precision_min_cells_per_donor_type,
    precision_min_donors = precision_min_donors,
    precision_raw_bounds = precision_raw_bounds,
    precision_multiplier_bounds = precision_multiplier_bounds,
    precision_shrinkage_target_donors =
      precision_shrinkage_target_donors,
    precision_shrinkage_weight_cap = precision_shrinkage_weight_cap,
    precision_profile_pseudocount = precision_profile_pseudocount,
    precision_profile_floor = precision_profile_floor,
    mcmc_iterations = mcmc_iterations,
    save_last_n_mc = save_last_n_mc
  )

  if (!is.null(aggregate_to)) {
    stop(
      "CDSeq4 donor-calibrated precision currently requires aggregate_to = NULL. ",
      "Estimate precision after constructing the final cell-type labels.",
      call. = FALSE
    )
  }

  resolved <- .cdseq4_resolve_reference_metadata(
    reference_scRNA_count = reference_scRNA_count,
    reference_scRNA_annotation = reference_scRNA_annotation,
    reference_metadata = reference_metadata,
    reference_scRNA_donor = reference_scRNA_donor,
    reference_cell_id_col = reference_cell_id_col,
    reference_donor_col = reference_donor_col,
    reference_cell_type_col = reference_cell_type_col
  )
  reference_scRNA_annotation <- resolved$annotation
  reference_scRNA_donor <- resolved$donor

  precision_gene_names <- .cdseq4_select_precision_genes(
    bulk_data = bulk_data,
    reference_scRNA_count = reference_scRNA_count,
    genes_use = genes_use,
    exclude_genes = exclude_genes,
    drop_all_zero_genes = drop_all_zero_genes
  )

  precision_fit <- .cdseq4_estimate_celltype_precisions(
    reference_scRNA_count = reference_scRNA_count[
      precision_gene_names,
      resolved$cell_ids,
      drop = FALSE
    ],
    reference_annotation = reference_scRNA_annotation,
    reference_donor = reference_scRNA_donor,
    global_precision = as.numeric(reference_precision),
    min_cells_per_donor_type =
      as.integer(precision_min_cells_per_donor_type),
    min_donors = as.integer(precision_min_donors),
    raw_bounds = as.numeric(precision_raw_bounds),
    multiplier_bounds = as.numeric(precision_multiplier_bounds),
    shrinkage_target_donors =
      as.numeric(precision_shrinkage_target_donors),
    shrinkage_weight_cap = as.numeric(precision_shrinkage_weight_cap),
    profile_pseudocount = as.numeric(precision_profile_pseudocount),
    profile_floor = as.numeric(precision_profile_floor),
    mean_aggregation = precision_mean_aggregation,
    fallback = precision_fallback
  )

  cell_type_precision <- precision_fit$cell_type_precision
  if (isTRUE(verbose)) {
    message("CDSeq4 donor-calibrated expression-prior precision:")
    for (cell_type in names(cell_type_precision)) {
      message(
        "  ", cell_type, ": ",
        format(cell_type_precision[[cell_type]], scientific = TRUE)
      )
    }
  }

  if (!requireNamespace("CDSeq2", quietly = TRUE)) {
    stop("The CDSeq2 package is required by cdseq4().", call. = FALSE)
  }
  if (!"cdseq2" %in% getNamespaceExports("CDSeq2")) {
    stop(
      "The installed CDSeq2 package does not export cdseq2().",
      call. = FALSE
    )
  }
  cdseq2_fn <- getExportedValue("CDSeq2", "cdseq2")
  supported_cdseq2_args <- names(formals(cdseq2_fn))
  if (!"reference_precision" %in% supported_cdseq2_args) {
    stop(
      "The installed CDSeq2::cdseq2() does not expose reference_precision.",
      call. = FALSE
    )
  }

  explicit_args <- list(
    bulk_data = bulk_data,
    mcmc_iterations = as.integer(mcmc_iterations),
    cpu_number = as.integer(cpu_number),
    reference_scRNA_count = reference_scRNA_count,
    reference_scRNA_annotation = reference_scRNA_annotation,
    reference_scRNA_donor = reference_scRNA_donor,
    reference_aggregation = reference_aggregation,
    reference_min_cells_per_donor_type =
      as.integer(reference_min_cells_per_donor_type),
    reference_proportion = reference_proportion,
    reference_precision = cell_type_precision,
    smooth_ref = smooth_ref,
    smoothing_alpha = smoothing_alpha,
    smoothing_scope = smoothing_scope,
    genes_use = genes_use,
    exclude_genes = exclude_genes,
    drop_all_zero_genes = drop_all_zero_genes,
    aggregate_to = aggregate_to,
    rng = as.integer(rng),
    save_last_n_mc = as.integer(save_last_n_mc),
    n_chains = as.integer(n_chains),
    seed = as.integer(seed),
    primary_gep_estimator = primary_gep_estimator,
    return_posterior_sd = return_posterior_sd,
    return_mcmc_arrays = return_mcmc_arrays,
    return_phi_draw_means = return_phi_draw_means,
    return_assignments = return_assignments,
    return_aggregated_arrays = return_aggregated_arrays,
    low_information_depth = low_information_depth,
    low_information_prior_weight = low_information_prior_weight,
    reference_correlation_warning = reference_correlation_warning,
    adjust_prior = adjust_prior,
    adjust_prior_weight = adjust_prior_weight,
    verbose = verbose
  )

  extra_args <- list(...)
  if (length(extra_args) > 0L) {
    if (is.null(names(extra_args)) || any(names(extra_args) == "")) {
      stop("All arguments supplied through ... must be named.", call. = FALSE)
    }
    if (anyDuplicated(names(extra_args))) {
      stop("Arguments supplied through ... contain duplicate names.", call. = FALSE)
    }
    duplicated_args <- intersect(names(extra_args), names(explicit_args))
    if (length(duplicated_args) > 0L) {
      stop(
        "Arguments supplied through ... duplicate explicit CDSeq4 arguments: ",
        paste(duplicated_args, collapse = ", "),
        call. = FALSE
      )
    }
  }

  call_args <- c(explicit_args, extra_args)
  supported <- names(formals(cdseq2_fn))
  unsupported <- setdiff(names(call_args), supported)
  unsupported_nonnull <- unsupported[
    !vapply(call_args[unsupported], is.null, logical(1))
  ]
  if (length(unsupported_nonnull) > 0L) {
    stop(
      "The installed CDSeq2::cdseq2() does not support: ",
      paste(unsupported_nonnull, collapse = ", "),
      call. = FALSE
    )
  }
  call_args <- call_args[intersect(names(call_args), supported)]

  fit_timing <- system.time({
    fit <- do.call(cdseq2_fn, call_args)
  })

  if (!is.list(fit)) {
    stop("CDSeq2::cdseq2() returned a non-list object.", call. = FALSE)
  }

  native_precision_audit <- .cdseq4_audit_native_precision(
    fit = fit,
    expected_precision = cell_type_precision
  )
  if (!isTRUE(native_precision_audit$verified)) {
    stop(
      "CDSeq4 passed a named cell-type-specific precision vector to ",
      "CDSeq2::cdseq2(), but the returned fit did not expose enough prior ",
      "information to verify that the vector was applied. ",
      "No CDSeq4 performance claim should be made from this fit.",
      call. = FALSE
    )
  }

  fit$global_reference_precision <- as.numeric(reference_precision)
  fit$cell_type_reference_precision <- cell_type_precision
  fit$selected_reference_precision <- cell_type_precision
  fit$precision_strategy <- "reference_donor_relative_dirichlet_multinomial"
  fit$precision_calibration_settings <- list(
    min_cells_per_donor_type =
      as.integer(precision_min_cells_per_donor_type),
    min_donors = as.integer(precision_min_donors),
    raw_bounds = as.numeric(precision_raw_bounds),
    multiplier_bounds = as.numeric(precision_multiplier_bounds),
    shrinkage_target_donors =
      as.numeric(precision_shrinkage_target_donors),
    shrinkage_weight_cap =
      as.numeric(precision_shrinkage_weight_cap),
    profile_pseudocount = as.numeric(precision_profile_pseudocount),
    profile_floor = as.numeric(precision_profile_floor),
    mean_aggregation = precision_mean_aggregation,
    fallback = precision_fallback,
    reference_aggregation = reference_aggregation
  )
  fit$precision_diagnostics <- precision_fit$diagnostics
  fit$weighted_geometric_mean_reference_precision <-
    precision_fit$weighted_geometric_mean_precision
  fit$precision_reference_profiles <- precision_fit$mean_profiles
  fit$precision_gene_names <- precision_gene_names
  fit$reference_metadata_audit <- resolved$audit
  fit$cdseq4_engine_audit <- c(
    list(
      adapter = "native named reference_precision vector passed to CDSeq2::cdseq2",
      final_fit_elapsed_seconds = unname(fit_timing[["elapsed"]])
    ),
    native_precision_audit
  )
  fit$cdseq4_call <- match.call()

  class(fit) <- unique(c("cdseq4_fit", class(fit)))
  fit
}

#' Print a CDSeq4 fit
#'
#' @param x A `cdseq4_fit` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.cdseq4_fit <- function(x, ...) {
  cat("CDSeq4 fit\n")
  cat("  Precision strategy: ", x$precision_strategy, "\n", sep = "")
  if (!is.null(x$cell_type_reference_precision)) {
    precision_text <- paste(
      names(x$cell_type_reference_precision),
      format(x$cell_type_reference_precision, scientific = TRUE),
      sep = "=",
      collapse = "; "
    )
    cat("  Cell-type precision: ", precision_text, "\n", sep = "")
  }
  if (!is.null(x$cdseq4_engine_audit$final_fit_elapsed_seconds)) {
    cat(
      "  Elapsed seconds: ",
      round(x$cdseq4_engine_audit$final_fit_elapsed_seconds, 2),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

.cdseq4_validate_scalar_settings <- function(
    reference_precision,
    precision_min_cells_per_donor_type,
    precision_min_donors,
    precision_raw_bounds,
    precision_multiplier_bounds,
    precision_shrinkage_target_donors,
    precision_shrinkage_weight_cap,
    precision_profile_pseudocount,
    precision_profile_floor,
    mcmc_iterations,
    save_last_n_mc
) {
  if (
    length(reference_precision) != 1L ||
    !is.finite(reference_precision) || reference_precision <= 0
  ) {
    stop("reference_precision must be one positive finite number.", call. = FALSE)
  }
  if (
    length(precision_min_cells_per_donor_type) != 1L ||
    !is.finite(precision_min_cells_per_donor_type) ||
    precision_min_cells_per_donor_type < 1
  ) {
    stop(
      "precision_min_cells_per_donor_type must be at least one.",
      call. = FALSE
    )
  }
  if (
    length(precision_min_donors) != 1L ||
    !is.finite(precision_min_donors) || precision_min_donors < 2
  ) {
    stop("precision_min_donors must be at least two.", call. = FALSE)
  }
  if (
    length(precision_raw_bounds) != 2L ||
    any(!is.finite(precision_raw_bounds)) ||
    any(precision_raw_bounds <= 0) ||
    precision_raw_bounds[[1L]] >= precision_raw_bounds[[2L]]
  ) {
    stop("precision_raw_bounds must be increasing and positive.", call. = FALSE)
  }
  if (
    length(precision_multiplier_bounds) != 2L ||
    any(!is.finite(precision_multiplier_bounds)) ||
    any(precision_multiplier_bounds <= 0) ||
    precision_multiplier_bounds[[1L]] > 1 ||
    precision_multiplier_bounds[[2L]] < 1
  ) {
    stop(
      "precision_multiplier_bounds must be positive and contain one.",
      call. = FALSE
    )
  }
  if (
    length(precision_shrinkage_target_donors) != 1L ||
    !is.finite(precision_shrinkage_target_donors) ||
    precision_shrinkage_target_donors <= 0
  ) {
    stop("precision_shrinkage_target_donors must be positive.", call. = FALSE)
  }
  if (
    length(precision_shrinkage_weight_cap) != 1L ||
    !is.finite(precision_shrinkage_weight_cap) ||
    precision_shrinkage_weight_cap < 0 ||
    precision_shrinkage_weight_cap > 1
  ) {
    stop("precision_shrinkage_weight_cap must be in [0, 1].", call. = FALSE)
  }
  if (
    length(precision_profile_pseudocount) != 1L ||
    !is.finite(precision_profile_pseudocount) ||
    precision_profile_pseudocount < 0
  ) {
    stop("precision_profile_pseudocount must be nonnegative.", call. = FALSE)
  }
  if (
    length(precision_profile_floor) != 1L ||
    !is.finite(precision_profile_floor) ||
    precision_profile_floor <= 0
  ) {
    stop("precision_profile_floor must be positive.", call. = FALSE)
  }
  if (
    length(mcmc_iterations) != 1L || length(save_last_n_mc) != 1L ||
    mcmc_iterations < 1 || save_last_n_mc < 1 ||
    save_last_n_mc > mcmc_iterations
  ) {
    stop("Require mcmc_iterations >= save_last_n_mc >= 1.", call. = FALSE)
  }
  invisible(TRUE)
}

.cdseq4_resolve_reference_metadata <- function(
    reference_scRNA_count,
    reference_scRNA_annotation,
    reference_metadata,
    reference_scRNA_donor,
    reference_cell_id_col,
    reference_donor_col,
    reference_cell_type_col
) {
  if (is.null(reference_scRNA_count) || length(dim(reference_scRNA_count)) != 2L) {
    stop("reference_scRNA_count must be a gene-by-cell matrix.", call. = FALSE)
  }
  cell_ids <- colnames(reference_scRNA_count)
  if (is.null(cell_ids) || anyNA(cell_ids) || any(cell_ids == "")) {
    stop("reference_scRNA_count must have nonempty cell column names.", call. = FALSE)
  }
  if (anyDuplicated(cell_ids)) {
    stop("Reference cell IDs must be unique.", call. = FALSE)
  }

  metadata_aligned <- NULL
  metadata_annotation <- NULL
  metadata_donor <- NULL

  if (!is.null(reference_metadata)) {
    metadata <- as.data.frame(reference_metadata, stringsAsFactors = FALSE)
    metadata_ids <- if (reference_cell_id_col %in% colnames(metadata)) {
      as.character(metadata[[reference_cell_id_col]])
    } else if (!is.null(rownames(metadata))) {
      as.character(rownames(metadata))
    } else {
      stop(
        "reference_metadata must contain '", reference_cell_id_col,
        "' or informative row names.",
        call. = FALSE
      )
    }
    if (anyNA(metadata_ids) || any(metadata_ids == "") || anyDuplicated(metadata_ids)) {
      stop("reference_metadata cell IDs must be complete and unique.", call. = FALSE)
    }
    match_index <- match(cell_ids, metadata_ids)
    if (anyNA(match_index)) {
      stop(
        "Some reference count columns are absent from reference_metadata: ",
        paste(utils::head(cell_ids[is.na(match_index)], 8L), collapse = ", "),
        call. = FALSE
      )
    }
    metadata_aligned <- metadata[match_index, , drop = FALSE]
    rownames(metadata_aligned) <- cell_ids

    if (!reference_donor_col %in% colnames(metadata_aligned)) {
      stop(
        "reference_metadata lacks donor column '", reference_donor_col, "'.",
        call. = FALSE
      )
    }
    if (!reference_cell_type_col %in% colnames(metadata_aligned)) {
      stop(
        "reference_metadata lacks cell-type column '",
        reference_cell_type_col, "'.",
        call. = FALSE
      )
    }
    metadata_donor <- as.character(metadata_aligned[[reference_donor_col]])
    metadata_annotation <-
      as.character(metadata_aligned[[reference_cell_type_col]])
  }

  annotation <- reference_scRNA_annotation
  if (!is.null(annotation)) {
    annotation <- as.character(annotation)
    if (!is.null(names(reference_scRNA_annotation))) {
      annotation_names <- names(reference_scRNA_annotation)
      if (setequal(annotation_names, cell_ids)) {
        annotation <- annotation[match(cell_ids, annotation_names)]
      }
    }
    if (length(annotation) != length(cell_ids)) {
      stop(
        "reference_scRNA_annotation must contain one label per reference cell.",
        call. = FALSE
      )
    }
  } else {
    annotation <- metadata_annotation
  }
  if (is.null(annotation)) {
    stop(
      "Supply reference_scRNA_annotation or reference_metadata cell types.",
      call. = FALSE
    )
  }
  if (anyNA(annotation) || any(annotation == "")) {
    stop("Reference cell-type labels must be complete.", call. = FALSE)
  }
  names(annotation) <- cell_ids

  if (!is.null(metadata_annotation)) {
    mismatch <- annotation != metadata_annotation
    if (any(mismatch)) {
      mismatch_ids <- cell_ids[mismatch]
      stop(
        "reference_scRNA_annotation disagrees with reference_metadata for: ",
        paste(utils::head(mismatch_ids, 8L), collapse = ", "),
        call. = FALSE
      )
    }
  }

  donor <- reference_scRNA_donor
  if (!is.null(donor)) {
    donor <- as.character(donor)
    if (!is.null(names(reference_scRNA_donor))) {
      donor_names <- names(reference_scRNA_donor)
      if (setequal(donor_names, cell_ids)) {
        donor <- donor[match(cell_ids, donor_names)]
      }
    }
    if (length(donor) != length(cell_ids)) {
      stop(
        "reference_scRNA_donor must contain one donor per reference cell.",
        call. = FALSE
      )
    }
  } else {
    donor <- metadata_donor
  }
  if (is.null(donor)) {
    stop(
      "Donor-calibrated precision requires reference_metadata or ",
      "reference_scRNA_donor.",
      call. = FALSE
    )
  }
  if (anyNA(donor) || any(donor == "")) {
    stop("Reference donor labels must be complete.", call. = FALSE)
  }
  names(donor) <- cell_ids

  if (!is.null(metadata_donor) && any(donor != metadata_donor)) {
    mismatch_ids <- cell_ids[donor != metadata_donor]
    stop(
      "reference_scRNA_donor disagrees with reference_metadata for: ",
      paste(utils::head(mismatch_ids, 8L), collapse = ", "),
      call. = FALSE
    )
  }

  audit <- data.frame(
    item = c(
      "reference_cells",
      "reference_donors",
      "reference_cell_types",
      "metadata_alignment"
    ),
    value = c(
      length(cell_ids),
      length(unique(donor)),
      length(unique(annotation)),
      if (is.null(reference_metadata)) "direct vectors" else "matched by cell ID"
    ),
    stringsAsFactors = FALSE
  )

  list(
    cell_ids = cell_ids,
    annotation = annotation,
    donor = donor,
    metadata = metadata_aligned,
    audit = audit
  )
}

.cdseq4_select_precision_genes <- function(
    bulk_data,
    reference_scRNA_count,
    genes_use,
    exclude_genes,
    drop_all_zero_genes
) {
  if (is.null(rownames(bulk_data)) || is.null(rownames(reference_scRNA_count))) {
    stop("Bulk and reference matrices must have gene row names.", call. = FALSE)
  }
  genes <- rownames(bulk_data)[
    rownames(bulk_data) %in% rownames(reference_scRNA_count)
  ]
  if (!is.null(genes_use)) {
    genes <- genes[genes %in% as.character(genes_use)]
  }
  if (!is.null(exclude_genes)) {
    genes <- setdiff(genes, as.character(exclude_genes))
  }
  if (isTRUE(drop_all_zero_genes) && length(genes) > 0L) {
    bulk_positive <- rowSums(as.matrix(bulk_data[genes, , drop = FALSE])) > 0
    reference_positive <-
      rowSums(as.matrix(reference_scRNA_count[genes, , drop = FALSE])) > 0
    genes <- genes[bulk_positive & reference_positive]
  }
  if (length(genes) < 2L) {
    stop("Fewer than two genes remain for precision calibration.", call. = FALSE)
  }
  genes
}

.cdseq4_estimate_celltype_precisions <- function(
    reference_scRNA_count,
    reference_annotation,
    reference_donor,
    global_precision,
    min_cells_per_donor_type,
    min_donors,
    raw_bounds,
    multiplier_bounds,
    shrinkage_target_donors,
    shrinkage_weight_cap,
    profile_pseudocount,
    profile_floor,
    mean_aggregation,
    fallback
) {
  counts <- as.matrix(reference_scRNA_count)
  storage.mode(counts) <- "double"
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("Reference counts must be finite and nonnegative.", call. = FALSE)
  }

  cell_types <- unique(as.character(reference_annotation))
  estimate_rows <- vector("list", length(cell_types))
  names(estimate_rows) <- cell_types
  mean_profiles <- vector("list", length(cell_types))
  names(mean_profiles) <- cell_types

  for (cell_type in cell_types) {
    cell_index <- which(reference_annotation == cell_type)
    donor_for_type <- reference_donor[cell_index]
    donor_cell_counts <- table(donor_for_type)
    eligible_donors <- names(donor_cell_counts)[
      donor_cell_counts >= min_cells_per_donor_type
    ]

    donor_count_matrix <- if (length(eligible_donors) > 0L) {
      donor_by_gene <- vapply(
        eligible_donors,
        function(donor_id) {
          selected_cells <- cell_index[donor_for_type == donor_id]
          rowSums(counts[, selected_cells, drop = FALSE])
        },
        numeric(nrow(counts))
      )
      donor_by_gene <- t(donor_by_gene)
      colnames(donor_by_gene) <- rownames(counts)
      rownames(donor_by_gene) <- eligible_donors
      donor_by_gene
    } else {
      matrix(numeric(0), nrow = 0L, ncol = nrow(counts))
    }

    positive_depth <- if (nrow(donor_count_matrix) > 0L) {
      rowSums(donor_count_matrix) > 0
    } else {
      logical(0)
    }
    donor_count_matrix <- donor_count_matrix[positive_depth, , drop = FALSE]
    eligible_donors <- rownames(donor_count_matrix)
    eligible_cell_counts <- as.numeric(donor_cell_counts[eligible_donors])

    status <- "ok"
    raw_precision <- NA_real_
    objective <- NA_real_
    lower_boundary <- FALSE
    upper_boundary <- FALSE
    mean_profile <- rep(NA_real_, nrow(counts))
    names(mean_profile) <- rownames(counts)

    if (length(eligible_donors) < min_donors) {
      status <- "insufficient_eligible_donors"
    } else {
      donor_depth <- rowSums(donor_count_matrix)
      donor_profiles <- sweep(
        donor_count_matrix + profile_pseudocount,
        1L,
        donor_depth + profile_pseudocount * ncol(donor_count_matrix),
        "/"
      )
      if (identical(mean_aggregation, "equal_donor")) {
        mean_profile <- colMeans(donor_profiles)
      } else {
        pooled_counts <- colSums(donor_count_matrix)
        mean_profile <- (pooled_counts + profile_pseudocount) /
          (sum(pooled_counts) + profile_pseudocount * length(pooled_counts))
      }
      mean_profile <- pmax(mean_profile, profile_floor)
      mean_profile <- mean_profile / sum(mean_profile)

      negative_log_likelihood <- function(log_precision) {
        precision <- exp(log_precision)
        alpha <- precision * mean_profile
        value <- 0
        for (donor_index in seq_len(nrow(donor_count_matrix))) {
          donor_counts <- donor_count_matrix[donor_index, ]
          donor_total <- sum(donor_counts)
          value <- value +
            lgamma(precision) - lgamma(precision + donor_total) +
            sum(lgamma(alpha + donor_counts) - lgamma(alpha))
        }
        if (!is.finite(value)) return(.Machine$double.xmax / 100)
        -value
      }

      optimization <- tryCatch(
        stats::optimize(
          f = negative_log_likelihood,
          interval = log(raw_bounds),
          tol = .Machine$double.eps^0.25
        ),
        error = function(e) e
      )
      if (inherits(optimization, "error")) {
        status <- paste0("optimization_error: ", conditionMessage(optimization))
      } else if (!is.finite(optimization$minimum)) {
        status <- "nonfinite_optimization_result"
      } else {
        raw_precision <- exp(optimization$minimum)
        objective <- optimization$objective
        log_bounds <- log(raw_bounds)
        boundary_tolerance <- 0.01 * diff(log_bounds)
        lower_boundary <-
          optimization$minimum <= log_bounds[[1L]] + boundary_tolerance
        upper_boundary <-
          optimization$minimum >= log_bounds[[2L]] - boundary_tolerance
      }
    }

    mean_profiles[[cell_type]] <- mean_profile
    estimate_rows[[cell_type]] <- data.frame(
      cell_type = cell_type,
      total_reference_cells = length(cell_index),
      total_reference_donors = length(unique(donor_for_type)),
      eligible_reference_donors = length(eligible_donors),
      min_cells_per_eligible_donor = if (length(eligible_cell_counts)) {
        min(eligible_cell_counts)
      } else {
        NA_real_
      },
      median_cells_per_eligible_donor = if (length(eligible_cell_counts)) {
        stats::median(eligible_cell_counts)
      } else {
        NA_real_
      },
      max_cells_per_eligible_donor = if (length(eligible_cell_counts)) {
        max(eligible_cell_counts)
      } else {
        NA_real_
      },
      median_pseudobulk_umi_depth = if (nrow(donor_count_matrix)) {
        stats::median(rowSums(donor_count_matrix))
      } else {
        NA_real_
      },
      raw_dm_precision = raw_precision,
      optimization_objective = objective,
      raw_at_lower_boundary = lower_boundary,
      raw_at_upper_boundary = upper_boundary,
      estimation_status = status,
      stringsAsFactors = FALSE
    )
  }

  diagnostics <- do.call(rbind, estimate_rows)
  rownames(diagnostics) <- NULL
  valid <- diagnostics$estimation_status == "ok" &
    is.finite(diagnostics$raw_dm_precision) &
    diagnostics$raw_dm_precision > 0

  if (!any(valid)) {
    if (identical(fallback, "error")) {
      stop(
        "No cell type had enough reference donors for precision estimation.",
        call. = FALSE
      )
    }
    diagnostics$reliability_weight <- 0
    diagnostics$relative_log_raw_precision <- NA_real_
    diagnostics$uncapped_calibrated_precision <- global_precision
    diagnostics$bounded_log_multiplier <- 0
    diagnostics$calibration_shift_after_bounding <- 0
    diagnostics$final_reference_precision <- global_precision
    diagnostics$precision_multiplier <- 1
    diagnostics$weighted_geometric_mean_multiplier <- 1
    diagnostics$fallback_to_global <- TRUE
  } else {
    center_log_precision <- stats::weighted.mean(
      log(diagnostics$raw_dm_precision[valid]),
      w = diagnostics$eligible_reference_donors[valid]
    )
    diagnostics$reliability_weight <- 0
    diagnostics$reliability_weight[valid] <- pmin(
      shrinkage_weight_cap,
      diagnostics$eligible_reference_donors[valid] /
        shrinkage_target_donors
    )
    diagnostics$relative_log_raw_precision <- NA_real_
    diagnostics$relative_log_raw_precision[valid] <-
      log(diagnostics$raw_dm_precision[valid]) - center_log_precision

    adjusted_relative_log <- rep(0, nrow(diagnostics))
    adjusted_relative_log[valid] <-
      diagnostics$reliability_weight[valid] *
      diagnostics$relative_log_raw_precision[valid]

    # Recenter after reliability shrinkage so the donor-calibrated strategy
    # redistributes, rather than globally weakens, the prespecified precision.
    adjusted_center <- stats::weighted.mean(
      adjusted_relative_log[valid],
      w = diagnostics$eligible_reference_donors[valid]
    )
    adjusted_relative_log[valid] <-
      adjusted_relative_log[valid] - adjusted_center

    diagnostics$uncapped_calibrated_precision <- global_precision
    diagnostics$uncapped_calibrated_precision[valid] <-
      global_precision * exp(adjusted_relative_log[valid])

    # Apply multiplier bounds while preserving the weighted geometric mean
    # of the valid cell-type multipliers at one. This prevents clipping from
    # silently shifting the complete strategy toward globally weaker or
    # stronger regularization.
    lower_log_multiplier <- log(multiplier_bounds[[1L]])
    upper_log_multiplier <- log(multiplier_bounds[[2L]])
    valid_weights <- diagnostics$eligible_reference_donors[valid]
    valid_relative_log <- adjusted_relative_log[valid]

    bounded_center_function <- function(shift) {
      bounded <- pmin(
        upper_log_multiplier,
        pmax(lower_log_multiplier, valid_relative_log + shift)
      )
      stats::weighted.mean(bounded, w = valid_weights)
    }

    lower_shift <- lower_log_multiplier - max(valid_relative_log) - 10
    upper_shift <- upper_log_multiplier - min(valid_relative_log) + 10
    lower_value <- bounded_center_function(lower_shift)
    upper_value <- bounded_center_function(upper_shift)

    if (abs(lower_value) < 1e-12) {
      calibration_shift <- lower_shift
    } else if (abs(upper_value) < 1e-12) {
      calibration_shift <- upper_shift
    } else if (lower_value < 0 && upper_value > 0) {
      calibration_shift <- stats::uniroot(
        bounded_center_function,
        interval = c(lower_shift, upper_shift),
        tol = 1e-10
      )$root
    } else {
      stop(
        "Could not center the bounded cell-type precision multipliers.",
        call. = FALSE
      )
    }

    bounded_log_multiplier <- rep(0, nrow(diagnostics))
    bounded_log_multiplier[valid] <- pmin(
      upper_log_multiplier,
      pmax(
        lower_log_multiplier,
        valid_relative_log + calibration_shift
      )
    )

    diagnostics$bounded_log_multiplier <- bounded_log_multiplier
    diagnostics$calibration_shift_after_bounding <- calibration_shift
    diagnostics$final_reference_precision <-
      global_precision * exp(bounded_log_multiplier)
    diagnostics$final_reference_precision <- signif(
      diagnostics$final_reference_precision,
      4L
    )
    diagnostics$precision_multiplier <-
      diagnostics$final_reference_precision / global_precision
    diagnostics$fallback_to_global <- !valid

    final_geometric_mean_multiplier <- exp(stats::weighted.mean(
      log(diagnostics$precision_multiplier[valid]),
      w = valid_weights
    ))
    diagnostics$weighted_geometric_mean_multiplier <-
      final_geometric_mean_multiplier
  }

  if (identical(fallback, "error") && any(diagnostics$fallback_to_global)) {
    failed <- diagnostics$cell_type[diagnostics$fallback_to_global]
    stop(
      "Precision estimation failed or lacked donor support for: ",
      paste(failed, collapse = ", "),
      call. = FALSE
    )
  }

  cell_type_precision <- diagnostics$final_reference_precision
  names(cell_type_precision) <- diagnostics$cell_type
  cell_type_precision <- cell_type_precision[cell_types]

  list(
    cell_type_precision = cell_type_precision,
    diagnostics = diagnostics,
    mean_profiles = mean_profiles,
    weighted_geometric_mean_precision = if (any(valid)) {
      exp(stats::weighted.mean(
        log(cell_type_precision[diagnostics$cell_type[valid]]),
        w = diagnostics$eligible_reference_donors[valid]
      ))
    } else {
      global_precision
    },
    strategy = "relative donor Dirichlet-multinomial dispersion centered at global precision"
  )
}

.cdseq4_audit_native_precision <- function(
    fit,
    expected_precision,
    tolerance = 1e-6
) {
  expected_precision <- as.numeric(expected_precision) |>
    stats::setNames(names(expected_precision))

  if (
    is.null(names(expected_precision)) ||
    any(names(expected_precision) == "") ||
    anyDuplicated(names(expected_precision))
  ) {
    stop(
      "expected_precision must be uniquely named by cell type.",
      call. = FALSE
    )
  }

  cell_types <- names(expected_precision)
  audit_sources <- character(0)
  observed_precision <- NULL
  beta_cell_type_axis <- NA_character_
  max_abs_precision_error <- NA_real_
  max_relative_precision_error <- NA_real_

  # First use an explicitly returned precision vector when available.
  parameter_precision <- NULL
  if (
    is.list(fit$parameters) &&
    !is.null(fit$parameters$reference_precision)
  ) {
    parameter_precision <- fit$parameters$reference_precision
  } else if (!is.null(fit$reference_precision)) {
    parameter_precision <- fit$reference_precision
  }

  if (!is.null(parameter_precision)) {
    parameter_precision <- as.numeric(parameter_precision) |>
      stats::setNames(names(parameter_precision))

    if (
      length(parameter_precision) == length(expected_precision) &&
      !is.null(names(parameter_precision)) &&
      setequal(names(parameter_precision), cell_types)
    ) {
      parameter_precision <- parameter_precision[cell_types]
      parameter_error <- max(abs(
        parameter_precision - expected_precision
      ))
      parameter_relative_error <- max(
        abs(parameter_precision - expected_precision) /
          pmax(1, abs(expected_precision))
      )
      if (parameter_relative_error <= tolerance) {
        audit_sources <- c(
          audit_sources,
          "returned reference_precision vector"
        )
        observed_precision <- parameter_precision
        max_abs_precision_error <- parameter_error
        max_relative_precision_error <- parameter_relative_error
      } else {
        stop(
          "The precision vector returned by CDSeq2 does not match the ",
          "CDSeq4 donor-calibrated vector.",
          call. = FALSE
        )
      }
    }
  }

  # Independently audit beta whenever CDSeq2 returns it. For a normalized
  # reference profile, the Dirichlet concentration for a cell type is the
  # sum of beta over genes.
  beta <- NULL
  if (is.list(fit$parameters) && !is.null(fit$parameters$beta)) {
    beta <- fit$parameters$beta
  } else if (!is.null(fit$beta)) {
    beta <- fit$beta
  }

  if (!is.null(beta) && length(dim(beta)) == 2L) {
    beta <- as.matrix(beta)
    beta_dimnames <- dimnames(beta)
    beta_rows <- if (is.null(beta_dimnames)) NULL else beta_dimnames[[1L]]
    beta_cols <- if (is.null(beta_dimnames)) NULL else beta_dimnames[[2L]]

    if (!is.null(beta_rows) && setequal(beta_rows, cell_types)) {
      beta_cell_type_axis <- "rows"
      beta_precision <- rowSums(beta[cell_types, , drop = FALSE])
    } else if (!is.null(beta_cols) && setequal(beta_cols, cell_types)) {
      beta_cell_type_axis <- "columns"
      beta_precision <- colSums(beta[, cell_types, drop = FALSE])
    } else if (
      nrow(beta) == length(cell_types) &&
      ncol(beta) != length(cell_types)
    ) {
      beta_cell_type_axis <- "rows_by_dimension"
      beta_precision <- rowSums(beta)
      names(beta_precision) <- cell_types
    } else if (
      ncol(beta) == length(cell_types) &&
      nrow(beta) != length(cell_types)
    ) {
      beta_cell_type_axis <- "columns_by_dimension"
      beta_precision <- colSums(beta)
      names(beta_precision) <- cell_types
    } else {
      beta_precision <- NULL
    }

    if (!is.null(beta_precision)) {
      if (is.null(names(beta_precision))) {
        names(beta_precision) <- cell_types
      }
      beta_precision <- beta_precision[cell_types]
      beta_error <- max(abs(beta_precision - expected_precision))
      beta_relative_error <- max(
        abs(beta_precision - expected_precision) /
          pmax(1, abs(expected_precision))
      )

      if (beta_relative_error > tolerance) {
        stop(
          "The beta matrix returned by CDSeq2 does not have the requested ",
          "CDSeq4 cell-type-specific concentrations. Maximum relative error: ",
          signif(beta_relative_error, 4),
          ".",
          call. = FALSE
        )
      }

      audit_sources <- c(audit_sources, "returned beta matrix")
      observed_precision <- beta_precision
      max_abs_precision_error <- beta_error
      max_relative_precision_error <- beta_relative_error
    }
  }

  verified <- length(audit_sources) > 0L

  list(
    verified = verified,
    verification_source = if (verified) {
      paste(unique(audit_sources), collapse = "; ")
    } else {
      NA_character_
    },
    beta_cell_type_axis = beta_cell_type_axis,
    requested_precision = expected_precision,
    observed_precision = observed_precision,
    max_abs_precision_error = max_abs_precision_error,
    max_relative_precision_error = max_relative_precision_error
  )
}
