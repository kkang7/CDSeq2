#' CDSeq3: direct posterior-profile estimation built on CDSeq2
#'
#' `cdseq3()` is a lightweight wrapper around `CDSeq2::cdseq2()` designed to
#' improve estimation of normalized, individual-level cell-type-specific gene
#' expression profiles (ctsGEPs) without changing the CDSeq2 likelihood or C++
#' Gibbs sampler.
#'
#' The key output change is that the normalized ctsGEP is estimated as the
#' posterior mean of the draw-normalized profile,
#'
#'     E[ Z_gts / sum_g Z_gts | data ],
#'
#' rather than by normalizing the posterior-mean allocated counts,
#'
#'     E[Z_gts | data] / sum_g E[Z_gts | data].
#'
#' These quantities need not be equal because the cell-type-specific profile
#' total can vary across retained MCMC draws. CDSeq3 returns both versions so
#' that the change can be audited directly.
#'
#' CDSeq3 also provides an optional data-only depth-adaptive strategy for
#' selecting one global reference precision. It uses a short pilot fit to
#' estimate the RNA depth assigned to each cell type and chooses a precision
#' that gives a prespecified approximate data weight at a lower quantile of
#' cell-type depths. Benchmark truth is never used for this selection.
#'
#' @param bulk_data Gene-by-sample raw count matrix.
#' @param mcmc_iterations Number of MCMC iterations for the final fit.
#' @param cpu_number Number of workers passed to `CDSeq2::cdseq2()`.
#' @param reference_scRNA_count Gene-by-cell reference count matrix.
#' @param reference_scRNA_annotation Cell-type label for every reference cell.
#' @param reference_scRNA_donor Optional donor label for every reference cell.
#' @param reference_aggregation Reference aggregation strategy. If `NULL`,
#'   `"equal_donor"` is used when donor labels are supplied and
#'   `"pooled_counts"` otherwise.
#' @param reference_min_cells_per_donor_type Minimum cells required for a
#'   donor-cell-type contribution under equal-donor aggregation.
#' @param reference_proportion Prior setting passed to `CDSeq2::cdseq2()`.
#' @param reference_precision Positive fixed or pilot reference precision.
#' @param precision_strategy Either `"fixed"` or `"depth_adaptive"`.
#' @param adaptive_target_data_weight Target heuristic data weight in `(0,1)`.
#'   If effective cell-type depth is `N_t`, the heuristic is
#'   `N_t / (rho + N_t)`.
#' @param adaptive_depth_quantile Quantile of median estimated cell-type RNA
#'   depths used to define the adaptive global precision.
#' @param adaptive_precision_bounds Positive lower and upper bounds for the
#'   selected adaptive precision.
#' @param pilot_mcmc_iterations,pilot_save_last_n_mc Pilot-chain settings.
#' @param pilot_seed Pilot seed. Defaults to `seed + 100000L`.
#' @param smooth_ref,smoothing_alpha,smoothing_scope Reference-smoothing options
#'   passed through when supported.
#' @param genes_use,exclude_genes,drop_all_zero_genes,aggregate_to Gene and
#'   annotation options passed through when supported.
#' @param rng RNG implementation selector passed through to CDSeq2.
#' @param save_last_n_mc Number of final retained MCMC draws.
#' @param n_chains Number of chains passed through when supported.
#' @param seed Final-fit seed.
#' @param proportion_component Posterior RNA-proportion component to extract.
#' @param primary_gep_estimator Base CDSeq2 estimator; defaults to
#'   `"rao_blackwell"`.
#' @param return_posterior_sd Whether to calculate posterior SDs of direct
#'   normalized ctsGEPs across retained draws.
#' @param return_assignments Whether to retain assignment-related output from
#'   `CDSeq2::cdseq2()` when supported.
#' @param return_aggregated_arrays Whether to retain aggregated posterior arrays
#'   from `CDSeq2::cdseq2()` when supported.
#' @param low_information_depth Optional depth threshold for low-information
#'   profiles, passed to `CDSeq2::cdseq2()` when supported.
#' @param low_information_prior_weight Optional prior weight for low-information
#'   profiles, passed to `CDSeq2::cdseq2()` when supported.
#' @param reference_correlation_warning Correlation threshold used to warn about
#'   highly similar reference profiles, passed to `CDSeq2::cdseq2()` when
#'   supported.
#' @param adjust_prior Whether to apply the optional CDSeq2 prior-adjustment
#'   procedure.
#' @param adjust_prior_weight Weight used by the optional CDSeq2 prior-adjustment
#'   procedure.
#' @param keep_cdseq2_fit Keep the final underlying `cdseq2()` fit.
#' @param keep_pilot_fit Keep the pilot fit under adaptive precision.
#' @param verbose Print progress messages.
#' @param ... Additional named arguments passed to `CDSeq2::cdseq2()`.
#'
#' @return An object of class `cdseq3_fit`. Important components are
#'   `ctsGEP`, `ctsGEP_normalized`, `ctsGEP_normalized_posthoc`,
#'   `ctsGEP_normalized_sd`, `rna_proportions`,
#'   `selected_reference_precision`, and `precision_diagnostics`.
#'
#' @examples
#' \dontrun{
#' fit3 <- cdseq3(
#'   bulk_data = bulk_data,
#'   reference_scRNA_count = reference_count,
#'   reference_scRNA_annotation = reference_annotation,
#'   reference_precision = 1e6,
#'   precision_strategy = "fixed",
#'   mcmc_iterations = 1000,
#'   save_last_n_mc = 500,
#'   cpu_number = 6
#' )
#' normalized_ctsGEP <- fit3$ctsGEP_normalized
#' }
#'
#' @export
cdseq3 <- function(
    bulk_data,
    mcmc_iterations = 700,
    cpu_number = 1,
    reference_scRNA_count = NULL,
    reference_scRNA_annotation = NULL,
    reference_scRNA_donor = NULL,
    reference_aggregation = NULL,
    reference_min_cells_per_donor_type = 1L,
    reference_proportion = 1,
    reference_precision = 1e10,
    precision_strategy = c("fixed", "depth_adaptive"),
    adaptive_target_data_weight = 0.50,
    adaptive_depth_quantile = 0.25,
    adaptive_precision_bounds = c(1e4, 1e6),
    pilot_mcmc_iterations = 200L,
    pilot_save_last_n_mc = 100L,
    pilot_seed = NULL,
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
    proportion_component = "estProp",
    primary_gep_estimator = c("rao_blackwell", "legacy"),
    return_posterior_sd = TRUE,
    return_assignments = TRUE,
    return_aggregated_arrays = FALSE,
    low_information_depth = NULL,
    low_information_prior_weight = NULL,
    reference_correlation_warning = 0.98,
    adjust_prior = FALSE,
    adjust_prior_weight = 0.2,
    keep_cdseq2_fit = TRUE,
    keep_pilot_fit = FALSE,
    verbose = FALSE,
    ...
) {
  precision_strategy <- match.arg(precision_strategy)
  primary_gep_estimator <- match.arg(primary_gep_estimator)
  smoothing_scope <- match.arg(smoothing_scope)

  .cdseq3_validate_inputs(
    bulk_data = bulk_data,
    reference_scRNA_count = reference_scRNA_count,
    reference_scRNA_annotation = reference_scRNA_annotation,
    reference_scRNA_donor = reference_scRNA_donor,
    reference_precision = reference_precision,
    mcmc_iterations = mcmc_iterations,
    save_last_n_mc = save_last_n_mc,
    adaptive_target_data_weight = adaptive_target_data_weight,
    adaptive_depth_quantile = adaptive_depth_quantile,
    adaptive_precision_bounds = adaptive_precision_bounds,
    pilot_mcmc_iterations = pilot_mcmc_iterations,
    pilot_save_last_n_mc = pilot_save_last_n_mc,
    proportion_component = proportion_component
  )

  if (is.null(reference_aggregation)) {
    reference_aggregation <- if (is.null(reference_scRNA_donor)) {
      "pooled_counts"
    } else {
      "equal_donor"
    }
  }
  reference_aggregation <- match.arg(
    reference_aggregation,
    c("pooled_counts", "equal_donor")
  )
  if (
    identical(reference_aggregation, "equal_donor") &&
    is.null(reference_scRNA_donor)
  ) {
    stop(
      "reference_scRNA_donor is required when ",
      "reference_aggregation = 'equal_donor'.",
      call. = FALSE
    )
  }

  if (is.null(pilot_seed)) pilot_seed <- as.integer(seed) + 100000L

  extra_args <- list(...)
  if (length(extra_args) > 0L) {
    if (is.null(names(extra_args)) || any(names(extra_args) == "")) {
      stop("All arguments supplied through ... must be named.", call. = FALSE)
    }
    if (anyDuplicated(names(extra_args))) {
      stop("Arguments supplied through ... contain duplicated names.", call. = FALSE)
    }
  }

  common_args <- list(
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
    # Direct posterior-profile estimation requires retained GEP arrays.
    return_mcmc_arrays = TRUE,
    # Ask package versions that support it to retain their own summary too.
    # CDSeq3 does not depend on an undocumented output name; it reconstructs
    # the direct profile from estGEP_array below.
    return_phi_draw_means = TRUE,
    return_assignments = return_assignments,
    return_aggregated_arrays = return_aggregated_arrays,
    low_information_depth = low_information_depth,
    low_information_prior_weight = low_information_prior_weight,
    reference_correlation_warning = reference_correlation_warning,
    adjust_prior = adjust_prior,
    adjust_prior_weight = adjust_prior_weight,
    verbose = verbose
  )

  duplicate_extra <- intersect(names(extra_args), names(common_args))
  if (length(duplicate_extra) > 0L) {
    stop(
      "Arguments supplied through ... duplicate explicit cdseq3 arguments: ",
      paste(duplicate_extra, collapse = ", "),
      call. = FALSE
    )
  }
  common_args <- c(common_args, extra_args)

  pilot_fit <- NULL
  selected_precision <- as.numeric(reference_precision)
  precision_diagnostics <- list(
    strategy = precision_strategy,
    requested_reference_precision = as.numeric(reference_precision),
    selected_reference_precision = as.numeric(reference_precision),
    target_data_weight = NA_real_,
    depth_quantile = NA_real_,
    median_estimated_cell_type_depth = NULL,
    estimated_data_weight_by_cell_type = NULL,
    depth_scale_used = NA_real_,
    unclipped_precision = as.numeric(reference_precision),
    precision_bounds = as.numeric(adaptive_precision_bounds),
    pilot_mcmc_iterations = 0L,
    pilot_save_last_n_mc = 0L,
    pilot_seed = NA_integer_,
    pilot_elapsed_seconds = 0
  )

  if (identical(precision_strategy, "depth_adaptive")) {
    pilot_args <- common_args
    pilot_args$reference_precision <- as.numeric(reference_precision)
    pilot_args$mcmc_iterations <- as.integer(pilot_mcmc_iterations)
    pilot_args$save_last_n_mc <- as.integer(pilot_save_last_n_mc)
    pilot_args$seed <- as.integer(pilot_seed)
    pilot_args$return_posterior_sd <- FALSE
    pilot_args$return_assignments <- FALSE
    pilot_args$return_aggregated_arrays <- FALSE

    if (isTRUE(verbose)) {
      message(
        "CDSeq3 pilot fit at reference precision ",
        format(reference_precision, scientific = TRUE), "."
      )
    }

    pilot_timing <- system.time({
      pilot_fit <- .cdseq3_call_cdseq2(pilot_args)
    })
    pilot_summary <- .cdseq3_summarize_fit(
      fit = pilot_fit,
      bulk_data = bulk_data,
      reference_scRNA_annotation = reference_scRNA_annotation,
      proportion_component = proportion_component,
      return_posterior_sd = FALSE
    )

    pilot_gene_names <- dimnames(pilot_summary$ctsGEP)[["gene"]]
    precision_choice <- .cdseq3_choose_depth_adaptive_precision(
      rna_proportions = pilot_summary$rna_proportions,
      bulk_data = bulk_data[pilot_gene_names, , drop = FALSE],
      target_data_weight = adaptive_target_data_weight,
      depth_quantile = adaptive_depth_quantile,
      precision_bounds = adaptive_precision_bounds
    )
    selected_precision <- precision_choice$selected_precision

    precision_diagnostics <- c(
      list(
        strategy = precision_strategy,
        requested_reference_precision = as.numeric(reference_precision),
        selected_reference_precision = selected_precision,
        target_data_weight = adaptive_target_data_weight,
        depth_quantile = adaptive_depth_quantile,
        pilot_mcmc_iterations = as.integer(pilot_mcmc_iterations),
        pilot_save_last_n_mc = as.integer(pilot_save_last_n_mc),
        pilot_seed = as.integer(pilot_seed),
        pilot_elapsed_seconds = unname(pilot_timing[["elapsed"]])
      ),
      precision_choice[setdiff(names(precision_choice), "selected_precision")]
    )

    if (isTRUE(verbose)) {
      message(
        "CDSeq3 selected reference precision ",
        format(selected_precision, scientific = TRUE), "."
      )
    }

    pilot_summary <- NULL
    if (!isTRUE(keep_pilot_fit)) {
      pilot_fit <- NULL
    }
    invisible(gc())
  }

  final_args <- common_args
  final_args$reference_precision <- selected_precision

  if (isTRUE(verbose)) {
    message(
      "CDSeq3 final fit at reference precision ",
      format(selected_precision, scientific = TRUE), "."
    )
  }
  final_timing <- system.time({
    final_fit <- .cdseq3_call_cdseq2(final_args)
  })

  standardized <- .cdseq3_summarize_fit(
    fit = final_fit,
    bulk_data = bulk_data,
    reference_scRNA_annotation = reference_scRNA_annotation,
    proportion_component = proportion_component,
    return_posterior_sd = return_posterior_sd
  )

  standardized$selected_reference_precision <- selected_precision
  standardized$precision_strategy <- precision_strategy
  standardized$precision_diagnostics <- precision_diagnostics
  standardized$normalization_estimator <-
    "posterior mean of within-draw normalized ctsGEP"

  result <- list(
    call = match.call(),
    method = "cdseq3",
    base_method = "CDSeq2::cdseq2",
    precision_strategy = precision_strategy,
    requested_reference_precision = as.numeric(reference_precision),
    selected_reference_precision = selected_precision,
    precision_diagnostics = precision_diagnostics,
    ctsGEP = standardized$ctsGEP,
    ctsGEP_normalized = standardized$ctsGEP_normalized,
    ctsGEP_normalized_posthoc = standardized$ctsGEP_normalized_posthoc,
    ctsGEP_normalized_sd = standardized$ctsGEP_normalized_sd,
    profile_total_posterior_mean = standardized$profile_total_posterior_mean,
    profile_total_posterior_sd = standardized$profile_total_posterior_sd,
    rna_proportions = standardized$rna_proportions,
    proportion_source = standardized$proportion_source,
    proportion_audit = standardized$proportion_audit,
    normalization_audit = standardized$normalization_audit,
    normalized_gep_source = standardized$normalized_gep_source,
    posterior_means = standardized,
    elapsed_seconds = list(
      pilot = precision_diagnostics$pilot_elapsed_seconds,
      final = unname(final_timing[["elapsed"]]),
      total = precision_diagnostics$pilot_elapsed_seconds +
        unname(final_timing[["elapsed"]])
    ),
    cdseq2_fit = if (isTRUE(keep_cdseq2_fit)) final_fit else NULL,
    pilot_fit = if (isTRUE(keep_pilot_fit)) pilot_fit else NULL
  )

  class(result) <- c("cdseq3_fit", "list")
  result
}

#' @export
print.cdseq3_fit <- function(x, ...) {
  cat("CDSeq3 fit\n")
  cat("  Base method: ", x$base_method, "\n", sep = "")
  cat("  Precision strategy: ", x$precision_strategy, "\n", sep = "")
  cat(
    "  Selected reference precision: ",
    format(x$selected_reference_precision, scientific = TRUE),
    "\n",
    sep = ""
  )
  cat(
    "  Dimensions: ", paste(dim(x$ctsGEP), collapse = " x "),
    " (gene x cell type x sample)\n",
    sep = ""
  )
  invisible(x)
}

.cdseq3_validate_inputs <- function(
    bulk_data,
    reference_scRNA_count,
    reference_scRNA_annotation,
    reference_scRNA_donor,
    reference_precision,
    mcmc_iterations,
    save_last_n_mc,
    adaptive_target_data_weight,
    adaptive_depth_quantile,
    adaptive_precision_bounds,
    pilot_mcmc_iterations,
    pilot_save_last_n_mc,
    proportion_component
) {
  if (length(dim(bulk_data)) != 2L) {
    stop("bulk_data must be a two-dimensional gene-by-sample matrix.", call. = FALSE)
  }
  if (is.null(rownames(bulk_data)) || is.null(colnames(bulk_data))) {
    stop("bulk_data must have gene row names and sample column names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(bulk_data)) || anyDuplicated(colnames(bulk_data))) {
    stop("bulk_data gene and sample names must be unique.", call. = FALSE)
  }
  if (is.null(reference_scRNA_count) || is.null(reference_scRNA_annotation)) {
    stop(
      "reference_scRNA_count and reference_scRNA_annotation are required.",
      call. = FALSE
    )
  }
  if (length(reference_scRNA_annotation) == 0L) {
    stop("reference_scRNA_annotation is empty.", call. = FALSE)
  }
  if (
    !is.null(reference_scRNA_donor) &&
    length(reference_scRNA_donor) != length(reference_scRNA_annotation)
  ) {
    stop(
      "reference_scRNA_donor must have one value per reference annotation.",
      call. = FALSE
    )
  }
  if (
    length(reference_precision) != 1L ||
    !is.finite(reference_precision) || reference_precision <= 0
  ) {
    stop("reference_precision must be one positive finite number.", call. = FALSE)
  }
  if (
    length(mcmc_iterations) != 1L || length(save_last_n_mc) != 1L ||
    mcmc_iterations < 1L || save_last_n_mc < 1L ||
    save_last_n_mc > mcmc_iterations
  ) {
    stop("Require mcmc_iterations >= save_last_n_mc >= 1.", call. = FALSE)
  }
  if (
    length(pilot_mcmc_iterations) != 1L ||
    length(pilot_save_last_n_mc) != 1L ||
    pilot_mcmc_iterations < 1L || pilot_save_last_n_mc < 1L ||
    pilot_save_last_n_mc > pilot_mcmc_iterations
  ) {
    stop(
      "Require pilot_mcmc_iterations >= pilot_save_last_n_mc >= 1.",
      call. = FALSE
    )
  }
  if (
    length(adaptive_target_data_weight) != 1L ||
    !is.finite(adaptive_target_data_weight) ||
    adaptive_target_data_weight <= 0 || adaptive_target_data_weight >= 1
  ) {
    stop("adaptive_target_data_weight must be strictly between 0 and 1.", call. = FALSE)
  }
  if (
    length(adaptive_depth_quantile) != 1L ||
    !is.finite(adaptive_depth_quantile) ||
    adaptive_depth_quantile < 0 || adaptive_depth_quantile > 1
  ) {
    stop("adaptive_depth_quantile must be between 0 and 1.", call. = FALSE)
  }
  if (
    length(adaptive_precision_bounds) != 2L ||
    any(!is.finite(adaptive_precision_bounds)) ||
    any(adaptive_precision_bounds <= 0) ||
    adaptive_precision_bounds[[1L]] > adaptive_precision_bounds[[2L]]
  ) {
    stop(
      "adaptive_precision_bounds must contain positive lower and upper bounds.",
      call. = FALSE
    )
  }
  if (
    length(proportion_component) != 1L || is.na(proportion_component) ||
    !nzchar(proportion_component)
  ) {
    stop("proportion_component must be one nonempty component name.", call. = FALSE)
  }
  invisible(TRUE)
}

.cdseq3_call_cdseq2 <- function(args) {
  if (!requireNamespace("CDSeq2", quietly = TRUE)) {
    stop("The CDSeq2 package is required by cdseq3().", call. = FALSE)
  }
  if (!"cdseq2" %in% getNamespaceExports("CDSeq2")) {
    stop(
      "The installed CDSeq2 package does not export cdseq2(). ",
      "Install the revised package version before using cdseq3().",
      call. = FALSE
    )
  }

  fn <- getExportedValue("CDSeq2", "cdseq2")
  supported <- names(formals(fn))
  unsupported <- setdiff(names(args), supported)
  unsupported_nonnull <- unsupported[
    !vapply(args[unsupported], is.null, logical(1))
  ]
  if (length(unsupported_nonnull) > 0L) {
    stop(
      "The installed CDSeq2::cdseq2() does not support: ",
      paste(unsupported_nonnull, collapse = ", "),
      call. = FALSE
    )
  }
  args <- args[intersect(names(args), supported)]
  do.call(fn, args)
}

.cdseq3_standardize_proportion_draws <- function(x, expected_cell_types) {
  if (is.null(x) || !is.numeric(x)) return(NULL)

  x <- drop(x)
  vector_names <- names(x)
  n_cell_types <- length(expected_cell_types)

  if (is.null(dim(x))) {
    if (length(x) != n_cell_types) return(NULL)
    x <- matrix(x, nrow = n_cell_types, ncol = 1L)
    rownames(x) <- vector_names
  } else if (length(dim(x)) == 2L) {
    if (nrow(x) == n_cell_types) {
      x <- as.matrix(x)
    } else if (ncol(x) == n_cell_types) {
      x <- t(x)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }

  if (!is.null(rownames(x)) && all(rownames(x) != "")) {
    if (!setequal(rownames(x), expected_cell_types)) return(NULL)
    x <- x[expected_cell_types, , drop = FALSE]
  } else {
    rownames(x) <- expected_cell_types
  }

  if (anyNA(x) || any(!is.finite(x)) || any(x < 0)) return(NULL)
  draw_sums <- colSums(x)
  if (any(!is.finite(draw_sums)) || any(draw_sums <= 0)) return(NULL)
  if (any(abs(draw_sums - 1) > 1e-6)) return(NULL)
  sweep(x, 2L, draw_sums, "/")
}

.cdseq3_find_proportion_draws <- function(
    sample_result,
    expected_cell_types,
    preferred_component = "estProp"
) {
  if (!is.list(sample_result) || is.null(names(sample_result))) return(NULL)

  known_names <- c(
    preferred_component,
    "estProp",
    "estProp_RB",
    "estProp_legacy",
    "theta_array",
    "estProp_array",
    "estProportion_array",
    "estProportions_array",
    "proportion_array",
    "theta",
    "estProportion",
    "cell_type_proportions"
  )
  candidate_names <- unique(c(
    intersect(known_names, names(sample_result)),
    grep(
      "prop|proportion|theta",
      names(sample_result),
      value = TRUE,
      ignore.case = TRUE
    )
  ))

  valid <- lapply(candidate_names, function(component_name) {
    .cdseq3_standardize_proportion_draws(
      sample_result[[component_name]],
      expected_cell_types
    )
  })
  names(valid) <- candidate_names
  valid <- valid[!vapply(valid, is.null, logical(1))]
  if (length(valid) == 0L) return(NULL)

  if (preferred_component %in% names(valid)) {
    selected_name <- preferred_component
  } else if ("estProp" %in% names(valid)) {
    selected_name <- "estProp"
  } else if (length(valid) == 1L) {
    selected_name <- names(valid)[[1L]]
  } else {
    stop(
      "The requested proportion component '", preferred_component,
      "' is unavailable and multiple valid alternatives remain: ",
      paste(names(valid), collapse = ", "),
      call. = FALSE
    )
  }

  selected <- valid[[selected_name]]
  selected_mean <- rowMeans(selected)
  comparison <- do.call(
    rbind,
    lapply(names(valid), function(component_name) {
      data.frame(
        component = component_name,
        selected = identical(component_name, selected_name),
        retained_draws = ncol(valid[[component_name]]),
        max_abs_mean_difference_from_selected = max(abs(
          rowMeans(valid[[component_name]]) - selected_mean
        )),
        stringsAsFactors = FALSE
      )
    })
  )

  list(
    draws = selected,
    source = selected_name,
    available_components = names(valid),
    comparison = comparison
  )
}

.cdseq3_infer_output_axes <- function(
    gep_array,
    sample_name,
    bulk_genes,
    reference_cell_types
) {
  if (is.null(gep_array) || length(dim(gep_array)) != 3L) {
    stop("Invalid estGEP_array for sample: ", sample_name, call. = FALSE)
  }
  dn <- dimnames(gep_array)
  if (is.null(dn) || is.null(dn[[1L]]) || is.null(dn[[2L]])) {
    stop(
      "estGEP_array lacks gene or cell-type names for sample: ",
      sample_name,
      call. = FALSE
    )
  }

  first_gene_score <- mean(dn[[1L]] %in% bulk_genes)
  second_gene_score <- mean(dn[[2L]] %in% bulk_genes)
  first_type_score <- mean(dn[[1L]] %in% reference_cell_types)
  second_type_score <- mean(dn[[2L]] %in% reference_cell_types)

  if (
    first_gene_score == 1 &&
    (second_type_score == 1 || second_gene_score < first_gene_score)
  ) {
    gene_axis <- 1L
    cell_type_axis <- 2L
  } else if (
    second_gene_score == 1 &&
    (first_type_score == 1 || first_gene_score < second_gene_score)
  ) {
    gene_axis <- 2L
    cell_type_axis <- 1L
  } else if (
    dim(gep_array)[1L] > dim(gep_array)[2L] &&
    first_gene_score >= second_gene_score
  ) {
    gene_axis <- 1L
    cell_type_axis <- 2L
  } else if (
    dim(gep_array)[2L] > dim(gep_array)[1L] &&
    second_gene_score >= first_gene_score
  ) {
    gene_axis <- 2L
    cell_type_axis <- 1L
  } else {
    stop(
      "Could not infer gene and cell-type axes for sample: ",
      sample_name,
      call. = FALSE
    )
  }

  output_genes <- dn[[gene_axis]]
  output_cell_types <- dn[[cell_type_axis]]
  if (!all(output_genes %in% bulk_genes)) {
    stop("The fitted GEP contains genes absent from bulk_data.", call. = FALSE)
  }
  output_genes <- bulk_genes[bulk_genes %in% output_genes]
  if (length(output_genes) == 0L || length(output_cell_types) == 0L) {
    stop("The fitted GEP has an empty gene or cell-type axis.", call. = FALSE)
  }

  list(
    genes = output_genes,
    cell_types = output_cell_types
  )
}

.cdseq3_orient_gep_draws <- function(
    gep_array,
    sample_name,
    expected_genes,
    expected_cell_types
) {
  if (is.null(gep_array) || length(dim(gep_array)) != 3L) {
    stop("Invalid estGEP_array for sample: ", sample_name, call. = FALSE)
  }
  dn <- dimnames(gep_array)
  if (is.null(dn) || is.null(dn[[1L]]) || is.null(dn[[2L]])) {
    stop(
      "estGEP_array lacks gene or cell-type names for sample: ",
      sample_name,
      call. = FALSE
    )
  }

  first_is_gene <- setequal(dn[[1L]], expected_genes)
  second_is_gene <- setequal(dn[[2L]], expected_genes)
  first_is_cell_type <- setequal(dn[[1L]], expected_cell_types)
  second_is_cell_type <- setequal(dn[[2L]], expected_cell_types)

  if (first_is_gene && second_is_cell_type) {
    oriented <- gep_array
  } else if (first_is_cell_type && second_is_gene) {
    oriented <- aperm(gep_array, c(2L, 1L, 3L))
  } else if (
    dim(gep_array)[1L] == length(expected_genes) &&
    dim(gep_array)[2L] == length(expected_cell_types)
  ) {
    oriented <- gep_array
  } else if (
    dim(gep_array)[2L] == length(expected_genes) &&
    dim(gep_array)[1L] == length(expected_cell_types)
  ) {
    oriented <- aperm(gep_array, c(2L, 1L, 3L))
  } else {
    stop(
      "Could not determine estGEP_array orientation for sample: ",
      sample_name,
      call. = FALSE
    )
  }

  if (!setequal(dimnames(oriented)[[1L]], expected_genes)) {
    stop("The fitted gene set does not match bulk_data.", call. = FALSE)
  }
  if (!setequal(dimnames(oriented)[[2L]], expected_cell_types)) {
    stop("The fitted cell-type set does not match the reference labels.", call. = FALSE)
  }
  oriented[expected_genes, expected_cell_types, , drop = FALSE]
}

.cdseq3_summarize_fit <- function(
    fit,
    bulk_data,
    reference_scRNA_annotation,
    proportion_component,
    return_posterior_sd = TRUE
) {
  result_list <- fit$cdseq_result
  if (is.null(result_list) || !is.list(result_list)) {
    stop("The cdseq2 fit does not contain a cdseq_result list.", call. = FALSE)
  }
  if (is.null(names(result_list)) || any(names(result_list) == "")) {
    stop("The cdseq_result list must be named by sample.", call. = FALSE)
  }

  sample_names <- colnames(bulk_data)
  bulk_gene_names <- rownames(bulk_data)
  reference_cell_types <- unique(as.character(reference_scRNA_annotation))
  if (!setequal(names(result_list), sample_names)) {
    stop("cdseq2 output samples do not match bulk_data.", call. = FALSE)
  }
  result_list <- result_list[sample_names]

  output_axes <- .cdseq3_infer_output_axes(
    gep_array = result_list[[sample_names[[1L]]]]$estGEP_array,
    sample_name = sample_names[[1L]],
    bulk_genes = bulk_gene_names,
    reference_cell_types = reference_cell_types
  )
  gene_names <- output_axes$genes
  cell_types <- output_axes$cell_types

  n_genes <- length(gene_names)
  n_cell_types <- length(cell_types)
  n_samples <- length(sample_names)

  count_mean_array <- array(
    NA_real_,
    dim = c(n_genes, n_cell_types, n_samples),
    dimnames = list(
      gene = gene_names,
      cell_type = cell_types,
      sample = sample_names
    )
  )
  phi_mean_array <- count_mean_array
  phi_posthoc_array <- count_mean_array
  phi_sd_array <- if (isTRUE(return_posterior_sd)) count_mean_array else NULL
  total_mean_matrix <- matrix(
    NA_real_,
    nrow = n_samples,
    ncol = n_cell_types,
    dimnames = list(sample_names, cell_types)
  )
  total_sd_matrix <- total_mean_matrix
  proportion_matrix <- total_mean_matrix
  proportion_source <- setNames(character(n_samples), sample_names)
  audit_rows <- vector("list", n_samples)

  for (sample_index in seq_along(sample_names)) {
    sample_name <- sample_names[[sample_index]]
    sample_result <- result_list[[sample_name]]
    gep_draws <- .cdseq3_orient_gep_draws(
      gep_array = sample_result$estGEP_array,
      sample_name = sample_name,
      expected_genes = gene_names,
      expected_cell_types = cell_types
    )

    if (anyNA(gep_draws) || any(!is.finite(gep_draws)) || any(gep_draws < 0)) {
      stop(
        "Nonfinite or negative values in estGEP_array for sample: ",
        sample_name,
        call. = FALSE
      )
    }

    count_mean <- apply(gep_draws, c(1L, 2L), mean)
    profile_totals <- apply(gep_draws, c(2L, 3L), sum)
    profile_totals <- matrix(
      profile_totals,
      nrow = n_cell_types,
      ncol = dim(gep_draws)[3L],
      dimnames = list(cell_types, dimnames(gep_draws)[[3L]])
    )
    positive_total <- is.finite(profile_totals) & profile_totals > 0
    denominator <- profile_totals
    denominator[!positive_total] <- NA_real_
    phi_draws <- sweep(gep_draws, c(2L, 3L), denominator, "/")

    phi_mean <- apply(phi_draws, c(1L, 2L), function(z) {
      z <- z[is.finite(z)]
      if (length(z) == 0L) NA_real_ else mean(z)
    })
    if (anyNA(phi_mean) || any(!is.finite(phi_mean))) {
      stop(
        "Direct posterior normalized GEP is undefined for sample: ",
        sample_name,
        call. = FALSE
      )
    }

    profile_mean_totals <- colSums(count_mean)
    if (any(!is.finite(profile_mean_totals)) || any(profile_mean_totals <= 0)) {
      stop("Posterior-mean profile total is nonpositive for ", sample_name)
    }
    phi_posthoc <- sweep(count_mean, 2L, profile_mean_totals, "/")

    if (isTRUE(return_posterior_sd)) {
      phi_sd_array[, , sample_index] <- apply(
        phi_draws,
        c(1L, 2L),
        function(z) {
          z <- z[is.finite(z)]
          if (length(z) < 2L) NA_real_ else stats::sd(z)
        }
      )
    }

    count_mean_array[, , sample_index] <- count_mean
    phi_mean_array[, , sample_index] <- phi_mean
    phi_posthoc_array[, , sample_index] <- phi_posthoc
    total_mean_matrix[sample_name, ] <- rowMeans(profile_totals, na.rm = TRUE)
    total_sd_matrix[sample_name, ] <- apply(
      profile_totals,
      1L,
      stats::sd,
      na.rm = TRUE
    )

    prop_result <- .cdseq3_find_proportion_draws(
      sample_result = sample_result,
      expected_cell_types = cell_types,
      preferred_component = proportion_component
    )
    if (is.null(prop_result)) {
      stop(
        "Could not identify posterior RNA proportions for sample ",
        sample_name,
        ". Available components: ",
        paste(names(sample_result), collapse = ", "),
        call. = FALSE
      )
    }
    proportion_matrix[sample_name, ] <- rowMeans(prop_result$draws)
    proportion_source[[sample_name]] <- prop_result$source

    estprop2_difference <- prop_result$comparison[
      prop_result$comparison$component == "estProp2",
      "max_abs_mean_difference_from_selected"
    ]
    if (length(estprop2_difference) == 0L) estprop2_difference <- NA_real_

    allocation_fraction <- total_mean_matrix[sample_name, ] /
      sum(total_mean_matrix[sample_name, ])

    audit_rows[[sample_index]] <- data.frame(
      sample = sample_name,
      selected_component = prop_result$source,
      available_valid_components = paste(
        prop_result$available_components,
        collapse = "; "
      ),
      selected_retained_draws = ncol(prop_result$draws),
      max_abs_difference_from_estProp2 = estprop2_difference[[1L]],
      max_abs_difference_from_allocation_fraction = max(abs(
        proportion_matrix[sample_name, ] - allocation_fraction
      )),
      direct_phi_max_sum_error = max(abs(colSums(phi_mean) - 1)),
      direct_vs_posthoc_max_abs_difference = max(abs(
        phi_mean - phi_posthoc
      )),
      zero_total_celltype_draws = sum(!positive_total),
      stringsAsFactors = FALSE
    )
  }

  if (any(abs(rowSums(proportion_matrix) - 1) > 1e-8)) {
    stop("Posterior-mean RNA proportions do not sum to one.", call. = FALSE)
  }
  direct_sum_error <- max(abs(apply(phi_mean_array, c(2L, 3L), sum) - 1))
  posthoc_sum_error <- max(abs(apply(phi_posthoc_array, c(2L, 3L), sum) - 1))
  if (direct_sum_error > 1e-8 || posthoc_sum_error > 1e-8) {
    stop("A normalized ctsGEP output does not sum to one.", call. = FALSE)
  }

  normalization_audit <- data.frame(
    max_direct_profile_sum_error = direct_sum_error,
    max_posthoc_profile_sum_error = posthoc_sum_error,
    max_direct_vs_posthoc_abs_difference = max(abs(
      phi_mean_array - phi_posthoc_array
    )),
    median_direct_vs_posthoc_abs_difference = stats::median(abs(
      phi_mean_array - phi_posthoc_array
    )),
    stringsAsFactors = FALSE
  )

  list(
    ctsGEP = count_mean_array,
    ctsGEP_normalized = phi_mean_array,
    ctsGEP_normalized_posthoc = phi_posthoc_array,
    ctsGEP_normalized_sd = phi_sd_array,
    profile_total_posterior_mean = total_mean_matrix,
    profile_total_posterior_sd = total_sd_matrix,
    rna_proportions = proportion_matrix,
    proportion_source = proportion_source,
    proportion_audit = do.call(rbind, audit_rows),
    normalization_audit = normalization_audit,
    extraction_contract = "cdseq3_direct_posterior_phi",
    normalized_gep_source = paste(
      "mean over retained draws after within-draw",
      "cell-type profile normalization"
    )
  )
}

.cdseq3_choose_depth_adaptive_precision <- function(
    rna_proportions,
    bulk_data,
    target_data_weight,
    depth_quantile,
    precision_bounds
) {
  bulk_depth <- colSums(bulk_data)
  if (!setequal(rownames(rna_proportions), names(bulk_depth))) {
    stop("RNA-proportion samples do not match bulk_data.", call. = FALSE)
  }
  rna_proportions <- rna_proportions[names(bulk_depth), , drop = FALSE]

  estimated_depth <- sweep(rna_proportions, 1L, bulk_depth, "*")
  median_depth <- apply(estimated_depth, 2L, stats::median, na.rm = TRUE)
  median_depth <- median_depth[is.finite(median_depth) & median_depth > 0]
  if (length(median_depth) == 0L) {
    stop("Could not estimate positive cell-type RNA depths from pilot fit.", call. = FALSE)
  }

  depth_scale <- unname(stats::quantile(
    median_depth,
    probs = depth_quantile,
    na.rm = TRUE,
    names = FALSE,
    type = 8
  ))
  unclipped <- depth_scale * (1 - target_data_weight) / target_data_weight
  selected <- min(
    max(unclipped, precision_bounds[[1L]]),
    precision_bounds[[2L]]
  )
  selected <- signif(selected, 3L)

  achieved_weight <- median_depth / (selected + median_depth)
  list(
    selected_precision = selected,
    median_estimated_cell_type_depth = median_depth,
    estimated_data_weight_by_cell_type = achieved_weight,
    depth_scale_used = depth_scale,
    unclipped_precision = unclipped,
    precision_bounds = as.numeric(precision_bounds)
  )
}
