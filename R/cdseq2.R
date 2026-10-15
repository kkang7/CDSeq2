# Internal utility: validate a raw-count matrix.
.cdseq2_validate_count_matrix <- function(x,
                                          name,
                                          require_colnames = FALSE) {
  if (is.null(x)) {
    stop(name, " is NULL.", call. = FALSE)
  }
  if (!is.matrix(x)) {
    stop(name, " must be a base R matrix.", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop(name, " must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " contains NA, NaN, or infinite values.", call. = FALSE)
  }
  if (any(x < 0)) {
    stop(name, " contains negative values.", call. = FALSE)
  }
  if (any(x != round(x))) {
    stop(name, " must contain unnormalized integer counts.", call. = FALSE)
  }
  if (is.null(rownames(x))) {
    stop(name, " must have gene names as row names.", call. = FALSE)
  }
  if (anyNA(rownames(x)) || any(!nzchar(rownames(x)))) {
    stop(name, " has missing or empty gene names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(x))) {
    stop(name, " must have unique gene names.", call. = FALSE)
  }
  if (require_colnames) {
    if (is.null(colnames(x))) {
      stop(name, " must have column names.", call. = FALSE)
    }
    if (anyNA(colnames(x)) || any(!nzchar(colnames(x)))) {
      stop(name, " has missing or empty column names.", call. = FALSE)
    }
    if (anyDuplicated(colnames(x))) {
      stop(name, " must have unique column names.", call. = FALSE)
    }
  }
  invisible(TRUE)
}


# Internal utility: expand a scalar, an ordered vector, or a named vector to
# one value per modeled cell type/state.
.cdseq2_expand_parameter <- function(x,
                                     cell_types,
                                     name,
                                     strictly_positive = TRUE,
                                     upper_bound = NULL) {
  if (is.null(x) || !is.numeric(x) || !is.atomic(x) || is.matrix(x)) {
    stop(name, " must be a numeric scalar or vector.", call. = FALSE)
  }

  if (length(x) == 1L) {
    ans <- rep(as.numeric(x), length(cell_types))
  } else if (!is.null(names(x))) {
    if (anyDuplicated(names(x))) {
      stop(name, " has duplicated names.", call. = FALSE)
    }
    missing_types <- setdiff(cell_types, names(x))
    if (length(missing_types) > 0L) {
      stop(
        name, " is missing values for: ",
        paste(missing_types, collapse = ", "),
        call. = FALSE
      )
    }
    ans <- as.numeric(x[cell_types])
  } else {
    if (length(x) != length(cell_types)) {
      stop(
        name, " must be scalar or have length equal to the number of ",
        "modeled cell types/states (", length(cell_types), ").",
        call. = FALSE
      )
    }
    ans <- as.numeric(x)
  }

  if (any(!is.finite(ans))) {
    stop(name, " contains a non-finite value.", call. = FALSE)
  }
  if (strictly_positive && any(ans <= 0)) {
    stop(name, " must contain only positive values.", call. = FALSE)
  }
  if (!strictly_positive && any(ans < 0)) {
    stop(name, " must contain only nonnegative values.", call. = FALSE)
  }
  if (!is.null(upper_bound) && any(ans > upper_bound)) {
    stop(name, " cannot exceed ", upper_bound, ".", call. = FALSE)
  }

  names(ans) <- cell_types
  ans
}


.cdseq2_row_normalize <- function(x, object_name = "matrix") {
  totals <- rowSums(x)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop(object_name, " contains a row with nonpositive total mass.", call. = FALSE)
  }
  sweep(x, 1L, totals, "/")
}


# Construct the T x G reference profile.  The corrected default adds the
# smoothing pseudocount after aggregation, rather than once to every cell.
.cdseq2_build_reference <- function(reference_scRNA_count,
                                    annotation,
                                    donor,
                                    cell_types,
                                    aggregation,
                                    min_cells_per_donor_type,
                                    smooth_ref,
                                    smoothing_alpha,
                                    smoothing_scope) {
  gene_symbols <- rownames(reference_scRNA_count)
  annotation <- as.character(annotation)
  cell_count_by_type <- as.integer(table(factor(annotation, levels = cell_types)))
  names(cell_count_by_type) <- cell_types

  donor_group_table <- NULL

  if (aggregation == "pooled_counts") {
    aggregated <- rowsum(
      t(reference_scRNA_count),
      group = annotation,
      reorder = TRUE
    )
    aggregated <- as.matrix(aggregated[cell_types, , drop = FALSE])
    rownames(aggregated) <- cell_types
    colnames(aggregated) <- gene_symbols

    if (smooth_ref) {
      multiplier <- if (smoothing_scope == "per_cell") {
        cell_count_by_type
      } else {
        rep(1, length(cell_types))
      }
      pseudo <- smoothing_alpha * multiplier
      aggregated <- sweep(aggregated, 1L, pseudo, "+")
    }

    reference_profile <- .cdseq2_row_normalize(
      aggregated,
      object_name = "aggregated reference"
    )
    donor_count_by_type <- rep(NA_integer_, length(cell_types))
    names(donor_count_by_type) <- cell_types

  } else {
    if (is.null(donor)) {
      stop(
        "reference_scRNA_donor is required when ",
        "reference_aggregation = 'equal_donor'.",
        call. = FALSE
      )
    }

    donor <- as.character(donor)
    group_key <- paste(annotation, donor, sep = "\u001f")
    unique_keys <- unique(group_key)
    group_id <- match(group_key, unique_keys)

    donor_pseudobulk <- rowsum(
      t(reference_scRNA_count),
      group = group_id,
      reorder = FALSE
    )
    donor_pseudobulk <- as.matrix(donor_pseudobulk)
    colnames(donor_pseudobulk) <- gene_symbols

    first_index <- match(seq_along(unique_keys), group_id)
    group_annotation <- annotation[first_index]
    group_donor <- donor[first_index]
    group_cell_count <- tabulate(group_id, nbins = length(unique_keys))

    keep_group <- group_cell_count >= min_cells_per_donor_type
    if (any(!keep_group)) {
      donor_pseudobulk <- donor_pseudobulk[keep_group, , drop = FALSE]
      group_annotation <- group_annotation[keep_group]
      group_donor <- group_donor[keep_group]
      group_cell_count <- group_cell_count[keep_group]
    }

    donor_count_by_type <- as.integer(
      table(factor(group_annotation, levels = cell_types))
    )
    names(donor_count_by_type) <- cell_types

    if (any(donor_count_by_type == 0L)) {
      missing_types <- names(donor_count_by_type)[donor_count_by_type == 0L]
      stop(
        "No eligible donor-specific reference pseudobulk remains for: ",
        paste(missing_types, collapse = ", "),
        call. = FALSE
      )
    }

    if (smooth_ref) {
      alpha_by_group <- smoothing_alpha[
        match(group_annotation, cell_types)
      ]
      multiplier <- if (smoothing_scope == "per_cell") {
        group_cell_count
      } else {
        rep(1, length(group_cell_count))
      }
      donor_pseudobulk <- sweep(
        donor_pseudobulk,
        1L,
        alpha_by_group * multiplier,
        "+"
      )
    }

    donor_profiles <- .cdseq2_row_normalize(
      donor_pseudobulk,
      object_name = "donor-by-cell-type reference pseudobulk"
    )

    summed_profiles <- rowsum(
      donor_profiles,
      group = group_annotation,
      reorder = TRUE
    )
    summed_profiles <- as.matrix(
      summed_profiles[cell_types, , drop = FALSE]
    )
    reference_profile <- sweep(
      summed_profiles,
      1L,
      donor_count_by_type,
      "/"
    )
    reference_profile <- .cdseq2_row_normalize(
      reference_profile,
      object_name = "equal-donor reference profile"
    )
    rownames(reference_profile) <- cell_types
    colnames(reference_profile) <- gene_symbols

    donor_group_table <- data.frame(
      cell_type = group_annotation,
      donor = group_donor,
      n_cells = group_cell_count,
      stringsAsFactors = FALSE
    )
  }

  if (any(reference_profile <= 0)) {
    stop(
      "The reference profile contains zero probabilities. Dirichlet ",
      "parameters must be strictly positive. Use smooth_ref = TRUE with ",
      "positive smoothing_alpha, or remove unsupported genes.",
      call. = FALSE
    )
  }

  list(
    profile = reference_profile,
    cell_count_by_type = cell_count_by_type,
    donor_count_by_type = donor_count_by_type,
    donor_group_table = donor_group_table
  )
}


.cdseq2_prepare_aggregation_map <- function(aggregate_to, cell_types) {
  if (is.null(aggregate_to)) {
    return(NULL)
  }

  if (!is.atomic(aggregate_to)) {
    stop("aggregate_to must be a character vector.", call. = FALSE)
  }
  aggregate_names <- names(aggregate_to)
  aggregate_to <- as.character(aggregate_to)
  names(aggregate_to) <- aggregate_names

  if (!is.null(names(aggregate_to))) {
    if (anyDuplicated(names(aggregate_to))) {
      stop("aggregate_to has duplicated names.", call. = FALSE)
    }
    missing_states <- setdiff(cell_types, names(aggregate_to))
    if (length(missing_states) > 0L) {
      stop(
        "aggregate_to is missing modeled states: ",
        paste(missing_states, collapse = ", "),
        call. = FALSE
      )
    }
    aggregate_to <- aggregate_to[cell_types]
  } else {
    if (length(aggregate_to) != length(cell_types)) {
      stop(
        "An unnamed aggregate_to vector must have one entry per modeled ",
        "cell type/state.",
        call. = FALSE
      )
    }
    names(aggregate_to) <- cell_types
  }

  if (anyNA(aggregate_to) || any(!nzchar(aggregate_to))) {
    stop("aggregate_to contains missing or empty major-type labels.", call. = FALSE)
  }

  aggregate_to
}


.cdseq2_combine_3d_arrays <- function(arrays, third_dim_names = NULL) {
  if (length(arrays) == 0L) {
    stop("No arrays were supplied for combination.", call. = FALSE)
  }

  first_dim <- dim(arrays[[1L]])
  if (length(first_dim) != 3L) {
    stop("All objects must be three-dimensional arrays.", call. = FALSE)
  }

  for (a in arrays) {
    if (length(dim(a)) != 3L || any(dim(a)[1:2] != first_dim[1:2])) {
      stop("MCMC arrays have incompatible dimensions.", call. = FALSE)
    }
  }

  total_draws <- sum(vapply(arrays, function(a) dim(a)[3L], integer(1L)))
  ans <- array(
    0,
    dim = c(first_dim[1L], first_dim[2L], total_draws),
    dimnames = list(
      dimnames(arrays[[1L]])[[1L]],
      dimnames(arrays[[1L]])[[2L]],
      NULL
    )
  )

  start <- 1L
  for (a in arrays) {
    end <- start + dim(a)[3L] - 1L
    ans[, , start:end] <- a
    start <- end + 1L
  }

  if (!is.null(third_dim_names)) {
    if (length(third_dim_names) != total_draws) {
      stop("third_dim_names has the wrong length.", call. = FALSE)
    }
    dimnames(ans)[[3L]] <- third_dim_names
  }

  ans
}


# Rao-Blackwellized posterior mean of phi from collapsed allocation draws.
.cdseq2_phi_mean_only <- function(z_array, beta) {
  n_draws <- dim(z_array)[3L]
  ans <- matrix(
    0,
    nrow = dim(z_array)[1L],
    ncol = dim(z_array)[2L],
    dimnames = dimnames(z_array)[1:2]
  )

  for (b in seq_len(n_draws)) {
    posterior_alpha <- z_array[, , b, drop = FALSE][, , 1L] + beta
    ans <- ans + sweep(
      posterior_alpha,
      1L,
      rowSums(posterior_alpha),
      "/"
    )
  }

  ans / n_draws
}


.cdseq2_basic_rhat <- function(draw_matrix) {
  draw_matrix <- as.matrix(draw_matrix)
  n_draws <- nrow(draw_matrix)
  n_chains <- ncol(draw_matrix)

  if (n_chains < 2L || n_draws < 2L) {
    return(NA_real_)
  }

  within_variances <- apply(draw_matrix, 2L, stats::var)
  within <- mean(within_variances)
  between <- n_draws * stats::var(colMeans(draw_matrix))

  if (!is.finite(within) || !is.finite(between)) {
    return(NA_real_)
  }
  if (within == 0 && between == 0) {
    return(1)
  }
  if (within == 0 && between > 0) {
    return(Inf)
  }

  variance_hat <- ((n_draws - 1) / n_draws) * within + between / n_draws
  sqrt(variance_hat / within)
}


.cdseq2_chain_diagnostics <- function(z_array,
                                      ssp_array,
                                      beta,
                                      alpha,
                                      chain_id,
                                      cell_types) {
  chain_levels <- unique(chain_id)
  n_chains <- length(chain_levels)

  empty_prop <- data.frame(
    cell_type = cell_types,
    rhat = NA_real_,
    chain_mean_range = NA_real_,
    stringsAsFactors = FALSE
  )
  empty_gep <- data.frame(
    cell_type = cell_types,
    mean_pairwise_correlation = NA_real_,
    min_pairwise_correlation = NA_real_,
    max_pairwise_total_variation = NA_real_,
    stringsAsFactors = FALSE
  )

  if (n_chains < 2L) {
    return(list(
      proportion = empty_prop,
      ctsGEP = empty_gep,
      draws_per_chain = as.integer(table(chain_id))
    ))
  }

  draws_per_chain <- as.integer(table(factor(chain_id, levels = chain_levels)))
  common_draws <- min(draws_per_chain)
  if (common_draws < 2L) {
    return(list(
      proportion = empty_prop,
      ctsGEP = empty_gep,
      draws_per_chain = draws_per_chain
    ))
  }

  n_types <- dim(ssp_array)[1L]
  n_total_draws <- dim(ssp_array)[3L]
  theta_conditional <- matrix(0, nrow = n_types, ncol = n_total_draws)
  for (b in seq_len(n_total_draws)) {
    a <- ssp_array[, 1L, b] + alpha
    theta_conditional[, b] <- a / sum(a)
  }

  rhat <- numeric(n_types)
  chain_range <- numeric(n_types)
  for (t in seq_len(n_types)) {
    chain_matrix <- vapply(
      chain_levels,
      function(ch) {
        idx <- which(chain_id == ch)[seq_len(common_draws)]
        theta_conditional[t, idx]
      },
      numeric(common_draws)
    )
    rhat[t] <- .cdseq2_basic_rhat(chain_matrix)
    chain_range[t] <- diff(range(colMeans(chain_matrix)))
  }

  phi_by_chain <- lapply(chain_levels, function(ch) {
    idx <- which(chain_id == ch)
    .cdseq2_phi_mean_only(z_array[, , idx, drop = FALSE], beta)
  })

  pair_index <- utils::combn(seq_len(n_chains), 2L)
  mean_cor <- min_cor <- max_tv <- rep(NA_real_, n_types)

  for (t in seq_len(n_types)) {
    correlations <- numeric(ncol(pair_index))
    total_variations <- numeric(ncol(pair_index))
    for (p in seq_len(ncol(pair_index))) {
      first <- phi_by_chain[[pair_index[1L, p]]][t, ]
      second <- phi_by_chain[[pair_index[2L, p]]][t, ]
      correlations[p] <- suppressWarnings(
        stats::cor(first, second, use = "pairwise.complete.obs")
      )
      total_variations[p] <- 0.5 * sum(abs(first - second))
    }
    mean_cor[t] <- mean(correlations, na.rm = TRUE)
    min_cor[t] <- min(correlations, na.rm = TRUE)
    max_tv[t] <- max(total_variations, na.rm = TRUE)
    if (!is.finite(mean_cor[t])) mean_cor[t] <- NA_real_
    if (!is.finite(min_cor[t])) min_cor[t] <- NA_real_
    if (!is.finite(max_tv[t])) max_tv[t] <- NA_real_
  }

  list(
    proportion = data.frame(
      cell_type = cell_types,
      rhat = rhat,
      chain_mean_range = chain_range,
      stringsAsFactors = FALSE
    ),
    ctsGEP = data.frame(
      cell_type = cell_types,
      mean_pairwise_correlation = mean_cor,
      min_pairwise_correlation = min_cor,
      max_pairwise_total_variation = max_tv,
      stringsAsFactors = FALSE
    ),
    draws_per_chain = draws_per_chain
  )
}


# Summarize allocation MCMC draws for one bulk sample.  The recommended ctsGEP
# is E[phi | X], computed as an average of E[phi | allocation, X].
.cdseq2_summarize_sample <- function(z_array,
                                     ssp_array,
                                     beta,
                                     alpha,
                                     cell_types,
                                     gene_symbols,
                                     return_posterior_sd,
                                     return_phi_draw_means,
                                     low_information_depth,
                                     low_information_prior_weight) {
  n_types <- dim(z_array)[1L]
  n_genes <- dim(z_array)[2L]
  n_draws <- dim(z_array)[3L]

  if (any(dim(ssp_array) != c(n_types, 1L, n_draws))) {
    stop("The saved SSP and ctsGEP arrays are dimensionally inconsistent.", call. = FALSE)
  }

  rho <- rowSums(beta)
  phi_sum <- matrix(0, nrow = n_types, ncol = n_genes)
  theta_sum <- numeric(n_types)
  z_sum <- matrix(0, nrow = n_types, ncol = n_genes)
  ssp_sum <- numeric(n_types)
  allocation_profile_sum <- matrix(0, nrow = n_types, ncol = n_genes)
  allocation_profile_n <- numeric(n_types)
  depth_draws <- matrix(0, nrow = n_types, ncol = n_draws)
  theta_conditional_draws <- matrix(0, nrow = n_types, ncol = n_draws)

  phi_second_sum <- if (return_posterior_sd) {
    matrix(0, nrow = n_types, ncol = n_genes)
  } else {
    NULL
  }
  theta_second_sum <- if (return_posterior_sd) numeric(n_types) else NULL

  phi_draw_means <- if (return_phi_draw_means) {
    array(
      0,
      dim = c(n_types, n_genes, n_draws),
      dimnames = dimnames(z_array)
    )
  } else {
    NULL
  }

  maximum_assignment_mismatch <- 0
  bulk_total_by_draw <- numeric(n_draws)

  for (b in seq_len(n_draws)) {
    z_b <- z_array[, , b, drop = FALSE][, , 1L]
    ssp_b <- as.numeric(ssp_array[, 1L, b])
    depth_b <- rowSums(z_b)

    depth_draws[, b] <- depth_b
    bulk_total_by_draw[b] <- sum(depth_b)
    maximum_assignment_mismatch <- max(
      maximum_assignment_mismatch,
      max(abs(depth_b - ssp_b))
    )

    posterior_alpha <- z_b + beta
    posterior_total <- rowSums(posterior_alpha)
    phi_b <- sweep(posterior_alpha, 1L, posterior_total, "/")

    theta_alpha <- ssp_b + alpha
    theta_total <- sum(theta_alpha)
    theta_b <- theta_alpha / theta_total

    phi_sum <- phi_sum + phi_b
    theta_sum <- theta_sum + theta_b
    z_sum <- z_sum + z_b
    ssp_sum <- ssp_sum + ssp_b
    theta_conditional_draws[, b] <- theta_b

    positive_depth <- which(depth_b > 0)
    if (length(positive_depth) > 0L) {
      allocation_profile_sum[positive_depth, ] <-
        allocation_profile_sum[positive_depth, , drop = FALSE] +
        sweep(
          z_b[positive_depth, , drop = FALSE],
          1L,
          depth_b[positive_depth],
          "/"
        )
      allocation_profile_n[positive_depth] <-
        allocation_profile_n[positive_depth] + 1
    }

    if (return_posterior_sd) {
      conditional_phi_variance <- sweep(
        phi_b * (1 - phi_b),
        1L,
        posterior_total + 1,
        "/"
      )
      phi_second_sum <- phi_second_sum + conditional_phi_variance + phi_b^2

      conditional_theta_variance <- theta_b * (1 - theta_b) /
        (theta_total + 1)
      theta_second_sum <- theta_second_sum +
        conditional_theta_variance + theta_b^2
    }

    if (return_phi_draw_means) {
      phi_draw_means[, , b] <- phi_b
    }
  }

  phi_mean <- phi_sum / n_draws
  theta_mean <- theta_sum / n_draws
  z_mean <- z_sum / n_draws
  ssp_mean <- ssp_sum / n_draws

  phi_sd <- theta_sd <- NULL
  if (return_posterior_sd) {
    phi_variance <- phi_second_sum / n_draws - phi_mean^2
    phi_variance[phi_variance < 0 & phi_variance > -1e-12] <- 0
    phi_variance <- pmax(phi_variance, 0)
    phi_sd <- sqrt(phi_variance)

    theta_variance <- theta_second_sum / n_draws - theta_mean^2
    theta_variance[theta_variance < 0 & theta_variance > -1e-12] <- 0
    theta_variance <- pmax(theta_variance, 0)
    theta_sd <- sqrt(theta_variance)
  }

  legacy_with_prior <- .cdseq2_row_normalize(
    z_mean + beta,
    object_name = "legacy posterior ctsGEP estimate"
  )

  legacy_allocation_only <- matrix(
    NA_real_,
    nrow = n_types,
    ncol = n_genes
  )
  positive_mean_depth <- rowSums(z_mean) > 0
  if (any(positive_mean_depth)) {
    legacy_allocation_only[positive_mean_depth, ] <- sweep(
      z_mean[positive_mean_depth, , drop = FALSE],
      1L,
      rowSums(z_mean)[positive_mean_depth],
      "/"
    )
  }

  draw_normalized_allocation <- matrix(
    NA_real_,
    nrow = n_types,
    ncol = n_genes
  )
  valid_allocation <- allocation_profile_n > 0
  if (any(valid_allocation)) {
    draw_normalized_allocation[valid_allocation, ] <- sweep(
      allocation_profile_sum[valid_allocation, , drop = FALSE],
      1L,
      allocation_profile_n[valid_allocation],
      "/"
    )
  }

  theta_allocation <- ssp_mean / sum(ssp_mean)

  rownames(phi_mean) <- rownames(legacy_with_prior) <-
    rownames(legacy_allocation_only) <-
    rownames(draw_normalized_allocation) <- cell_types
  colnames(phi_mean) <- colnames(legacy_with_prior) <-
    colnames(legacy_allocation_only) <-
    colnames(draw_normalized_allocation) <- gene_symbols
  rownames(z_mean) <- cell_types
  colnames(z_mean) <- gene_symbols

  if (!is.null(phi_sd)) {
    rownames(phi_sd) <- cell_types
    colnames(phi_sd) <- gene_symbols
  }

  depth_quantiles <- t(apply(
    depth_draws,
    1L,
    stats::quantile,
    probs = c(0.05, 0.50, 0.95),
    names = FALSE,
    type = 8
  ))
  colnames(depth_quantiles) <- c("depth_q05", "depth_q50", "depth_q95")
  rownames(depth_quantiles) <- cell_types

  rho_matrix <- matrix(rho, nrow = n_types, ncol = n_draws)
  prior_weight_draws <- rho_matrix / (rho_matrix + depth_draws)
  mean_prior_weight <- rowMeans(prior_weight_draws)

  depth_flag <- if (is.null(low_information_depth)) {
    rep(NA, n_types)
  } else {
    rowMeans(depth_draws) < low_information_depth
  }
  prior_flag <- if (is.null(low_information_prior_weight)) {
    rep(NA, n_types)
  } else {
    mean_prior_weight > low_information_prior_weight
  }
  low_information <-
    ifelse(is.na(depth_flag), FALSE, depth_flag) |
    ifelse(is.na(prior_flag), FALSE, prior_flag)

  information <- data.frame(
    cell_type = cell_types,
    estimated_rna_fraction = theta_mean,
    estimated_rna_fraction_sd = if (is.null(theta_sd)) NA_real_ else theta_sd,
    estimated_rna_fraction_cv = if (is.null(theta_sd)) {
      NA_real_
    } else {
      ifelse(theta_mean > 0, theta_sd / theta_mean, NA_real_)
    },
    mean_allocated_reads = rowMeans(depth_draws),
    sd_allocated_reads = apply(depth_draws, 1L, stats::sd),
    depth_q05 = depth_quantiles[, "depth_q05"],
    depth_q50 = depth_quantiles[, "depth_q50"],
    depth_q95 = depth_quantiles[, "depth_q95"],
    reference_precision = rho,
    mean_reference_prior_weight = mean_prior_weight,
    fraction_draws_with_positive_depth = allocation_profile_n / n_draws,
    below_depth_threshold = depth_flag,
    above_prior_weight_threshold = prior_flag,
    low_information = low_information,
    stringsAsFactors = FALSE
  )
  rownames(information) <- cell_types

  if (!is.null(phi_draw_means)) {
    dimnames(phi_draw_means) <- list(
      cell_types,
      gene_symbols,
      dimnames(z_array)[[3L]]
    )
  }

  list(
    estGEP_RB = t(phi_mean),
    estGEP_RB_sd = if (is.null(phi_sd)) NULL else t(phi_sd),
    estGEP_legacy = t(legacy_with_prior),
    estGEP_allocation = t(legacy_allocation_only),
    estGEP_draw_normalized_allocation = t(draw_normalized_allocation),
    estProp_RB = matrix(
      theta_mean,
      nrow = n_types,
      ncol = 1L,
      dimnames = list(cell_types, NULL)
    ),
    estProp_RB_sd = if (is.null(theta_sd)) {
      NULL
    } else {
      matrix(
        theta_sd,
        nrow = n_types,
        ncol = 1L,
        dimnames = list(cell_types, NULL)
      )
    },
    estProp_allocation = matrix(
      theta_allocation,
      nrow = n_types,
      ncol = 1L,
      dimnames = list(cell_types, NULL)
    ),
    allocation_mean = t(z_mean),
    assignment_mean = matrix(
      ssp_mean,
      nrow = n_types,
      ncol = 1L,
      dimnames = list(cell_types, NULL)
    ),
    phi_conditional_mean_array = phi_draw_means,
    information = information,
    internal = list(
      depth_draws = depth_draws,
      theta_conditional_draws = theta_conditional_draws,
      maximum_assignment_mismatch = maximum_assignment_mismatch,
      bulk_total_by_draw = bulk_total_by_draw
    )
  )
}


# Exact conditional first and second moments for a broad ctsGEP obtained by
# aggregating fitted cell states.  This keeps the fitted model unchanged.
.cdseq2_aggregate_states <- function(z_array,
                                     ssp_array,
                                     beta,
                                     alpha,
                                     state_to_major,
                                     gene_symbols,
                                     return_posterior_sd,
                                     return_arrays,
                                     low_information_depth,
                                     low_information_prior_weight) {
  states <- dimnames(z_array)[[1L]]
  major_types <- unique(as.character(state_to_major[states]))
  n_major <- length(major_types)
  n_genes <- dim(z_array)[2L]
  n_draws <- dim(z_array)[3L]
  rho <- rowSums(beta)

  phi_sum <- matrix(0, nrow = n_major, ncol = n_genes)
  theta_sum <- numeric(n_major)
  z_sum <- matrix(0, nrow = n_major, ncol = n_genes)
  ssp_sum <- numeric(n_major)
  depth_draws <- matrix(0, nrow = n_major, ncol = n_draws)
  weighted_prior_draws <- matrix(0, nrow = n_major, ncol = n_draws)
  allocation_profile_sum <- matrix(0, nrow = n_major, ncol = n_genes)
  allocation_profile_n <- numeric(n_major)

  phi_second_sum <- if (return_posterior_sd) {
    matrix(0, nrow = n_major, ncol = n_genes)
  } else {
    NULL
  }
  theta_second_sum <- if (return_posterior_sd) numeric(n_major) else NULL

  z_major_array <- if (return_arrays) {
    array(
      0,
      dim = c(n_major, n_genes, n_draws),
      dimnames = list(major_types, gene_symbols, dimnames(z_array)[[3L]])
    )
  } else {
    NULL
  }
  ssp_major_array <- if (return_arrays) {
    array(
      0,
      dim = c(n_major, 1L, n_draws),
      dimnames = list(major_types, dimnames(ssp_array)[[2L]], dimnames(z_array)[[3L]])
    )
  } else {
    NULL
  }

  for (b in seq_len(n_draws)) {
    z_b <- z_array[, , b, drop = FALSE][, , 1L]
    ssp_b <- as.numeric(ssp_array[, 1L, b])
    state_depth <- rowSums(z_b)

    state_phi_alpha <- z_b + beta
    state_phi_total <- rowSums(state_phi_alpha)
    state_phi_mean <- sweep(
      state_phi_alpha,
      1L,
      state_phi_total,
      "/"
    )
    state_phi_variance <- sweep(
      state_phi_mean * (1 - state_phi_mean),
      1L,
      state_phi_total + 1,
      "/"
    )

    state_theta_alpha <- ssp_b + alpha
    theta_total <- sum(state_theta_alpha)

    for (k in seq_len(n_major)) {
      idx <- which(state_to_major[states] == major_types[k])
      group_theta_alpha <- state_theta_alpha[idx]
      group_theta_total <- sum(group_theta_alpha)
      group_theta_mean <- group_theta_total / theta_total

      state_weights <- group_theta_alpha / group_theta_total
      group_phi_mean <- colSums(
        sweep(
          state_phi_mean[idx, , drop = FALSE],
          1L,
          state_weights,
          "*"
        )
      )

      group_z <- colSums(z_b[idx, , drop = FALSE])
      group_ssp <- sum(ssp_b[idx])
      group_depth <- sum(state_depth[idx])

      phi_sum[k, ] <- phi_sum[k, ] + group_phi_mean
      theta_sum[k] <- theta_sum[k] + group_theta_mean
      z_sum[k, ] <- z_sum[k, ] + group_z
      ssp_sum[k] <- ssp_sum[k] + group_ssp
      depth_draws[k, b] <- group_depth
      weighted_prior_draws[k, b] <- sum(
        state_weights * rho[idx] / (rho[idx] + state_depth[idx])
      )

      if (group_depth > 0) {
        allocation_profile_sum[k, ] <-
          allocation_profile_sum[k, ] + group_z / group_depth
        allocation_profile_n[k] <- allocation_profile_n[k] + 1
      }

      if (return_posterior_sd) {
        term_phi_uncertainty <- colSums(
          sweep(
            state_phi_variance[idx, , drop = FALSE],
            1L,
            group_theta_alpha * (group_theta_alpha + 1),
            "*"
          )
        )
        term_state_weight_uncertainty <- colSums(
          sweep(
            state_phi_mean[idx, , drop = FALSE]^2,
            1L,
            group_theta_alpha,
            "*"
          )
        )
        group_phi_variance <-
          (
            term_phi_uncertainty +
              term_state_weight_uncertainty -
              group_theta_total * group_phi_mean^2
          ) /
          (group_theta_total * (group_theta_total + 1))
        group_phi_variance[
          group_phi_variance < 0 & group_phi_variance > -1e-12
        ] <- 0
        group_phi_variance <- pmax(group_phi_variance, 0)

        group_theta_variance <- group_theta_mean *
          (1 - group_theta_mean) / (theta_total + 1)

        phi_second_sum[k, ] <- phi_second_sum[k, ] +
          group_phi_variance + group_phi_mean^2
        theta_second_sum[k] <- theta_second_sum[k] +
          group_theta_variance + group_theta_mean^2
      }

      if (return_arrays) {
        z_major_array[k, , b] <- group_z
        ssp_major_array[k, 1L, b] <- group_ssp
      }
    }
  }

  phi_mean <- phi_sum / n_draws
  theta_mean <- theta_sum / n_draws
  z_mean <- z_sum / n_draws
  ssp_mean <- ssp_sum / n_draws

  phi_sd <- theta_sd <- NULL
  if (return_posterior_sd) {
    phi_variance <- phi_second_sum / n_draws - phi_mean^2
    phi_variance[phi_variance < 0 & phi_variance > -1e-12] <- 0
    phi_variance <- pmax(phi_variance, 0)
    phi_sd <- sqrt(phi_variance)

    theta_variance <- theta_second_sum / n_draws - theta_mean^2
    theta_variance[theta_variance < 0 & theta_variance > -1e-12] <- 0
    theta_variance <- pmax(theta_variance, 0)
    theta_sd <- sqrt(theta_variance)
  }

  allocation_profile <- matrix(
    NA_real_,
    nrow = n_major,
    ncol = n_genes
  )
  positive_mean_depth <- rowSums(z_mean) > 0
  if (any(positive_mean_depth)) {
    allocation_profile[positive_mean_depth, ] <- sweep(
      z_mean[positive_mean_depth, , drop = FALSE],
      1L,
      rowSums(z_mean)[positive_mean_depth],
      "/"
    )
  }

  draw_normalized_allocation <- matrix(
    NA_real_,
    nrow = n_major,
    ncol = n_genes
  )
  valid_allocation <- allocation_profile_n > 0
  if (any(valid_allocation)) {
    draw_normalized_allocation[valid_allocation, ] <- sweep(
      allocation_profile_sum[valid_allocation, , drop = FALSE],
      1L,
      allocation_profile_n[valid_allocation],
      "/"
    )
  }

  theta_allocation <- ssp_mean / sum(ssp_mean)

  rownames(phi_mean) <- rownames(allocation_profile) <-
    rownames(draw_normalized_allocation) <- rownames(z_mean) <- major_types
  colnames(phi_mean) <- colnames(allocation_profile) <-
    colnames(draw_normalized_allocation) <- colnames(z_mean) <- gene_symbols
  if (!is.null(phi_sd)) {
    rownames(phi_sd) <- major_types
    colnames(phi_sd) <- gene_symbols
  }

  depth_quantiles <- t(apply(
    depth_draws,
    1L,
    stats::quantile,
    probs = c(0.05, 0.50, 0.95),
    names = FALSE,
    type = 8
  ))
  colnames(depth_quantiles) <- c("depth_q05", "depth_q50", "depth_q95")
  rownames(depth_quantiles) <- major_types

  mean_prior_weight <- rowMeans(weighted_prior_draws)
  depth_flag <- if (is.null(low_information_depth)) {
    rep(NA, n_major)
  } else {
    rowMeans(depth_draws) < low_information_depth
  }
  prior_flag <- if (is.null(low_information_prior_weight)) {
    rep(NA, n_major)
  } else {
    mean_prior_weight > low_information_prior_weight
  }
  low_information <-
    ifelse(is.na(depth_flag), FALSE, depth_flag) |
    ifelse(is.na(prior_flag), FALSE, prior_flag)

  information <- data.frame(
    cell_type = major_types,
    estimated_rna_fraction = theta_mean,
    estimated_rna_fraction_sd = if (is.null(theta_sd)) NA_real_ else theta_sd,
    estimated_rna_fraction_cv = if (is.null(theta_sd)) {
      NA_real_
    } else {
      ifelse(theta_mean > 0, theta_sd / theta_mean, NA_real_)
    },
    mean_allocated_reads = rowMeans(depth_draws),
    sd_allocated_reads = apply(depth_draws, 1L, stats::sd),
    depth_q05 = depth_quantiles[, "depth_q05"],
    depth_q50 = depth_quantiles[, "depth_q50"],
    depth_q95 = depth_quantiles[, "depth_q95"],
    mean_reference_prior_weight = mean_prior_weight,
    fraction_draws_with_positive_depth = allocation_profile_n / n_draws,
    below_depth_threshold = depth_flag,
    above_prior_weight_threshold = prior_flag,
    low_information = low_information,
    stringsAsFactors = FALSE
  )
  rownames(information) <- major_types

  list(
    estGEP = t(phi_mean),
    estGEP_RB = t(phi_mean),
    estGEP_RB_sd = if (is.null(phi_sd)) NULL else t(phi_sd),
    estGEP2 = t(allocation_profile),
    estGEP_allocation = t(allocation_profile),
    estGEP_draw_normalized_allocation = t(draw_normalized_allocation),
    estProp = matrix(
      theta_mean,
      nrow = n_major,
      ncol = 1L,
      dimnames = list(major_types, NULL)
    ),
    estProp_sd = if (is.null(theta_sd)) {
      NULL
    } else {
      matrix(
        theta_sd,
        nrow = n_major,
        ncol = 1L,
        dimnames = list(major_types, NULL)
      )
    },
    estProp2 = matrix(
      theta_allocation,
      nrow = n_major,
      ncol = 1L,
      dimnames = list(major_types, NULL)
    ),
    estGEP_mat = t(z_mean),
    estSSp_mat = matrix(
      ssp_mean,
      nrow = n_major,
      ncol = 1L,
      dimnames = list(major_types, NULL)
    ),
    estGEP_array = z_major_array,
    estSSP_array = ssp_major_array,
    information = information,
    state_to_major = state_to_major,
    major_cell_types = major_types
  )
}


#' CDSeq2: Reference-informed individual-level complete deconvolution
#'
#' Performs reference-informed complete deconvolution of bulk RNA-seq count
#' data. Each bulk sample is analyzed independently under the current CDSeq2
#' Dirichlet-multinomial model, while the single-cell reference supplies the
#' prior mean for cell-type-specific gene-expression profiles (ctsGEPs).
#'
#' @details
#' `cdseq2()` preserves the existing CDSeq2 generative model and Gibbs sampler,
#' but adds model-preserving improvements to reference construction, posterior
#' summarization, reproducibility, diagnostics, and optional post-hoc
#' aggregation of fitted cell states.
#'
#' For a retained allocation draw, let `Z[t, g]` be the number of reads assigned
#' to gene `g` and modeled cell type or state `t`, and let `beta[t, g]` be the
#' corresponding Dirichlet prior parameter. The recommended normalized ctsGEP
#' estimator is the Rao-Blackwellized posterior mean obtained by averaging
#'
#' `E(phi[t, g] | Z) = (Z[t, g] + beta[t, g]) /
#' (sum_g Z[t, g] + sum_g beta[t, g])`
#'
#' over all retained draws and chains. This estimator is returned as
#' `estGEP_RB`. The posterior mean latent allocation counts are returned
#' separately as `estGEP_count`; they are not normalized ctsGEP parameters.
#'
#' All modeled proportions are RNA/read fractions under the multinomial
#' likelihood. They should not be interpreted automatically as cell-count
#' fractions when cell types differ in total RNA content.
#'
#' @param bulk_data A numeric base R matrix of unnormalized, nonnegative integer
#'   bulk RNA-seq counts with genes in rows and samples in columns. Unique,
#'   nonempty gene names and sample names are required. Each retained sample
#'   must have positive total count over the modeled gene panel.
#'
#' @param mcmc_iterations A positive integer greater than one giving the total
#'   number of Gibbs iterations run for each sample and each chain. Only the
#'   final `save_last_n_mc` allocation states are retained.
#'
#' @param cpu_number `NULL` or a positive integer giving the number of parallel
#'   worker processes used across bulk samples. `NULL` and `1` use sequential
#'   execution. Parallelization is across samples, not across chains within a
#'   sample.
#'
#' @param reference_scRNA_count A numeric base R matrix of unnormalized,
#'   nonnegative integer single-cell or single-nucleus RNA-seq counts with genes
#'   in rows and reference cells in columns. Gene row names must be unique and
#'   must match the bulk genes after applying `genes_use` and `exclude_genes`.
#'
#' @param reference_scRNA_annotation A character or factor vector of length
#'   `ncol(reference_scRNA_count)` containing the modeled cell-type or cell-state
#'   label for every reference cell. The order of first appearance determines
#'   the modeled-label order used for unnamed parameter vectors.
#'
#' @param reference_scRNA_donor An optional character or factor vector of length
#'   `ncol(reference_scRNA_count)` containing the donor identifier for every
#'   reference cell. It is required when
#'   `reference_aggregation = "equal_donor"` and ignored for pooled-count
#'   aggregation.
#'
#' @param reference_aggregation A character string specifying how the
#'   reference prior mean is constructed. `"pooled_counts"` sums all reference
#'   counts within each modeled label and then normalizes. `"equal_donor"`
#'   constructs one donor-by-label pseudobulk profile, normalizes each profile,
#'   and averages eligible donor profiles with equal donor weight. The default
#'   is `"pooled_counts"`.
#'
#' @param reference_min_cells_per_donor_type A positive integer giving the
#'   minimum number of reference cells required for a donor-by-modeled-label
#'   pseudobulk when `reference_aggregation = "equal_donor"`. Donor-label groups
#'   below the threshold are discarded. Every modeled label must retain at
#'   least one eligible donor.
#'
#' @param reference_proportion A positive numeric scalar, an unnamed vector in
#'   modeled-label order, or a named vector indexed by modeled labels. The
#'   expanded values are the Dirichlet concentration parameters for cell-type
#'   RNA proportions. Their normalized values define the prior mean, and their
#'   sum defines the total prior concentration.
#'
#' @param reference_precision A positive numeric scalar, an unnamed vector in
#'   modeled-label order, or a named vector indexed by modeled labels. For label
#'   `t`, the value is the total Dirichlet concentration
#'   `rho[t] = sum_g beta[t, g]` for its normalized ctsGEP. Larger values impose
#'   stronger shrinkage toward the reference profile; smaller values allow more
#'   sample-specific deviation.
#'
#' @param smooth_ref A single logical value indicating whether positive
#'   pseudocount smoothing is applied when constructing the reference profile.
#'   Smoothing is generally required when any modeled gene has zero reference
#'   count in a modeled label because all Dirichlet parameters must be strictly
#'   positive.
#'
#' @param smoothing_alpha A nonnegative numeric scalar, an unnamed vector in
#'   modeled-label order, or a named vector indexed by modeled labels. It gives
#'   the additive pseudocount used for reference smoothing. Values must be
#'   strictly positive when `smooth_ref = TRUE`.
#'
#' @param smoothing_scope A character string controlling where smoothing is
#'   applied. `"aggregate"` adds one pseudocount per gene after pooling within a
#'   modeled label and is recommended. `"per_cell"` reproduces the legacy
#'   behavior in which the effective aggregate pseudocount is multiplied by the
#'   number of contributing reference cells. The default is `"aggregate"`.
#'
#' @param genes_use `NULL` or a unique character vector specifying the genes to
#'   model and their order. When `NULL`, bulk and reference matrices must contain
#'   the same gene set; the reference is reordered to the bulk row order. When
#'   supplied, every selected gene must occur in both matrices.
#'
#' @param exclude_genes `NULL` or a character vector of genes to remove after
#'   alignment and optional `genes_use` selection. Genes absent from the aligned
#'   matrices are ignored.
#'
#' @param drop_all_zero_genes A single logical value. When `TRUE`, genes having
#'   zero total count in both the complete bulk matrix and the complete
#'   reference matrix are removed before reference construction and sampling.
#'
#' @param aggregate_to `NULL`, an unnamed character vector with one entry per
#'   modeled cell state, or a named character vector whose names are modeled
#'   states and whose values are final major-cell-type labels. When supplied,
#'   sampling is still performed at the original state resolution and posterior
#'   results are subsequently aggregated to the requested major types. This does
#'   not alter the fitted model. The mapping must cover every modeled state.
#'
#' @param rng A scalar equal to `1` or `2`, selecting the random-number generator
#'   used by `gibbsSampler`: `1` for Mersenne Twister and `2` for
#'   Xoshiro256plus. `set.seed()` is applied before each chain when `seed` is not
#'   `NULL`; exact reproducibility for `rng = 2` also depends on how the compiled
#'   Xoshiro implementation is seeded.
#'
#' @param save_last_n_mc A positive integer not exceeding `mcmc_iterations`.
#'   This many final collapsed allocation states are retained from each chain.
#'   The implicit burn-in is `mcmc_iterations - save_last_n_mc`. A warning is
#'   issued when all iterations are retained.
#'
#' @param n_chains A positive integer giving the number of independent Gibbs
#'   chains fitted for every bulk sample. Retained draws from all chains are
#'   combined for posterior summaries, and between-chain diagnostics are
#'   returned when `n_chains > 1`.
#'
#' @param seed `NULL` or one finite numeric value used as the base seed for
#'   deterministic sample-by-chain seeds. The derived seeds are returned in
#'   `parameters$chain_seeds`. `NULL` leaves the R random-number seed unset by
#'   this function.
#'
#' @param primary_gep_estimator A character string controlling only which
#'   normalized estimate is exposed under the backward-compatible name
#'   `estGEP`. `"rao_blackwell"` uses `estGEP_RB` and is recommended.
#'   `"legacy"` uses the legacy sampler-level normalized estimate. Explicitly
#'   named versions of both estimators remain available in the returned object.
#'
#' @param return_posterior_sd A single logical value. When `TRUE`, return
#'   posterior standard deviations for normalized ctsGEPs and RNA proportions.
#'   The calculation uses the law of total variance and therefore combines
#'   conditional Dirichlet uncertainty with variation across retained allocation
#'   states.
#'
#' @param return_mcmc_arrays A single logical value indicating whether retained
#'   collapsed allocation arrays are included in each sample result. These
#'   arrays can require substantial memory for large gene panels, many samples,
#'   many chains, or large `save_last_n_mc`.
#'
#' @param return_phi_draw_means A single logical value. When `TRUE`, return the
#'   conditional posterior mean of each normalized ctsGEP for every retained
#'   allocation state as `estPhi_conditional_mean_array`. This may require
#'   substantial additional memory.
#'
#' @param return_assignments A single logical value. When `TRUE`, retain
#'   `cellTypeAssignSplit` and `cellTypeBag` from the final fitted chain for
#'   compatibility with the existing sampler interface. These objects are not
#'   combined across chains.
#'
#' @param return_aggregated_arrays A single logical value. When `TRUE` and
#'   `aggregate_to` is supplied, include major-type allocation arrays in the
#'   post-hoc aggregated result. This can substantially increase memory use.
#'
#' @param low_information_depth `NULL`, a nonnegative scalar, or a
#'   label-specific nonnegative vector. A sample-by-label estimate is flagged
#'   when its posterior mean allocated read depth is below this threshold. When
#'   `aggregate_to` is used, a scalar threshold is the safest choice because the
#'   same argument is evaluated at both state and major-type resolutions.
#'
#' @param low_information_prior_weight `NULL`, a scalar in `[0, 1]`, or a
#'   label-specific vector in `[0, 1]`. A sample-by-label estimate is flagged
#'   when its mean reference-prior weight, approximately
#'   `rho / (rho + allocated depth)`, exceeds this threshold. When
#'   `aggregate_to` is used, a scalar threshold is the safest choice.
#'
#' @param reference_correlation_warning `NULL` or a scalar in `[0, 1]`.
#'   Pairwise Pearson correlations between modeled reference profiles are
#'   calculated across genes. A warning is issued for pairs with correlation at
#'   least this value because highly similar profiles may produce spillover.
#'   Set to `NULL` to disable the warning.
#'
#' @param adjust_prior A single logical value enabling an experimental,
#'   phenotype-blind two-pass reference adaptation. A first fit to all supplied
#'   bulk samples generates a cohort-level profile, which is blended with the
#'   external reference before the individual sample fits. This option does not
#'   use phenotype labels, but each bulk sample contributes to the adapted prior;
#'   it should therefore be benchmarked carefully and used cautiously in formal
#'   validation analyses.
#'
#' @param adjust_prior_weight A scalar in `[0, 1]` giving the convex weight of
#'   the first-pass cohort profile when `adjust_prior = TRUE`. Zero leaves the
#'   external reference unchanged; one replaces it with the first-pass profile.
#'   The default is `0.2`.
#'
#' @param verbose A single logical value controlling progress messages from
#'   `cdseq2()` and the underlying Gibbs sampler.
#'
#' @section Reference construction:
#' With `reference_aggregation = "pooled_counts"`, cells contribute in
#' proportion to their total observed counts. With
#' `reference_aggregation = "equal_donor"`, eligible donor-by-label profiles
#' contribute equally after within-profile normalization. The latter can reduce
#' domination by donors having many cells or unusually deep libraries, but it
#' requires valid donor identifiers and sufficient cells for every modeled
#' label.
#'
#' @section ctsGEP output semantics:
#' `estGEP_RB` is the recommended genes-by-modeled-label posterior mean of the
#' normalized ctsGEP parameter, and each column sums to one. `estGEP_count` is
#' the genes-by-modeled-label posterior mean of latent allocated counts and is
#' not normalized. `estGEP_allocation` normalizes the posterior mean allocation
#' counts. `estGEP_draw_normalized_allocation` averages draw-wise normalized
#' allocations over draws having positive depth. `estGEP_legacy` preserves the
#' legacy normalized estimate for controlled comparisons. `estGEP` is an alias
#' selected by `primary_gep_estimator`.
#'
#' @section Low-information diagnostics:
#' The per-sample `information` table reports estimated RNA fraction, posterior
#' uncertainty, allocated-depth summaries, the effective reference-prior weight,
#' and threshold flags. Thresholds are diagnostic only: the function does not
#' remove, mask, or reweight estimates automatically.
#'
#' @return A named list with the following top-level components:
#'
#' \describe{
#'   \item{cdseq_result}{A named list with one element per bulk sample. Important
#'   components include `estGEP_RB`, `estGEP_RB_sd`, `estGEP_count`, `estProp_RB`,
#'   `estProp_RB_sd`, alternative and legacy estimators, allocation arrays,
#'   low-information diagnostics, chain diagnostics, sampler consistency checks,
#'   and optional state-to-major-type aggregated results. Normalized ctsGEP and
#'   allocated-count matrices have genes in rows and modeled labels in columns.}
#'   \item{parameters}{The final Dirichlet parameters, original and final
#'   reference profiles, expanded precision and proportion parameters, reference
#'   construction diagnostics, selected genes, state-to-major mapping, MCMC
#'   settings, and deterministic chain seeds.}
#'   \item{running_time}{A list containing the Gibbs-running-time object and
#'   total elapsed time in hours.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- cdseq2(
#'   bulk_data = bulk_counts,
#'   reference_scRNA_count = reference_counts,
#'   reference_scRNA_annotation = reference_cell_type,
#'   reference_precision = 1e6,
#'   mcmc_iterations = 1000,
#'   save_last_n_mc = 500,
#'   n_chains = 3,
#'   seed = 123,
#'   primary_gep_estimator = "rao_blackwell"
#' )
#'
#' # Recommended normalized ctsGEP for the first sample:
#' phi_hat <- fit$cdseq_result[[1]]$estGEP_RB
#'
#' # Posterior mean cell-type-specific allocated counts:
#' z_hat <- fit$cdseq_result[[1]]$estGEP_count
#'
#' # Equal-donor reference construction:
#' fit_equal_donor <- cdseq2(
#'   bulk_data = bulk_counts,
#'   reference_scRNA_count = reference_counts,
#'   reference_scRNA_annotation = reference_cell_type,
#'   reference_scRNA_donor = reference_donor,
#'   reference_aggregation = "equal_donor",
#'   reference_min_cells_per_donor_type = 10
#' )
#'
#' # Fit moderate-resolution states and aggregate afterward:
#' state_to_major <- c(
#'   Ast_homeostatic = "Astrocytes",
#'   Ast_reactive = "Astrocytes",
#'   Exc_upper = "Excitatory",
#'   Exc_deep = "Excitatory"
#' )
#' fit_states <- cdseq2(
#'   bulk_data = bulk_counts,
#'   reference_scRNA_count = reference_counts,
#'   reference_scRNA_annotation = reference_state,
#'   aggregate_to = state_to_major
#' )
#' major_phi <- fit_states$cdseq_result[[1]]$aggregated$estGEP_RB
#' }
#'
#' @seealso \code{\link{cdseq}} for the original user-facing implementation.
#'
#' @md
#' @md
#' @importFrom foreach %dopar%
#' @export
cdseq2 <- function(bulk_data,
                   mcmc_iterations = 700,
                   cpu_number = 1,
                   reference_scRNA_count = NULL,
                   reference_scRNA_annotation = NULL,
                   reference_scRNA_donor = NULL,
                   reference_aggregation = c("pooled_counts", "equal_donor"),
                   reference_min_cells_per_donor_type = 1L,
                   reference_proportion = 1,
                   reference_precision = 1e10,
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
                   verbose = FALSE) {

  cat(
    "cdseq2 - reference-based complete deconvolution and ",
    "cell-type-specific expression analysis.\n",
    sep = ""
  )
  start_time <- Sys.time()

  reference_aggregation <- match.arg(reference_aggregation)
  smoothing_scope <- match.arg(smoothing_scope)
  primary_gep_estimator <- match.arg(primary_gep_estimator)

  .cdseq2_validate_count_matrix(
    bulk_data,
    "bulk_data",
    require_colnames = TRUE
  )
  .cdseq2_validate_count_matrix(
    reference_scRNA_count,
    "reference_scRNA_count",
    require_colnames = FALSE
  )

  if (is.null(reference_scRNA_annotation)) {
    stop("Please provide reference_scRNA_annotation.", call. = FALSE)
  }
  if (length(reference_scRNA_annotation) != ncol(reference_scRNA_count)) {
    stop(
      "reference_scRNA_annotation must have length ncol(reference_scRNA_count).",
      call. = FALSE
    )
  }
  if (anyNA(reference_scRNA_annotation) ||
      any(!nzchar(as.character(reference_scRNA_annotation)))) {
    stop("reference_scRNA_annotation contains missing or empty labels.", call. = FALSE)
  }

  if (!is.null(reference_scRNA_donor)) {
    if (length(reference_scRNA_donor) != ncol(reference_scRNA_count)) {
      stop(
        "reference_scRNA_donor must have length ncol(reference_scRNA_count).",
        call. = FALSE
      )
    }
    if (anyNA(reference_scRNA_donor) ||
        any(!nzchar(as.character(reference_scRNA_donor)))) {
      stop("reference_scRNA_donor contains missing or empty labels.", call. = FALSE)
    }
  }

  if (!is.numeric(mcmc_iterations) || length(mcmc_iterations) != 1L ||
      !is.finite(mcmc_iterations) || mcmc_iterations < 2) {
    stop("mcmc_iterations must be an integer greater than 1.", call. = FALSE)
  }
  mcmc_iterations <- as.integer(ceiling(mcmc_iterations))

  if (!is.numeric(save_last_n_mc) || length(save_last_n_mc) != 1L ||
      !is.finite(save_last_n_mc) || save_last_n_mc < 1) {
    stop("save_last_n_mc must be a positive integer.", call. = FALSE)
  }
  save_last_n_mc <- as.integer(ceiling(save_last_n_mc))
  if (save_last_n_mc > mcmc_iterations) {
    stop("save_last_n_mc cannot exceed mcmc_iterations.", call. = FALSE)
  }
  if (save_last_n_mc == mcmc_iterations) {
    warning(
      "All MCMC iterations are being retained; no iterations are excluded as burn-in.",
      call. = FALSE
    )
  }

  if (!is.numeric(n_chains) || length(n_chains) != 1L ||
      !is.finite(n_chains) || n_chains < 1) {
    stop("n_chains must be a positive integer.", call. = FALSE)
  }
  n_chains <- as.integer(ceiling(n_chains))

  if (!is.numeric(rng) || length(rng) != 1L || !(rng %in% c(1, 2))) {
    stop("rng must be 1 (Mersenne Twister) or 2 (Xoshiro256plus).", call. = FALSE)
  }
  rng <- as.integer(rng)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or one finite numeric value.", call. = FALSE)
    }
    seed <- as.integer(abs(seed) %% (.Machine$integer.max - 1L))
    if (seed == 0L) seed <- 1L
  }
  if (rng == 2L && !is.null(seed)) {
    warning(
      "set.seed() is applied before each chain, but reproducibility for rng = 2 ",
      "depends on how the package's Xoshiro generator is seeded internally.",
      call. = FALSE
    )
  }

  if (!is.numeric(reference_min_cells_per_donor_type) ||
      length(reference_min_cells_per_donor_type) != 1L ||
      !is.finite(reference_min_cells_per_donor_type) ||
      reference_min_cells_per_donor_type < 1) {
    stop(
      "reference_min_cells_per_donor_type must be a positive integer.",
      call. = FALSE
    )
  }
  reference_min_cells_per_donor_type <- as.integer(
    ceiling(reference_min_cells_per_donor_type)
  )

  logical_arguments <- list(
    smooth_ref = smooth_ref,
    drop_all_zero_genes = drop_all_zero_genes,
    return_posterior_sd = return_posterior_sd,
    return_mcmc_arrays = return_mcmc_arrays,
    return_phi_draw_means = return_phi_draw_means,
    return_assignments = return_assignments,
    return_aggregated_arrays = return_aggregated_arrays,
    adjust_prior = adjust_prior,
    verbose = verbose
  )
  invalid_logical <- names(logical_arguments)[
    !vapply(
      logical_arguments,
      function(x) is.logical(x) && length(x) == 1L && !is.na(x),
      logical(1L)
    )
  ]
  if (length(invalid_logical) > 0L) {
    stop(
      "The following arguments must be TRUE or FALSE: ",
      paste(invalid_logical, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(adjust_prior_weight) || length(adjust_prior_weight) != 1L ||
      !is.finite(adjust_prior_weight) || adjust_prior_weight < 0 ||
      adjust_prior_weight > 1) {
    stop("adjust_prior_weight must lie in [0, 1].", call. = FALSE)
  }

  if (!is.null(reference_correlation_warning)) {
    if (!is.numeric(reference_correlation_warning) ||
        length(reference_correlation_warning) != 1L ||
        !is.finite(reference_correlation_warning) ||
        reference_correlation_warning < 0 ||
        reference_correlation_warning > 1) {
      stop("reference_correlation_warning must be NULL or lie in [0,1].", call. = FALSE)
    }
  }

  # Align and optionally subset the modeling genes.
  bulk_genes_original <- rownames(bulk_data)
  reference_genes_original <- rownames(reference_scRNA_count)

  if (is.null(genes_use)) {
    if (!setequal(bulk_genes_original, reference_genes_original)) {
      stop(
        "bulk_data and reference_scRNA_count must contain the same genes, ",
        "or genes_use must explicitly select a shared panel.",
        call. = FALSE
      )
    }
    reference_scRNA_count <- reference_scRNA_count[
      match(bulk_genes_original, reference_genes_original),
      ,
      drop = FALSE
    ]
  } else {
    if (!is.character(genes_use) || anyNA(genes_use) ||
        any(!nzchar(genes_use)) || anyDuplicated(genes_use)) {
      stop("genes_use must be a unique, non-missing character vector.", call. = FALSE)
    }
    missing_bulk <- setdiff(genes_use, bulk_genes_original)
    missing_reference <- setdiff(genes_use, reference_genes_original)
    if (length(missing_bulk) > 0L || length(missing_reference) > 0L) {
      stop(
        "genes_use contains genes absent from bulk or reference data.",
        call. = FALSE
      )
    }
    bulk_data <- bulk_data[genes_use, , drop = FALSE]
    reference_scRNA_count <- reference_scRNA_count[genes_use, , drop = FALSE]
  }

  if (!is.null(exclude_genes)) {
    if (!is.character(exclude_genes) || anyNA(exclude_genes)) {
      stop("exclude_genes must be NULL or a character vector.", call. = FALSE)
    }
    keep <- !(rownames(bulk_data) %in% exclude_genes)
    bulk_data <- bulk_data[keep, , drop = FALSE]
    reference_scRNA_count <- reference_scRNA_count[keep, , drop = FALSE]
  }

  all_zero_removed <- character(0)
  if (drop_all_zero_genes) {
    all_zero <- rowSums(bulk_data) == 0 & rowSums(reference_scRNA_count) == 0
    if (any(all_zero)) {
      all_zero_removed <- rownames(bulk_data)[all_zero]
      bulk_data <- bulk_data[!all_zero, , drop = FALSE]
      reference_scRNA_count <- reference_scRNA_count[!all_zero, , drop = FALSE]
    }
  }

  if (nrow(bulk_data) < 2L) {
    stop("Fewer than two genes remain after filtering.", call. = FALSE)
  }
  if (any(colSums(bulk_data) <= 0)) {
    bad_samples <- colnames(bulk_data)[colSums(bulk_data) <= 0]
    stop(
      "The selected gene panel has zero total bulk counts for: ",
      paste(bad_samples, collapse = ", "),
      call. = FALSE
    )
  }

  gene_symbols <- rownames(bulk_data)
  sample_names <- colnames(bulk_data)
  annotation <- as.character(reference_scRNA_annotation)
  cell_types <- unique(annotation)
  cell_type_number <- length(cell_types)

  if (cell_type_number < 2L) {
    stop("At least two modeled cell types/states are required.", call. = FALSE)
  }

  reference_precision_vector <- .cdseq2_expand_parameter(
    reference_precision,
    cell_types,
    "reference_precision",
    strictly_positive = TRUE
  )
  reference_proportion_vector <- .cdseq2_expand_parameter(
    reference_proportion,
    cell_types,
    "reference_proportion",
    strictly_positive = TRUE
  )
  smoothing_alpha_vector <- .cdseq2_expand_parameter(
    smoothing_alpha,
    cell_types,
    "smoothing_alpha",
    strictly_positive = smooth_ref
  )

  low_information_depth_vector <- if (is.null(low_information_depth)) {
    NULL
  } else {
    .cdseq2_expand_parameter(
      low_information_depth,
      cell_types,
      "low_information_depth",
      strictly_positive = FALSE
    )
  }
  low_information_prior_weight_vector <-
    if (is.null(low_information_prior_weight)) {
      NULL
    } else {
      .cdseq2_expand_parameter(
        low_information_prior_weight,
        cell_types,
        "low_information_prior_weight",
        strictly_positive = FALSE,
        upper_bound = 1
      )
    }

  reference_object <- .cdseq2_build_reference(
    reference_scRNA_count = reference_scRNA_count,
    annotation = annotation,
    donor = reference_scRNA_donor,
    cell_types = cell_types,
    aggregation = reference_aggregation,
    min_cells_per_donor_type = reference_min_cells_per_donor_type,
    smooth_ref = smooth_ref,
    smoothing_alpha = smoothing_alpha_vector,
    smoothing_scope = smoothing_scope
  )

  reference_scRNA_original <- reference_object$profile
  reference_scRNA <- reference_scRNA_original
  beta <- sweep(reference_scRNA, 1L, reference_precision_vector, "*")
  alpha_initial <- reference_proportion_vector

  reference_correlation <- stats::cor(
    t(reference_scRNA),
    use = "pairwise.complete.obs"
  )
  if (!is.null(reference_correlation_warning) && cell_type_number > 1L) {
    high_pair_index <- which(
      upper.tri(reference_correlation) &
        reference_correlation >= reference_correlation_warning,
      arr.ind = TRUE
    )
    if (nrow(high_pair_index) > 0L) {
      high_pairs <- apply(high_pair_index, 1L, function(ix) {
        paste0(
          rownames(reference_correlation)[ix[1L]],
          " vs ",
          colnames(reference_correlation)[ix[2L]],
          " (r=",
          formatC(reference_correlation[ix[1L], ix[2L]], digits = 3, format = "f"),
          ")"
        )
      })
      warning(
        "Highly correlated modeled reference profiles may cause spillover: ",
        paste(high_pairs, collapse = "; "),
        call. = FALSE
      )
    }
  }

  state_to_major <- .cdseq2_prepare_aggregation_map(
    aggregate_to,
    cell_types
  )
  major_types <- if (is.null(state_to_major)) {
    NULL
  } else {
    unique(as.character(state_to_major[cell_types]))
  }

  low_information_depth_major <- if (is.null(state_to_major) ||
                                     is.null(low_information_depth)) {
    NULL
  } else {
    .cdseq2_expand_parameter(
      low_information_depth,
      major_types,
      "low_information_depth for aggregate_to",
      strictly_positive = FALSE
    )
  }
  low_information_prior_weight_major <-
    if (is.null(state_to_major) || is.null(low_information_prior_weight)) {
      NULL
    } else {
      .cdseq2_expand_parameter(
        low_information_prior_weight,
        major_types,
        "low_information_prior_weight for aggregate_to",
        strictly_positive = FALSE,
        upper_bound = 1
      )
    }

  verbose_int <- if (verbose) 1L else 0L
  printout <- 1L

  # Optional conservative two-pass reference update. This is an input/prior
  # update, not a change to the per-sample model.
  prior_adjustment <- NULL
  if (adjust_prior && adjust_prior_weight > 0) {
    if (verbose) message("Running experimental phenotype-blind prior adjustment...")
    if (!is.null(seed)) {
      adjustment_seed <- as.integer(
        ((as.double(seed) + 17) %% (.Machine$integer.max - 1)) + 1
      )
      set.seed(adjustment_seed)
    }

    alpha_full <- matrix(
      alpha_initial,
      nrow = cell_type_number,
      ncol = ncol(bulk_data)
    )
    adjust_result <- gibbsSampler(
      alpha_full,
      beta,
      bulk_data,
      cell_type_number,
      mcmc_iterations,
      printout,
      0,
      0,
      verbose_int,
      rng,
      min(5L, save_last_n_mc)
    )

    if (length(adjust_result$csGEP_vec) !=
        cell_type_number * nrow(bulk_data)) {
      stop(
        "gibbsSampler returned an unexpected ctsGEP length during prior adjustment.",
        call. = FALSE
      )
    }

    adjusted_gxt <- t(matrix(
      adjust_result$csGEP_vec,
      nrow = cell_type_number,
      ncol = nrow(bulk_data)
    ))
    adjusted_tg <- t(adjusted_gxt)

    adjusted_profile <- reference_scRNA
    positive_type <- rowSums(adjusted_tg) > 0
    adjusted_profile[positive_type, ] <- sweep(
      adjusted_tg[positive_type, , drop = FALSE],
      1L,
      rowSums(adjusted_tg)[positive_type],
      "/"
    )

    reference_scRNA <-
      (1 - adjust_prior_weight) * reference_scRNA_original +
      adjust_prior_weight * adjusted_profile
    reference_scRNA <- .cdseq2_row_normalize(
      reference_scRNA,
      object_name = "adjusted reference profile"
    )
    beta <- sweep(reference_scRNA, 1L, reference_precision_vector, "*")

    prior_adjustment <- list(
      weight = adjust_prior_weight,
      first_pass_profile = adjusted_profile,
      final_profile = reference_scRNA
    )
  }

  individual_sample <- lapply(seq_len(ncol(bulk_data)), function(i) {
    bulk_data[, i, drop = FALSE]
  })
  names(individual_sample) <- sample_names

  # Deterministic, sample-by-chain seeds for reproducible parallel execution.
  seed_matrix <- matrix(
    NA_integer_,
    nrow = length(sample_names),
    ncol = n_chains,
    dimnames = list(sample_names, paste0("chain", seq_len(n_chains)))
  )
  if (!is.null(seed)) {
    for (i in seq_along(sample_names)) {
      for (ch in seq_len(n_chains)) {
        proposed <- (
          as.double(seed) +
            1000003 * i +
            9176 * ch
        ) %% (.Machine$integer.max - 1)
        seed_matrix[i, ch] <- as.integer(proposed + 1)
      }
    }
  }

  if (is.null(cpu_number)) cpu_number <- 1L
  if (!is.numeric(cpu_number) || length(cpu_number) != 1L ||
      !is.finite(cpu_number) || cpu_number < 1) {
    stop("cpu_number must be NULL or a positive integer.", call. = FALSE)
  }
  cpu_number <- as.integer(ceiling(cpu_number))

  if (cpu_number > 1L) {
    doParallel::registerDoParallel(cpu_number)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  helper_exports <- c(
    ".cdseq2_row_normalize",
    ".cdseq2_combine_3d_arrays",
    ".cdseq2_phi_mean_only",
    ".cdseq2_basic_rhat",
    ".cdseq2_chain_diagnostics",
    ".cdseq2_summarize_sample",
    ".cdseq2_aggregate_states",
    "gibbsSampler",
    "logpost",
    "loglikelihood"
  )
  helper_exports <- helper_exports[vapply(
    helper_exports,
    exists,
    logical(1L),
    envir = environment(),
    inherits = TRUE
  )]

  gibbs_start_time <- Sys.time()

  # Keep the file sourceable outside an installed package namespace.
  `%dopar%` <- foreach::`%dopar%`

  allresult <- foreach::foreach(
    i = seq_along(individual_sample),
    .inorder = TRUE,
    .verbose = FALSE,
    .export = helper_exports
  ) %dopar% {
    sample_start_time <- Sys.time()
    sample_data <- individual_sample[[i]]
    sample_name <- colnames(sample_data)

    alpha <- matrix(
      alpha_initial,
      nrow = cell_type_number,
      ncol = 1L,
      dimnames = list(cell_types, sample_name)
    )

    chain_z <- vector("list", n_chains)
    chain_ssp <- vector("list", n_chains)
    chain_sampler_gep <- vector("list", n_chains)
    chain_sampler_ssp <- vector("list", n_chains)
    last_sampler_result <- NULL

    for (ch in seq_len(n_chains)) {
      if (!is.null(seed)) set.seed(seed_matrix[i, ch])

      sampler_result <- gibbsSampler(
        alpha,
        beta,
        sample_data,
        cell_type_number,
        mcmc_iterations,
        printout,
        0,
        0,
        verbose_int,
        rng,
        save_last_n_mc
      )

      expected_gep_length <- cell_type_number * nrow(sample_data)
      expected_gep_all_length <- expected_gep_length * save_last_n_mc
      expected_ssp_length <- cell_type_number
      expected_ssp_all_length <- expected_ssp_length * save_last_n_mc

      if (length(sampler_result$csGEP_vec) != expected_gep_length ||
          length(sampler_result$csGEP_all) != expected_gep_all_length ||
          length(sampler_result$SSP_vec) != expected_ssp_length ||
          length(sampler_result$SSP_all) != expected_ssp_all_length) {
        stop(
          "gibbsSampler returned an object with unexpected dimensions for sample ",
          sample_name,
          ".",
          call. = FALSE
        )
      }

      chain_sampler_gep[[ch]] <- t(matrix(
        sampler_result$csGEP_vec,
        nrow = cell_type_number,
        ncol = nrow(sample_data)
      ))
      chain_sampler_ssp[[ch]] <- matrix(
        sampler_result$SSP_vec,
        nrow = cell_type_number,
        ncol = 1L
      )

      chain_z[[ch]] <- array(
        sampler_result$csGEP_all,
        dim = c(cell_type_number, nrow(sample_data), save_last_n_mc),
        dimnames = list(cell_types, gene_symbols, NULL)
      )
      chain_ssp[[ch]] <- array(
        sampler_result$SSP_all,
        dim = c(cell_type_number, 1L, save_last_n_mc),
        dimnames = list(cell_types, sample_name, NULL)
      )
      last_sampler_result <- sampler_result
    }

    chain_id <- rep(seq_len(n_chains), each = save_last_n_mc)
    draw_within_chain <- rep(seq_len(save_last_n_mc), times = n_chains)
    draw_labels <- paste0(
      "chain", chain_id,
      "_draw", draw_within_chain
    )

    z_array <- .cdseq2_combine_3d_arrays(chain_z, draw_labels)
    ssp_array <- .cdseq2_combine_3d_arrays(chain_ssp, draw_labels)

    posterior_summary <- .cdseq2_summarize_sample(
      z_array = z_array,
      ssp_array = ssp_array,
      beta = beta,
      alpha = alpha_initial,
      cell_types = cell_types,
      gene_symbols = gene_symbols,
      return_posterior_sd = return_posterior_sd,
      return_phi_draw_means = return_phi_draw_means,
      low_information_depth = low_information_depth_vector,
      low_information_prior_weight = low_information_prior_weight_vector
    )

    expected_bulk_total <- sum(sample_data)
    maximum_total_mismatch <- max(abs(
      posterior_summary$internal$bulk_total_by_draw - expected_bulk_total
    ))
    if (posterior_summary$internal$maximum_assignment_mismatch > 1e-8 ||
        maximum_total_mismatch > 1e-8) {
      stop(
        "Saved MCMC arrays are not valid collapsed read-allocation states for sample ",
        sample_name,
        ". Check the semantics and storage order of csGEP_all and SSP_all.",
        call. = FALSE
      )
    }

    chain_diagnostics <- .cdseq2_chain_diagnostics(
      z_array = z_array,
      ssp_array = ssp_array,
      beta = beta,
      alpha = alpha_initial,
      chain_id = chain_id,
      cell_types = cell_types
    )

    # Preserve sampler-level legacy summaries for backward compatibility.
    sampler_estGEP_mat <- Reduce("+", chain_sampler_gep) / n_chains
    sampler_estSSp_mat <- Reduce("+", chain_sampler_ssp) / n_chains
    colnames(sampler_estGEP_mat) <- cell_types
    rownames(sampler_estGEP_mat) <- gene_symbols
    rownames(sampler_estSSp_mat) <- cell_types
    colnames(sampler_estSSp_mat) <- sample_name

    estGEP_legacy_sampler <- t(t(sampler_estGEP_mat + t(beta)) /
                                 colSums(sampler_estGEP_mat + t(beta)))
    estGEP2_sampler <- matrix(
      NA_real_,
      nrow = nrow(sampler_estGEP_mat),
      ncol = ncol(sampler_estGEP_mat),
      dimnames = dimnames(sampler_estGEP_mat)
    )
    positive_sampler_depth <- colSums(sampler_estGEP_mat) > 0
    if (any(positive_sampler_depth)) {
      estGEP2_sampler[, positive_sampler_depth] <- t(
        t(sampler_estGEP_mat[, positive_sampler_depth, drop = FALSE]) /
          colSums(sampler_estGEP_mat)[positive_sampler_depth]
      )
    }

    estProp_legacy_sampler <- t(
      t(sampler_estSSp_mat + alpha) /
        colSums(sampler_estSSp_mat + alpha)
    )
    estProp2_sampler <- t(
      t(sampler_estSSp_mat) /
        colSums(sampler_estSSp_mat)
    )

    recommended_estGEP <- posterior_summary$estGEP_RB
    estGEP <- if (primary_gep_estimator == "rao_blackwell") {
      recommended_estGEP
    } else {
      estGEP_legacy_sampler
    }
    estProp <- posterior_summary$estProp_RB
    colnames(estProp) <- sample_name
    if (!is.null(posterior_summary$estProp_RB_sd)) {
      colnames(posterior_summary$estProp_RB_sd) <- sample_name
    }
    colnames(posterior_summary$estProp_allocation) <- sample_name

    lgpst <- logpost(estProp, estGEP, sample_data, alpha, beta)
    lglike <- loglikelihood(estProp, estGEP, sample_data, alpha, beta)
    if (all(is.finite(estGEP2_sampler)) &&
        all(is.finite(estProp2_sampler))) {
      lgpst2 <- logpost(
        estProp2_sampler,
        estGEP2_sampler,
        sample_data,
        alpha,
        beta
      )
      lglike2 <- loglikelihood(
        estProp2_sampler,
        estGEP2_sampler,
        sample_data,
        alpha,
        beta
      )
    } else {
      lgpst2 <- NA_real_
      lglike2 <- NA_real_
    }

    aggregated_result <- NULL
    if (!is.null(state_to_major)) {
      aggregated_result <- .cdseq2_aggregate_states(
        z_array = z_array,
        ssp_array = ssp_array,
        beta = beta,
        alpha = alpha_initial,
        state_to_major = state_to_major,
        gene_symbols = gene_symbols,
        return_posterior_sd = return_posterior_sd,
        return_arrays = return_aggregated_arrays,
        low_information_depth = low_information_depth_major,
        low_information_prior_weight = low_information_prior_weight_major
      )
      colnames(aggregated_result$estProp) <- sample_name
      colnames(aggregated_result$estProp2) <- sample_name
      colnames(aggregated_result$estSSp_mat) <- sample_name
      if (!is.null(aggregated_result$estProp_sd)) {
        colnames(aggregated_result$estProp_sd) <- sample_name
      }
    }

    cell_type_assign_split <- if (return_assignments) {
      last_sampler_result$cellTypeAssignSplit
    } else {
      NULL
    }
    cell_type_bag <- if (return_assignments) {
      last_sampler_result$cellTypeBag
    } else {
      NULL
    }

    if (!is.null(cell_type_assign_split) &&
        length(dim(cell_type_assign_split)) == 3L) {
      dimnames(cell_type_assign_split) <- list(
        gene_symbols,
        sample_name,
        cell_types
      )
    }

    output <- list(
      # Backward-compatible primary names.
      estProp = estProp,
      estGEP = estGEP,
      logpost = lgpst,
      loglikelihood = lglike,
      estProp2 = estProp2_sampler,
      estGEP2 = estGEP2_sampler,
      logpost2 = lgpst2,
      loglikelihood2 = lglike2,

      # Recommended, explicitly named posterior summaries.
      estProp_RB = estProp,
      estProp_RB_sd = posterior_summary$estProp_RB_sd,
      estGEP_RB = posterior_summary$estGEP_RB,
      estGEP_RB_sd = posterior_summary$estGEP_RB_sd,

      # Alternative estimators retained for audits and sensitivity analyses.
      estProp_legacy = estProp_legacy_sampler,
      estGEP_legacy = estGEP_legacy_sampler,
      estProp_allocation = posterior_summary$estProp_allocation,
      estGEP_allocation = posterior_summary$estGEP_allocation,
      estGEP_draw_normalized_allocation =
        posterior_summary$estGEP_draw_normalized_allocation,

      # Unnormalized posterior mean allocation counts.
      estGEP_count = posterior_summary$allocation_mean,
      estGEP_mat = sampler_estGEP_mat,
      estSSp_mat = sampler_estSSp_mat,

      # Saved collapsed allocation states. estGEP_array is a legacy alias.
      estAllocation_array = if (return_mcmc_arrays) z_array else NULL,
      estGEP_array = if (return_mcmc_arrays) z_array else NULL,
      estSSP_array = if (return_mcmc_arrays) ssp_array else NULL,
      estPhi_conditional_mean_array =
        posterior_summary$phi_conditional_mean_array,
      draw_chain = chain_id,
      draw_within_chain = draw_within_chain,

      information = posterior_summary$information,
      mcmc_diagnostics = chain_diagnostics,
      allocation_checks = list(
        maximum_abs_difference_between_GEP_and_SSP_totals =
          posterior_summary$internal$maximum_assignment_mismatch,
        maximum_abs_difference_from_observed_bulk_total =
          maximum_total_mismatch,
        observed_bulk_total = expected_bulk_total,
        bulk_total_by_draw = posterior_summary$internal$bulk_total_by_draw
      ),
      aggregated = aggregated_result,

      cell_type_assignment = NULL,
      cellTypeAssignSplit = cell_type_assign_split,
      cellTypeBag = cell_type_bag,
      sample_names = sample_name,
      sample_running_time = Sys.time() - sample_start_time,
      primary_gep_estimator = primary_gep_estimator
    )

    output
  }

  if (cpu_number > 1L) {
    doParallel::stopImplicitCluster()
  }

  gibbs_end_time <- Sys.time()
  gibbs_running_time <- gibbs_end_time - gibbs_start_time
  names(allresult) <- sample_names
  cdseq_result <- lapply(allresult, function(x) {
    x$gibbsRunningTime <- gibbs_running_time
    x
  })

  parameters <- list(
    version = "cdseq2-performance-extraction-1.0",
    beta = beta,
    alpha = alpha_initial,
    reference_scRNA = reference_scRNA,
    reference_scRNA_original = reference_scRNA_original,
    reference_precision = reference_precision_vector,
    reference_proportion = reference_proportion_vector,
    reference_aggregation = reference_aggregation,
    reference_min_cells_per_donor_type = reference_min_cells_per_donor_type,
    reference_cell_count_by_type = reference_object$cell_count_by_type,
    reference_donor_count_by_type = reference_object$donor_count_by_type,
    reference_donor_group_table = reference_object$donor_group_table,
    reference_correlation = reference_correlation,
    smooth_ref = smooth_ref,
    smoothing_alpha = smoothing_alpha_vector,
    smoothing_scope = smoothing_scope,
    prior_adjustment = prior_adjustment,
    genes = gene_symbols,
    genes_use = genes_use,
    exclude_genes = exclude_genes,
    all_zero_genes_removed = all_zero_removed,
    cell_types = cell_types,
    state_to_major = state_to_major,
    major_cell_types = major_types,
    mcmc_iterations = mcmc_iterations,
    save_last_n_mc = save_last_n_mc,
    n_chains = n_chains,
    rng = rng,
    seed = seed,
    chain_seeds = seed_matrix,
    cpu_number = cpu_number,
    primary_gep_estimator = primary_gep_estimator,
    return_posterior_sd = return_posterior_sd
  )

  end_time <- Sys.time()
  cdseq2_running_time <- as.numeric(end_time - start_time, units = "hours")
  message(sprintf(
    "cdseq2 completed successfully using %.4f hours\n",
    cdseq2_running_time
  ))

  list(
    cdseq_result = cdseq_result,
    parameters = parameters,
    running_time = list(
      gibbs = gibbs_running_time,
      total_hours = cdseq2_running_time
    )
  )
}
