# CDSeq5: two-stage, near-fixed-proportion CDSeq inference
#
# This implementation deliberately reuses the existing CDSeq2 C++ Gibbs sampler.
# Stage 1 supplies a cell-type-specific expression prior (BETA) and, unless the
# user supplies proportions directly, a plug-in RNA-proportion matrix. Stage 2
# uses a high-concentration, sample-specific Dirichlet prior
#
#   ALPHA[, m] = kappa_m * theta_star[, m]
#   kappa_m     = lock_ratio_m * sum_g bulk_data[g, m]
#
# to make the composition term in the collapsed Gibbs update approach the fixed
# vector theta_star[, m]. This is near-fixed/proportion-anchored inference, not
# exact conditioning on theta unless the C++ transition kernel is modified.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.cdseq5_stop <- function(...) {
  stop(..., call. = FALSE)
}

.cdseq5_warn <- function(...) {
  warning(..., call. = FALSE, immediate. = TRUE)
}

.cdseq5_is_integer_like <- function(x, tolerance = 1e-8) {
  all(is.finite(x)) && all(abs(x - round(x)) <= tolerance)
}

.cdseq5_validate_bulk <- function(bulk_data) {
  if (inherits(bulk_data, "Matrix")) {
    bulk_data <- as.matrix(bulk_data)
  }
  if (!is.matrix(bulk_data) || length(dim(bulk_data)) != 2L) {
    .cdseq5_stop("bulk_data must be a two-dimensional gene-by-sample matrix.")
  }
  if (!is.numeric(bulk_data)) {
    .cdseq5_stop("bulk_data must be numeric.")
  }
  if (is.null(rownames(bulk_data)) || is.null(colnames(bulk_data))) {
    .cdseq5_stop("bulk_data must have gene row names and sample column names.")
  }
  if (anyNA(rownames(bulk_data)) || any(rownames(bulk_data) == "") ||
      anyDuplicated(rownames(bulk_data))) {
    .cdseq5_stop("bulk_data gene names must be complete and unique.")
  }
  if (anyNA(colnames(bulk_data)) || any(colnames(bulk_data) == "") ||
      anyDuplicated(colnames(bulk_data))) {
    .cdseq5_stop("bulk_data sample names must be complete and unique.")
  }
  if (anyNA(bulk_data) || any(!is.finite(bulk_data)) || any(bulk_data < 0)) {
    .cdseq5_stop("bulk_data must contain finite, nonnegative values.")
  }
  if (!.cdseq5_is_integer_like(bulk_data)) {
    .cdseq5_stop("bulk_data must contain untransformed integer-like read counts.")
  }
  if (any(colSums(bulk_data) <= 0)) {
    bad <- colnames(bulk_data)[colSums(bulk_data) <= 0]
    .cdseq5_stop(
      "Every bulk sample must have positive total counts. Invalid sample(s): ",
      paste(bad, collapse = ", ")
    )
  }
  storage.mode(bulk_data) <- "double"
  bulk_data
}

.cdseq5_validate_named_list <- function(x, label) {
  if (is.null(x)) return(list())
  if (!is.list(x)) .cdseq5_stop(label, " must be a list.")
  nm <- names(x)
  if (length(x) > 0L && (is.null(nm) || anyNA(nm) || any(nm == ""))) {
    .cdseq5_stop(label, " must contain only named elements.")
  }
  x
}

.cdseq5_merge_named_lists <- function(...) {
  inputs <- list(...)
  out <- list()
  for (input in inputs) {
    if (length(input) == 0L) next
    if (is.null(names(input)) || any(names(input) == "")) {
      .cdseq5_stop("All arguments forwarded to the Stage-1 function must be named.")
    }
    for (nm in names(input)) out[[nm]] <- input[[nm]]
  }
  out
}

.cdseq5_resolve_stage1_function <- function(stage1_function = NULL) {
  if (is.function(stage1_function)) return(stage1_function)

  if (is.character(stage1_function) && length(stage1_function) == 1L) {
    if (!exists(stage1_function, mode = "function", inherits = TRUE)) {
      .cdseq5_stop("Could not find Stage-1 function named '", stage1_function, "'.")
    }
    return(get(stage1_function, mode = "function", inherits = TRUE))
  }

  if (!is.null(stage1_function)) {
    .cdseq5_stop("stage1_function must be NULL, a function, or one function name.")
  }

  if (exists("cdseq4", mode = "function", inherits = TRUE)) {
    return(get("cdseq4", mode = "function", inherits = TRUE))
  }

  if (requireNamespace("CDSeq2", quietly = TRUE) &&
      "cdseq4" %in% getNamespaceExports("CDSeq2")) {
    return(getExportedValue("CDSeq2", "cdseq4"))
  }

  .cdseq5_stop(
    "No Stage-1 function was supplied and cdseq4() was not found. ",
    "Pass stage1_function = cdseq4, supply stage1_fit, or provide both ",
    "fixed_proportions and reference_beta."
  )
}

.cdseq5_call_stage1 <- function(stage1_function, args, strict = TRUE) {
  fn_formals <- names(formals(stage1_function))
  accepts_dots <- "..." %in% fn_formals
  unsupported <- setdiff(names(args), fn_formals)

  if (!accepts_dots && length(unsupported) > 0L) {
    if (isTRUE(strict)) {
      .cdseq5_stop(
        "The Stage-1 function does not accept argument(s): ",
        paste(unsupported, collapse = ", "),
        ". Its formals are: ", paste(fn_formals, collapse = ", "), "."
      )
    }
    .cdseq5_warn(
      "Dropping unsupported Stage-1 argument(s): ",
      paste(unsupported, collapse = ", ")
    )
    args <- args[intersect(names(args), fn_formals)]
  }

  do.call(stage1_function, args)
}

.cdseq5_get_path <- function(x, path) {
  current <- x
  for (nm in path) {
    if (!is.list(current) || is.null(current[[nm]])) return(NULL)
    current <- current[[nm]]
  }
  current
}

.cdseq5_find_first_path <- function(x, paths, predicate = function(z) !is.null(z)) {
  for (path in paths) {
    value <- .cdseq5_get_path(x, path)
    if (predicate(value)) {
      return(list(value = value, path = paste(path, collapse = "$")))
    }
  }
  NULL
}

.cdseq5_find_beta <- function(stage1_fit) {
  paths <- list(
    c("parameters", "beta"),
    c("parameters", "BETA"),
    c("reference_beta"),
    c("beta"),
    c("cdseq2_fit", "parameters", "beta"),
    c("cdseq2_fit", "parameters", "BETA"),
    c("final_fit", "parameters", "beta"),
    c("underlying_fit", "parameters", "beta"),
    c("fit", "parameters", "beta"),
    c("posterior_means", "reference_beta")
  )
  found <- .cdseq5_find_first_path(
    stage1_fit,
    paths,
    predicate = function(z) {
      (is.matrix(z) || inherits(z, "Matrix") || is.data.frame(z)) &&
        is.numeric(as.matrix(z))
    }
  )
  if (is.null(found)) {
    .cdseq5_stop(
      "Could not locate the Stage-1 BETA matrix. Supply reference_beta ",
      "explicitly or retain fit$parameters$beta in the Stage-1 object."
    )
  }
  found
}

.cdseq5_locate_result_list <- function(stage1_fit) {
  paths <- list(
    c("cdseq_result"),
    c("cdseq2_fit", "cdseq_result"),
    c("final_fit", "cdseq_result"),
    c("underlying_fit", "cdseq_result"),
    c("fit", "cdseq_result")
  )
  found <- .cdseq5_find_first_path(
    stage1_fit,
    paths,
    predicate = function(z) is.list(z) && length(z) > 0L
  )
  if (is.null(found)) {
    .cdseq5_stop("Could not locate a named cdseq_result list in the Stage-1 fit.")
  }
  found
}

.cdseq5_infer_cell_types <- function(
    cell_types,
    beta_raw,
    genes,
    reference_scRNA_annotation = NULL,
    fixed_proportions = NULL
) {
  if (!is.null(cell_types)) {
    out <- as.character(cell_types)
  } else if (!is.null(rownames(beta_raw)) && ncol(beta_raw) == length(genes) &&
             setequal(colnames(beta_raw) %||% genes, genes)) {
    out <- rownames(beta_raw)
  } else if (!is.null(colnames(beta_raw)) && nrow(beta_raw) == length(genes) &&
             setequal(rownames(beta_raw) %||% genes, genes)) {
    out <- colnames(beta_raw)
  } else if (!is.null(reference_scRNA_annotation)) {
    out <- unique(as.character(reference_scRNA_annotation))
  } else if (is.matrix(fixed_proportions) || is.data.frame(fixed_proportions)) {
    fixed_proportions <- as.matrix(fixed_proportions)
    if (!is.null(rownames(fixed_proportions))) {
      out <- rownames(fixed_proportions)
    } else if (!is.null(colnames(fixed_proportions))) {
      out <- colnames(fixed_proportions)
    } else {
      .cdseq5_stop("cell_types could not be inferred from unnamed fixed_proportions.")
    }
  } else {
    .cdseq5_stop(
      "cell_types could not be inferred. Supply cell_types explicitly or provide ",
      "named BETA/reference annotations."
    )
  }

  if (length(out) < 2L || anyNA(out) || any(out == "") || anyDuplicated(out)) {
    .cdseq5_stop("cell_types must contain at least two complete, unique labels.")
  }
  out
}

.cdseq5_orient_beta <- function(beta, cell_types, genes) {
  beta <- as.matrix(beta)
  if (!is.numeric(beta)) .cdseq5_stop("reference_beta/BETA must be numeric.")

  nr <- nrow(beta)
  nc <- ncol(beta)
  Tn <- length(cell_types)
  Gn <- length(genes)

  row_named_as_ct <- !is.null(rownames(beta)) && all(cell_types %in% rownames(beta))
  col_named_as_gene <- !is.null(colnames(beta)) && all(genes %in% colnames(beta))
  row_named_as_gene <- !is.null(rownames(beta)) && all(genes %in% rownames(beta))
  col_named_as_ct <- !is.null(colnames(beta)) && all(cell_types %in% colnames(beta))

  if (row_named_as_ct && col_named_as_gene) {
    beta <- beta[cell_types, genes, drop = FALSE]
  } else if (row_named_as_gene && col_named_as_ct) {
    beta <- t(beta[genes, cell_types, drop = FALSE])
  } else if (nr == Tn && nc == Gn && !(nr == Gn && nc == Tn)) {
    .cdseq5_warn(
      "BETA lacked sufficient dimnames; assuming rows are cell types and columns are genes."
    )
    dimnames(beta) <- list(cell_types, genes)
  } else if (nr == Gn && nc == Tn && !(nr == Tn && nc == Gn)) {
    .cdseq5_warn(
      "BETA lacked sufficient dimnames; assuming rows are genes and columns are cell types."
    )
    beta <- t(beta)
    dimnames(beta) <- list(cell_types, genes)
  } else {
    .cdseq5_stop(
      "Could not orient BETA as cell type x gene. Observed dimension ",
      nr, " x ", nc, "; expected ", Tn, " x ", Gn, " or its transpose."
    )
  }

  if (anyNA(beta) || any(!is.finite(beta)) || any(beta <= 0)) {
    .cdseq5_stop("Every BETA entry must be finite and strictly positive.")
  }
  storage.mode(beta) <- "double"
  beta
}

.cdseq5_standardize_proportion_matrix <- function(x, cell_types, sample_names, label) {
  if (is.data.frame(x)) x <- as.matrix(x)

  if (is.list(x) && !is.matrix(x)) {
    if (is.null(names(x)) || !setequal(names(x), sample_names)) {
      .cdseq5_stop(label, " list must be named by all bulk samples.")
    }
    out <- matrix(
      NA_real_,
      nrow = length(cell_types),
      ncol = length(sample_names),
      dimnames = list(cell_types, sample_names)
    )
    for (sample_name in sample_names) {
      v <- x[[sample_name]]
      if (is.data.frame(v)) v <- as.matrix(v)
      if (is.matrix(v)) {
        if (ncol(v) == 1L && !is.null(rownames(v)) &&
            all(cell_types %in% rownames(v))) {
          v <- v[cell_types, 1L]
        } else if (nrow(v) == 1L && !is.null(colnames(v)) &&
                   all(cell_types %in% colnames(v))) {
          v <- v[1L, cell_types]
        } else if (length(v) == length(cell_types)) {
          v <- as.numeric(v)
        } else {
          .cdseq5_stop(label, " entry for ", sample_name, " has the wrong shape.")
        }
      }
      if (!is.null(names(v)) && all(cell_types %in% names(v))) {
        v <- v[cell_types]
      } else if (length(v) != length(cell_types)) {
        .cdseq5_stop(label, " entry for ", sample_name, " has the wrong length.")
      }
      out[, sample_name] <- as.numeric(v)
    }
    x <- out
  }

  if (is.atomic(x) && is.null(dim(x))) {
    if (length(sample_names) != 1L) {
      .cdseq5_stop(label, " may be a vector only for a single bulk sample.")
    }
    if (!is.null(names(x)) && all(cell_types %in% names(x))) x <- x[cell_types]
    if (length(x) != length(cell_types)) {
      .cdseq5_stop(label, " vector has the wrong length.")
    }
    x <- matrix(
      as.numeric(x),
      nrow = length(cell_types),
      ncol = 1L,
      dimnames = list(cell_types, sample_names)
    )
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    .cdseq5_stop(label, " must be a numeric matrix, data frame, named list, or vector.")
  }

  row_is_ct <- !is.null(rownames(x)) && all(cell_types %in% rownames(x))
  col_is_sample <- !is.null(colnames(x)) && all(sample_names %in% colnames(x))
  row_is_sample <- !is.null(rownames(x)) && all(sample_names %in% rownames(x))
  col_is_ct <- !is.null(colnames(x)) && all(cell_types %in% colnames(x))

  if (row_is_ct && col_is_sample) {
    out <- x[cell_types, sample_names, drop = FALSE]
  } else if (row_is_sample && col_is_ct) {
    out <- t(x[sample_names, cell_types, drop = FALSE])
  } else if (identical(dim(x), c(length(cell_types), length(sample_names)))) {
    .cdseq5_warn(
      label, " lacked sufficient dimnames; assuming cell type x sample orientation."
    )
    out <- x
    dimnames(out) <- list(cell_types, sample_names)
  } else if (identical(dim(x), c(length(sample_names), length(cell_types)))) {
    .cdseq5_warn(
      label, " lacked sufficient dimnames; assuming sample x cell type orientation."
    )
    out <- t(x)
    dimnames(out) <- list(cell_types, sample_names)
  } else {
    .cdseq5_stop(
      "Could not orient ", label, ". Observed dimension ",
      paste(dim(x), collapse = " x "), "; expected cell type x sample or its transpose."
    )
  }

  if (anyNA(out) || any(!is.finite(out)) || any(out < 0)) {
    .cdseq5_stop(label, " must contain finite, nonnegative values.")
  }
  if (any(colSums(out) <= 0)) {
    bad <- colnames(out)[colSums(out) <= 0]
    .cdseq5_stop(label, " has nonpositive column sums for: ", paste(bad, collapse = ", "))
  }
  storage.mode(out) <- "double"
  out
}

.cdseq5_extract_stage1_proportions <- function(
    stage1_fit,
    cell_types,
    sample_names,
    component = "estProp"
) {
  # Prefer the explicitly requested sample-level component whenever a
  # cdseq_result list is available. High-level rna_proportions is a fallback for
  # wrappers that no longer retain the sample-level objects.
  result_paths <- list(
    c("cdseq_result"),
    c("cdseq2_fit", "cdseq_result"),
    c("final_fit", "cdseq_result"),
    c("underlying_fit", "cdseq_result"),
    c("fit", "cdseq_result")
  )
  located <- .cdseq5_find_first_path(
    stage1_fit,
    result_paths,
    predicate = function(z) is.list(z) && length(z) > 0L
  )

  if (!is.null(located)) {
    result_list <- located$value
    if (is.null(names(result_list)) || !setequal(names(result_list), sample_names)) {
      .cdseq5_stop("Stage-1 cdseq_result must be named by all bulk samples.")
    }
    result_list <- result_list[sample_names]

    out <- matrix(
      NA_real_,
      nrow = length(cell_types),
      ncol = length(sample_names),
      dimnames = list(cell_types, sample_names)
    )

    for (sample_name in sample_names) {
      sample_result <- result_list[[sample_name]]
      value <- sample_result[[component]]
      if (is.null(value)) {
        .cdseq5_stop(
          "Stage-1 sample '", sample_name, "' does not contain proportion component '",
          component, "'. Available components: ",
          paste(names(sample_result), collapse = ", "), "."
        )
      }
      value <- as.matrix(value)
      if (ncol(value) == 1L && !is.null(rownames(value)) &&
          all(cell_types %in% rownames(value))) {
        out[, sample_name] <- value[cell_types, 1L]
      } else if (nrow(value) == 1L && !is.null(colnames(value)) &&
                 all(cell_types %in% colnames(value))) {
        out[, sample_name] <- value[1L, cell_types]
      } else if (length(value) == length(cell_types)) {
        v <- as.numeric(value)
        candidate_names <- if (ncol(value) == 1L) {
          rownames(value)
        } else if (nrow(value) == 1L) {
          colnames(value)
        } else {
          NULL
        }
        if (!is.null(candidate_names)) names(v) <- candidate_names
        if (!is.null(names(v)) && all(cell_types %in% names(v))) v <- v[cell_types]
        out[, sample_name] <- v
      } else {
        .cdseq5_stop(
          "Could not interpret Stage-1 proportion component for sample '",
          sample_name, "'."
        )
      }
    }

    return(list(
      value = out,
      path = paste0(located$path, "[[sample]]$", component)
    ))
  }

  direct_paths <- list(
    c("rna_proportions"),
    c("posterior_means", "rna_proportions")
  )
  direct <- .cdseq5_find_first_path(
    stage1_fit,
    direct_paths,
    predicate = function(z) is.matrix(z) || is.data.frame(z)
  )
  if (!is.null(direct)) {
    return(list(
      value = .cdseq5_standardize_proportion_matrix(
        direct$value, cell_types, sample_names, "Stage-1 proportions"
      ),
      path = direct$path
    ))
  }

  .cdseq5_stop(
    "Could not locate Stage-1 proportions in either cdseq_result or a high-level ",
    "rna_proportions component."
  )
}

.cdseq5_extract_stage1_assignments <- function(stage1_fit, sample_names) {
  located <- .cdseq5_locate_result_list(stage1_fit)
  result_list <- located$value
  if (is.null(names(result_list)) || !setequal(names(result_list), sample_names)) {
    .cdseq5_stop("Stage-1 cdseq_result must be named by all bulk samples.")
  }
  result_list <- result_list[sample_names]

  out <- lapply(sample_names, function(sample_name) {
    value <- result_list[[sample_name]]$cellTypeBag
    if (is.null(value)) {
      .cdseq5_stop(
        "initialization = 'stage1' requested, but cellTypeBag is absent for sample '",
        sample_name, "'."
      )
    }
    value <- as.integer(value)
    if (length(value) == 0L || anyNA(value) || any(value < 0L)) {
      .cdseq5_stop("Invalid Stage-1 cellTypeBag for sample '", sample_name, "'.")
    }
    value
  })
  names(out) <- sample_names
  out
}

.cdseq5_prepare_fixed_proportions <- function(
    x,
    cell_types,
    sample_names,
    proportion_floor
) {
  raw <- .cdseq5_standardize_proportion_matrix(
    x, cell_types, sample_names, "fixed_proportions"
  )
  normalized_raw <- sweep(raw, 2L, colSums(raw), "/")
  adjusted <- pmax(normalized_raw, proportion_floor)
  adjusted <- sweep(adjusted, 2L, colSums(adjusted), "/")

  list(
    raw = raw,
    normalized_raw = normalized_raw,
    adjusted = adjusted,
    max_abs_floor_change = max(abs(adjusted - normalized_raw))
  )
}

.cdseq5_expand_sample_parameter <- function(x, sample_names, label) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
    .cdseq5_stop(label, " must contain finite positive values.")
  }
  if (length(x) == 1L) {
    out <- rep(as.numeric(x), length(sample_names))
    names(out) <- sample_names
    return(out)
  }
  if (!is.null(names(x))) {
    if (!all(sample_names %in% names(x))) {
      .cdseq5_stop(label, " is missing sample(s): ",
                   paste(setdiff(sample_names, names(x)), collapse = ", "))
    }
    out <- as.numeric(x[sample_names])
    names(out) <- sample_names
    return(out)
  }
  if (length(x) == length(sample_names)) {
    out <- as.numeric(x)
    names(out) <- sample_names
    return(out)
  }
  .cdseq5_stop(label, " must be scalar, sample-named, or have one value per sample.")
}

.cdseq5_resolve_gibbs_sampler <- function(gibbs_sampler = NULL) {
  if (is.function(gibbs_sampler)) return(gibbs_sampler)
  if (!is.null(gibbs_sampler)) {
    .cdseq5_stop("gibbs_sampler must be NULL or a function.")
  }
  if (exists("gibbsSampler", mode = "function", inherits = TRUE)) {
    return(get("gibbsSampler", mode = "function", inherits = TRUE))
  }
  if (!requireNamespace("CDSeq2", quietly = TRUE)) {
    .cdseq5_stop("CDSeq2 must be installed to resolve its compiled Gibbs sampler.")
  }
  tryCatch(
    utils::getFromNamespace("gibbsSampler", "CDSeq2"),
    error = function(e) {
      .cdseq5_stop(
        "Could not resolve CDSeq2's internal gibbsSampler(). Original error: ",
        conditionMessage(e)
      )
    }
  )
}

.cdseq5_safe_column_normalize <- function(x, label) {
  totals <- colSums(x)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    .cdseq5_stop(label, " contains a nonpositive column total.")
  }
  sweep(x, 2L, totals, "/")
}

.cdseq5_run_one_sample <- function(
    sample_index,
    bulk_data,
    beta,
    fixed_proportions,
    fixed_proportions_before_floor,
    kappa_by_sample,
    lock_ratio_by_sample,
    mcmc_iterations,
    save_last_n_mc,
    output_control,
    initialization,
    stage1_assignments,
    verbose_int,
    rng,
    seed,
    return_mcmc_arrays,
    return_assignments,
    gibbs_sampler
) {
  sample_name <- colnames(bulk_data)[[sample_index]]
  genes <- rownames(bulk_data)
  cell_types <- rownames(beta)
  Tn <- length(cell_types)
  Gn <- length(genes)

  sample_counts <- bulk_data[, sample_index, drop = FALSE]
  colnames(sample_counts) <- sample_name
  sample_depth <- sum(sample_counts)
  theta_star <- fixed_proportions[, sample_name]
  theta_before_floor <- fixed_proportions_before_floor[, sample_name]
  kappa <- kappa_by_sample[[sample_name]]
  alpha <- matrix(
    kappa * theta_star,
    nrow = Tn,
    ncol = 1L,
    dimnames = list(cell_types, sample_name)
  )
  if (anyNA(alpha) || any(!is.finite(alpha)) || any(alpha <= 0)) {
    .cdseq5_stop(
      "Invalid Stage-2 ALPHA values for sample '", sample_name,
      "'. Reduce lock_ratio/lock_concentration or inspect fixed proportions."
    )
  }

  initialization_flag <- if (identical(initialization, "stage1")) 1L else 0L
  ct_initial <- if (initialization_flag == 1L) {
    stage1_assignments[[sample_name]]
  } else {
    as.integer(0L)
  }

  withr::local_seed(as.integer(seed + sample_index - 1L))
  sampler <- .cdseq5_resolve_gibbs_sampler(gibbs_sampler)

  timing <- system.time({
    sampler_result <- sampler(
      alpha,
      beta,
      sample_counts,
      as.integer(Tn),
      as.integer(mcmc_iterations),
      as.integer(output_control),
      as.integer(initialization_flag),
      ct_initial,
      as.integer(verbose_int),
      as.integer(rng),
      as.integer(save_last_n_mc)
    )
  })

  required <- c("csGEP_vec", "SSP_vec", "csGEP_all", "SSP_all")
  missing <- setdiff(required, names(sampler_result))
  if (length(missing) > 0L) {
    .cdseq5_stop(
      "gibbsSampler output for sample '", sample_name,
      "' is missing: ", paste(missing, collapse = ", "), "."
    )
  }

  if (length(sampler_result$csGEP_vec) != Tn * Gn) {
    .cdseq5_stop("Unexpected csGEP_vec length for sample '", sample_name, "'.")
  }
  if (length(sampler_result$SSP_vec) != Tn) {
    .cdseq5_stop("Unexpected SSP_vec length for sample '", sample_name, "'.")
  }
  if (length(sampler_result$csGEP_all) != Tn * Gn * save_last_n_mc) {
    .cdseq5_stop("Unexpected csGEP_all length for sample '", sample_name, "'.")
  }
  if (length(sampler_result$SSP_all) != Tn * save_last_n_mc) {
    .cdseq5_stop("Unexpected SSP_all length for sample '", sample_name, "'.")
  }

  numeric_outputs <- sampler_result[required]
  invalid_output <- vapply(
    numeric_outputs,
    function(x) {
      !is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x < 0)
    },
    logical(1)
  )
  if (any(invalid_output)) {
    .cdseq5_stop(
      "gibbsSampler returned nonnumeric, nonfinite, or negative values for sample '",
      sample_name, "' in: ",
      paste(names(invalid_output)[invalid_output], collapse = ", "), "."
    )
  }

  estGEP_mat <- t(matrix(
    as.numeric(sampler_result$csGEP_vec),
    nrow = Tn,
    ncol = Gn
  ))
  dimnames(estGEP_mat) <- list(genes, cell_types)

  estSSp_mat <- matrix(
    as.numeric(sampler_result$SSP_vec),
    nrow = Tn,
    ncol = 1L,
    dimnames = list(cell_types, sample_name)
  )

  estGEP_array <- array(
    as.numeric(sampler_result$csGEP_all),
    dim = c(Tn, Gn, save_last_n_mc),
    dimnames = list(
      cell_type = cell_types,
      gene = genes,
      draw = as.character(seq_len(save_last_n_mc))
    )
  )
  estSSP_array <- array(
    as.numeric(sampler_result$SSP_all),
    dim = c(Tn, 1L, save_last_n_mc),
    dimnames = list(
      cell_type = cell_types,
      sample = sample_name,
      draw = as.character(seq_len(save_last_n_mc))
    )
  )

  count_draw_mean_tg <- apply(estGEP_array, c(1L, 2L), mean)
  count_draw_mean <- t(count_draw_mean_tg)
  dimnames(count_draw_mean) <- list(genes, cell_types)
  ctsGEP_normalized_posthoc <- .cdseq5_safe_column_normalize(
    count_draw_mean,
    paste0("posterior mean allocated counts for ", sample_name)
  )

  phi_draw_sum <- matrix(
    0,
    nrow = Tn,
    ncol = Gn,
    dimnames = list(cell_types, genes)
  )
  valid_phi_draws <- setNames(integer(Tn), cell_types)
  zero_profile_draws <- setNames(integer(Tn), cell_types)
  profile_total_draws <- matrix(
    NA_real_,
    nrow = Tn,
    ncol = save_last_n_mc,
    dimnames = list(cell_types, as.character(seq_len(save_last_n_mc)))
  )
  effective_draws <- profile_total_draws
  allocation_draws <- profile_total_draws

  for (draw_index in seq_len(save_last_n_mc)) {
    draw_counts <- estGEP_array[, , draw_index, drop = FALSE][, , 1L]
    totals <- rowSums(draw_counts)
    positive_profiles <- is.finite(totals) & totals > 0
    if (any(positive_profiles)) {
      phi_draw_sum[positive_profiles, ] <-
        phi_draw_sum[positive_profiles, , drop = FALSE] +
        sweep(
          draw_counts[positive_profiles, , drop = FALSE],
          1L,
          totals[positive_profiles],
          "/"
        )
      valid_phi_draws[positive_profiles] <-
        valid_phi_draws[positive_profiles] + 1L
    }
    zero_profile_draws[!positive_profiles] <-
      zero_profile_draws[!positive_profiles] + 1L
    profile_total_draws[, draw_index] <- totals

    assignment_counts <- estSSP_array[, 1L, draw_index]
    assignment_total <- sum(assignment_counts)
    if (assignment_total <= 0) {
      .cdseq5_stop("A retained SSP draw has zero total for sample '", sample_name, "'.")
    }
    allocation_draws[, draw_index] <- assignment_counts / assignment_total
    effective_draws[, draw_index] <-
      (assignment_counts + alpha[, 1L]) /
      (assignment_total + sum(alpha[, 1L]))
  }

  ctsGEP_normalized_draw_mean_tg <- matrix(
    NA_real_,
    nrow = Tn,
    ncol = Gn,
    dimnames = list(cell_types, genes)
  )
  valid_cell_types <- valid_phi_draws > 0L
  ctsGEP_normalized_draw_mean_tg[valid_cell_types, ] <-
    sweep(
      phi_draw_sum[valid_cell_types, , drop = FALSE],
      1L,
      valid_phi_draws[valid_cell_types],
      "/"
    )
  ctsGEP_normalized_draw_mean <- t(ctsGEP_normalized_draw_mean_tg)
  dimnames(ctsGEP_normalized_draw_mean) <- list(genes, cell_types)

  allocation_mass_tolerance <- max(1e-6, 1e-10 * sample_depth)
  if (abs(sum(estSSp_mat) - sample_depth) > allocation_mass_tolerance) {
    .cdseq5_stop(
      "SSP_vec does not preserve the bulk count total for sample '",
      sample_name, "'."
    )
  }
  if (abs(sum(estGEP_mat) - sample_depth) > allocation_mass_tolerance) {
    .cdseq5_stop(
      "csGEP_vec does not preserve the bulk count total for sample '",
      sample_name, "'."
    )
  }
  ssp_draw_totals <- colSums(matrix(
    as.numeric(estSSP_array),
    nrow = Tn,
    ncol = save_last_n_mc
  ))
  gep_draw_totals <- vapply(
    seq_len(save_last_n_mc),
    function(draw_index) sum(estGEP_array[, , draw_index]),
    numeric(1)
  )
  if (any(abs(ssp_draw_totals - sample_depth) > allocation_mass_tolerance) ||
      any(abs(gep_draw_totals - sample_depth) > allocation_mass_tolerance)) {
    .cdseq5_stop(
      "One or more retained Gibbs draws do not preserve the bulk count total for sample '",
      sample_name, "'."
    )
  }

  posterior_gep <- .cdseq5_safe_column_normalize(
    estGEP_mat + t(beta),
    paste0("BETA-adjusted GEP for ", sample_name)
  )
  allocation_gep <- .cdseq5_safe_column_normalize(
    estGEP_mat,
    paste0("allocation GEP for ", sample_name)
  )

  allocation_proportion <- estSSp_mat / sum(estSSp_mat)
  effective_proportion <- (estSSp_mat + alpha) / sum(estSSp_mat + alpha)
  fixed_proportion <- matrix(
    theta_star,
    nrow = Tn,
    ncol = 1L,
    dimnames = list(cell_types, sample_name)
  )

  effective_abs_dev <- abs(effective_draws - theta_star)
  allocation_abs_dev <- abs(allocation_draws - theta_star)

  lock_diagnostic <- data.frame(
    sample = sample_name,
    sample_depth = sample_depth,
    lock_ratio = lock_ratio_by_sample[[sample_name]],
    kappa = kappa,
    theoretical_assignment_feedback_weight = sample_depth / (sample_depth + kappa),
    max_abs_effective_draw_deviation = max(effective_abs_dev),
    median_abs_effective_draw_deviation = stats::median(effective_abs_dev),
    p95_abs_effective_draw_deviation = unname(
      stats::quantile(effective_abs_dev, 0.95, names = FALSE)
    ),
    max_abs_effective_mean_deviation = max(
      abs(rowMeans(effective_draws) - theta_star)
    ),
    max_abs_allocation_mean_deviation = max(
      abs(rowMeans(allocation_draws) - theta_star)
    ),
    floor_adjustment_applied = max(abs(theta_star - theta_before_floor)) > 0,
    cell_types_with_any_zero_profile_draw = sum(zero_profile_draws > 0L),
    maximum_zero_profile_draw_fraction = max(
      zero_profile_draws / save_last_n_mc
    ),
    elapsed_seconds = unname(timing[["elapsed"]]),
    stringsAsFactors = FALSE
  )

  output <- list(
    sample_name = sample_name,
    fixed_proportion = fixed_proportion,
    allocation_proportion = allocation_proportion,
    effective_proportion = effective_proportion,
    estProp = effective_proportion,
    estProp2 = allocation_proportion,
    estGEP = posterior_gep,
    estGEP2 = allocation_gep,
    estGEP_mat = estGEP_mat,
    estSSp_mat = estSSp_mat,
    ctsGEP_counts = count_draw_mean,
    ctsGEP_normalized = ctsGEP_normalized_posthoc,
    ctsGEP_normalized_posthoc = ctsGEP_normalized_posthoc,
    ctsGEP_normalized_draw_mean = ctsGEP_normalized_draw_mean,
    profile_total_posterior_mean = rowMeans(profile_total_draws),
    profile_total_posterior_sd = apply(profile_total_draws, 1L, stats::sd),
    effective_proportion_draw_mean = rowMeans(effective_draws),
    allocation_proportion_draw_mean = rowMeans(allocation_draws),
    alpha = alpha,
    kappa = kappa,
    lock_ratio = lock_ratio_by_sample[[sample_name]],
    valid_normalized_draws_by_cell_type = valid_phi_draws,
    zero_profile_draws_by_cell_type = zero_profile_draws,
    lock_diagnostic = lock_diagnostic,
    gibbs_running_time_seconds = unname(timing[["elapsed"]]),
    gibbsRunningTime = unname(timing[["elapsed"]])
  )

  if (isTRUE(return_mcmc_arrays)) {
    output$estGEP_array <- estGEP_array
    output$estSSP_array <- estSSP_array
    output$effective_proportion_draws <- effective_draws
    output$allocation_proportion_draws <- allocation_draws
    output$profile_total_draws <- profile_total_draws
  }

  if (isTRUE(return_assignments)) {
    output$cellTypeAssignSplit <- if (!is.null(sampler_result$cellTypeAssignSplit)) {
      sampler_result$cellTypeAssignSplit
    } else {
      NULL
    }
    output$cellTypeBag <- if (!is.null(sampler_result$cellTypeBag)) {
      sampler_result$cellTypeBag
    } else {
      NULL
    }
  }

  output
}

#' Two-stage CDSeq with strongly anchored plug-in RNA proportions
#'
#' `cdseq5()` first obtains (or receives) a standard CDSeq/CDSeq4 fit and then
#' reruns the existing compiled Gibbs sampler using a high-concentration,
#' sample-specific Dirichlet prior centered on plug-in RNA proportions. The
#' second stage does not require a change to the C++ Gibbs sampler.
#'
#' The Stage-2 hyperparameter for sample `m` is
#'
#' `ALPHA[, m] = kappa_m * theta_star[, m]`,
#'
#' with `kappa_m = lock_ratio * sum(bulk_data[, m])` by default. The composition
#' term is therefore approximately 99.01% anchored to `theta_star` when
#' `lock_ratio = 100`. This is near-fixed rather than exact fixed-proportion
#' inference because the finite-ALPHA Gibbs kernel still contains assignment
#' counts.
#'
#' @param bulk_data Gene-by-sample matrix of untransformed, nonnegative,
#'   integer-like bulk RNA-seq counts.
#' @param mcmc_iterations Number of Stage-2 Gibbs iterations.
#' @param cpu_number Number of sample-level parallel workers.
#' @param reference_scRNA_count Optional reference count matrix passed to the
#'   Stage-1 function when `stage1_fit` is not supplied.
#' @param reference_scRNA_annotation Optional reference cell-type labels passed
#'   to Stage 1 and used to infer cell-type order when needed.
#' @param reference_metadata Optional reference metadata passed to Stage 1.
#' @param reference_scRNA_donor Optional reference-donor vector passed to Stage 1.
#' @param reference_proportion Stage-1 proportion prior.
#' @param reference_precision Stage-1 expression-prior precision or global
#'   precision center. CDSeq5 reuses the resulting Stage-1 BETA matrix exactly.
#' @param smooth_ref,smoothing_alpha Standard Stage-1 reference smoothing options.
#' @param rng Compiled Gibbs RNG implementation: 1 or 2.
#' @param save_last_n_mc Number of retained Stage-2 draws.
#' @param n_chains Currently must equal 1 because the existing C++ function is
#'   called directly once per sample.
#' @param seed R-level seed used before Stage 1 and to initialize sample workers.
#' @param stage1_function Stage-1 function object or function name. Defaults to
#'   `cdseq4()` if it is available.
#' @param stage1_fit Existing Stage-1 fit. Supplying this is recommended for a
#'   controlled benchmark because the identical fit and BETA are reused.
#' @param stage1_extra_args Named list of additional Stage-1 arguments.
#' @param strict_stage1_args Stop, rather than silently drop, arguments not
#'   supported by a Stage-1 function without `...`.
#' @param fixed_proportions Optional plug-in RNA proportions. Supply a cell-type
#'   by sample matrix, its transpose, or a sample-named list. If NULL, proportions
#'   are extracted from `stage1_fit` using `proportion_component`.
#' @param proportion_component Stage-1 sample-level component, normally `estProp`.
#' @param cell_types Optional explicit cell-type order.
#' @param reference_beta Optional explicit Stage-1 BETA matrix. If NULL it is
#'   extracted from `stage1_fit`.
#' @param lock_ratio Positive scalar or sample-specific vector. With default
#'   depth-relative locking, `kappa_m = lock_ratio_m * N_m`.
#' @param lock_concentration Optional absolute `kappa_m`. When supplied, it
#'   overrides `lock_ratio`; the realized ratio is reported as `kappa_m / N_m`.
#' @param proportion_floor Positive floor applied before renormalizing plug-in
#'   proportions, required because Dirichlet parameters must be positive.
#' @param initialization Either `"random"` or `"stage1"`. Stage-1 initialization
#'   is only a warm start; it is not the locking mechanism.
#' @param output_control Integer passed to the compiled Gibbs sampler.
#' @param keep_stage1_fit Retain the Stage-1 fit in the returned object.
#' @param return_mcmc_arrays Retain ctsGEP/SSP arrays and lock-audit draws.
#' @param return_assignments Retain compiled assignment outputs when available.
#' @param verbose Print progress messages.
#' @param gibbs_sampler Optional sampler function for testing. By default,
#'   `CDSeq2`'s internal `gibbsSampler()` is resolved at runtime.
#' @param ... Additional named arguments forwarded to Stage 1.
#'
#' @return An object of class `cdseq5_fit`. Main components include
#'   `ctsGEP`, `ctsGEP_normalized`, `fixed_proportions`,
#'   `allocation_proportions`, `effective_proportions`, `lock_diagnostics`,
#'   and the sample-level `cdseq_result` list.
#'
#' @export
cdseq5 <- function(
    bulk_data,
    mcmc_iterations = 700L,
    cpu_number = 1L,
    reference_scRNA_count = NULL,
    reference_scRNA_annotation = NULL,
    reference_metadata = NULL,
    reference_scRNA_donor = NULL,
    reference_proportion = 1,
    reference_precision = 1e6,
    smooth_ref = TRUE,
    smoothing_alpha = 1,
    rng = 1L,
    save_last_n_mc = 100L,
    n_chains = 1L,
    seed = 123L,
    stage1_function = NULL,
    stage1_fit = NULL,
    stage1_extra_args = list(),
    strict_stage1_args = TRUE,
    fixed_proportions = NULL,
    proportion_component = "estProp",
    cell_types = NULL,
    reference_beta = NULL,
    lock_ratio = 100,
    lock_concentration = NULL,
    proportion_floor = 1e-8,
    initialization = c("random", "stage1"),
    output_control = 1L,
    keep_stage1_fit = TRUE,
    return_mcmc_arrays = TRUE,
    return_assignments = TRUE,
    verbose = FALSE,
    gibbs_sampler = NULL,
    ...
) {
  call <- match.call()
  start_time <- Sys.time()
  initialization <- match.arg(initialization)
  bulk_data <- .cdseq5_validate_bulk(bulk_data)
  sample_names <- colnames(bulk_data)
  genes <- rownames(bulk_data)

  if (!is.numeric(mcmc_iterations) || length(mcmc_iterations) != 1L ||
      !is.finite(mcmc_iterations) || mcmc_iterations < 2 ||
      abs(mcmc_iterations - round(mcmc_iterations)) > 1e-8) {
    .cdseq5_stop("mcmc_iterations must be an integer greater than one.")
  }
  mcmc_iterations <- as.integer(round(mcmc_iterations))
  if (!is.numeric(save_last_n_mc) || length(save_last_n_mc) != 1L ||
      !is.finite(save_last_n_mc) || save_last_n_mc < 1 ||
      save_last_n_mc > mcmc_iterations ||
      abs(save_last_n_mc - round(save_last_n_mc)) > 1e-8) {
    .cdseq5_stop(
      "save_last_n_mc must be an integer between 1 and mcmc_iterations."
    )
  }
  save_last_n_mc <- as.integer(round(save_last_n_mc))
  if (!is.numeric(cpu_number) || length(cpu_number) != 1L ||
      !is.finite(cpu_number) || cpu_number < 1) {
    .cdseq5_stop("cpu_number must be a positive integer.")
  }
  cpu_number <- as.integer(ceiling(cpu_number))
  if (!is.numeric(n_chains) || length(n_chains) != 1L ||
      !is.finite(n_chains) || abs(n_chains - 1) > 1e-8) {
    .cdseq5_stop(
      "This direct-C++ CDSeq5 implementation currently supports n_chains = 1 only."
    )
  }
  if (!is.numeric(rng) || length(rng) != 1L || is.na(rng) ||
      !is.finite(rng) || !rng %in% c(1, 2)) {
    .cdseq5_stop("rng must be the scalar 1 or 2.")
  }
  rng <- as.integer(rng)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || abs(seed - round(seed)) > 1e-8) {
    .cdseq5_stop("seed must be one finite integer.")
  }
  seed <- as.integer(round(seed))
  if (!is.numeric(output_control) || length(output_control) != 1L ||
      is.na(output_control) || !is.finite(output_control) ||
      abs(output_control - round(output_control)) > 1e-8) {
    .cdseq5_stop("output_control must be one finite integer.")
  }
  output_control <- as.integer(round(output_control))
  if (!is.numeric(proportion_floor) || length(proportion_floor) != 1L ||
      !is.finite(proportion_floor) || proportion_floor <= 0 ||
      proportion_floor >= 1) {
    .cdseq5_stop("proportion_floor must be a finite value strictly between 0 and 1.")
  }
  if (!is.character(proportion_component) || length(proportion_component) != 1L ||
      is.na(proportion_component) || proportion_component == "") {
    .cdseq5_stop("proportion_component must be one nonempty character string.")
  }

  stage1_extra_args <- .cdseq5_validate_named_list(
    stage1_extra_args, "stage1_extra_args"
  )
  dots <- .cdseq5_validate_named_list(list(...), "Additional Stage-1 arguments")

  need_stage1 <- is.null(stage1_fit) &&
    (is.null(reference_beta) || is.null(fixed_proportions) ||
       identical(initialization, "stage1"))

  if (need_stage1) {
    stage1_function <- .cdseq5_resolve_stage1_function(stage1_function)
    common_stage1_args <- list(
      bulk_data = bulk_data,
      mcmc_iterations = mcmc_iterations,
      cpu_number = cpu_number,
      reference_scRNA_count = reference_scRNA_count,
      reference_scRNA_annotation = reference_scRNA_annotation,
      reference_metadata = reference_metadata,
      reference_scRNA_donor = reference_scRNA_donor,
      reference_proportion = reference_proportion,
      reference_precision = reference_precision,
      smooth_ref = smooth_ref,
      smoothing_alpha = smoothing_alpha,
      rng = rng,
      save_last_n_mc = save_last_n_mc,
      n_chains = 1L,
      seed = seed,
      return_mcmc_arrays = TRUE,
      return_assignments = identical(initialization, "stage1") || return_assignments,
      verbose = verbose
    )
    common_stage1_args <- common_stage1_args[
      !vapply(common_stage1_args, is.null, logical(1))
    ]
    stage1_args <- .cdseq5_merge_named_lists(
      common_stage1_args,
      stage1_extra_args,
      dots
    )

    if (isTRUE(verbose)) message("CDSeq5 Stage 1: fitting plug-in proportions and BETA.")
    withr::local_seed(seed)
    stage1_fit <- .cdseq5_call_stage1(
      stage1_function,
      stage1_args,
      strict = strict_stage1_args
    )
  }

  beta_found <- if (!is.null(reference_beta)) {
    list(value = reference_beta, path = "reference_beta argument")
  } else {
    if (is.null(stage1_fit)) {
      .cdseq5_stop("A Stage-1 fit is required when reference_beta is not supplied.")
    }
    .cdseq5_find_beta(stage1_fit)
  }

  cell_types <- .cdseq5_infer_cell_types(
    cell_types = cell_types,
    beta_raw = beta_found$value,
    genes = genes,
    reference_scRNA_annotation = reference_scRNA_annotation,
    fixed_proportions = fixed_proportions
  )
  beta <- .cdseq5_orient_beta(beta_found$value, cell_types, genes)

  if (is.null(fixed_proportions)) {
    if (is.null(stage1_fit)) {
      .cdseq5_stop("A Stage-1 fit is required when fixed_proportions is NULL.")
    }
    proportion_found <- .cdseq5_extract_stage1_proportions(
      stage1_fit = stage1_fit,
      cell_types = cell_types,
      sample_names = sample_names,
      component = proportion_component
    )
    proportion_source <- paste0("stage1:", proportion_found$path)
    fixed_input <- proportion_found$value
  } else {
    proportion_source <- "provided fixed_proportions"
    fixed_input <- fixed_proportions
  }

  fixed_prepared <- .cdseq5_prepare_fixed_proportions(
    fixed_input,
    cell_types,
    sample_names,
    proportion_floor
  )
  fixed_matrix <- fixed_prepared$adjusted

  stage1_assignments <- NULL
  if (identical(initialization, "stage1")) {
    if (is.null(stage1_fit)) {
      .cdseq5_stop("initialization = 'stage1' requires stage1_fit.")
    }
    stage1_assignments <- .cdseq5_extract_stage1_assignments(
      stage1_fit, sample_names
    )
  }

  sample_depths <- colSums(bulk_data)
  absolute_kappa <- .cdseq5_expand_sample_parameter(
    lock_concentration, sample_names, "lock_concentration"
  )

  if (is.null(absolute_kappa)) {
    lock_ratio_by_sample <- .cdseq5_expand_sample_parameter(
      lock_ratio, sample_names, "lock_ratio"
    )
    kappa_by_sample <- lock_ratio_by_sample * sample_depths
  } else {
    kappa_by_sample <- absolute_kappa
    lock_ratio_by_sample <- kappa_by_sample / sample_depths
  }

  if (anyNA(kappa_by_sample) || any(!is.finite(kappa_by_sample)) ||
      any(kappa_by_sample <= 0)) {
    .cdseq5_stop(
      "Stage-2 lock concentrations are nonfinite or nonpositive. Reduce the lock ratio."
    )
  }

  if (any(kappa_by_sample > 1e12)) {
    .cdseq5_warn(
      "Some Stage-2 Dirichlet concentrations exceed 1e12. This is usually ",
      "unnecessary and may create numerical problems."
    )
  }

  # Resolve once in the main process to fail early. Parallel workers resolve the
  # package-internal function again when gibbs_sampler is NULL.
  invisible(.cdseq5_resolve_gibbs_sampler(gibbs_sampler))

  if (isTRUE(verbose)) {
    message(
      "CDSeq5 Stage 2: strongly anchored Gibbs runs for ",
      length(sample_names), " sample(s); median lock ratio = ",
      format(stats::median(lock_ratio_by_sample), digits = 5), "."
    )
  }

  worker_args <- list(
    bulk_data = bulk_data,
    beta = beta,
    fixed_proportions = fixed_matrix,
    fixed_proportions_before_floor = fixed_prepared$normalized_raw,
    kappa_by_sample = kappa_by_sample,
    lock_ratio_by_sample = lock_ratio_by_sample,
    mcmc_iterations = mcmc_iterations,
    save_last_n_mc = save_last_n_mc,
    output_control = output_control,
    initialization = initialization,
    stage1_assignments = stage1_assignments,
    verbose_int = if (isTRUE(verbose)) 1L else 0L,
    rng = rng,
    seed = seed,
    return_mcmc_arrays = return_mcmc_arrays,
    return_assignments = return_assignments,
    gibbs_sampler = gibbs_sampler
  )

  run_index <- function(sample_index) {
    do.call(
      .cdseq5_run_one_sample,
      c(list(sample_index = sample_index), worker_args)
    )
  }

  if (cpu_number == 1L || length(sample_names) == 1L) {
    all_results <- lapply(seq_along(sample_names), run_index)
  } else {
    if (!requireNamespace("foreach", quietly = TRUE) ||
        !requireNamespace("doParallel", quietly = TRUE)) {
      .cdseq5_stop(
        "foreach and doParallel are required when cpu_number > 1."
      )
    }
    cluster <- parallel::makePSOCKcluster(min(cpu_number, length(sample_names)))
    on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
    doParallel::registerDoParallel(cluster)
    `%dopar%` <- foreach::`%dopar%`

    package_vector <- if (is.null(gibbs_sampler)) "CDSeq2" else character(0)
    # foreach() creates this iteration variable through non-standard
    # evaluation. Defining it locally keeps R CMD check/codetools from
    # treating it as an unresolved global binding.
    sample_index <- NULL
    all_results <- foreach::foreach(
      sample_index = seq_along(sample_names),
      .inorder = TRUE,
      .errorhandling = "stop",
      .packages = package_vector,
      .export = c(
        ".cdseq5_run_one_sample",
        ".cdseq5_resolve_gibbs_sampler",
        ".cdseq5_safe_column_normalize",
        ".cdseq5_stop"
      )
    ) %dopar% {
      do.call(
        .cdseq5_run_one_sample,
        c(list(sample_index = sample_index), worker_args)
      )
    }
    parallel::stopCluster(cluster)
    on.exit(NULL, add = FALSE)
  }

  names(all_results) <- vapply(all_results, `[[`, character(1), "sample_name")
  all_results <- all_results[sample_names]

  Gn <- length(genes)
  Tn <- length(cell_types)
  Mn <- length(sample_names)
  ctsGEP <- array(
    NA_real_, c(Gn, Tn, Mn),
    dimnames = list(gene = genes, cell_type = cell_types, sample = sample_names)
  )
  ctsGEP_normalized <- ctsGEP
  ctsGEP_normalized_draw_mean <- ctsGEP
  profile_total_mean <- matrix(
    NA_real_, Mn, Tn, dimnames = list(sample_names, cell_types)
  )
  profile_total_sd <- profile_total_mean
  allocation_proportions <- matrix(
    NA_real_, Mn, Tn, dimnames = list(sample_names, cell_types)
  )
  effective_proportions <- allocation_proportions

  for (sample_index in seq_along(sample_names)) {
    sample_name <- sample_names[[sample_index]]
    sample_result <- all_results[[sample_name]]
    ctsGEP[, , sample_name] <- sample_result$ctsGEP_counts
    ctsGEP_normalized[, , sample_name] <- sample_result$ctsGEP_normalized
    ctsGEP_normalized_draw_mean[, , sample_name] <-
      sample_result$ctsGEP_normalized_draw_mean
    profile_total_mean[sample_name, ] <-
      sample_result$profile_total_posterior_mean[cell_types]
    profile_total_sd[sample_name, ] <-
      sample_result$profile_total_posterior_sd[cell_types]
    allocation_proportions[sample_name, ] <-
      sample_result$allocation_proportion[cell_types, 1L]
    effective_proportions[sample_name, ] <-
      sample_result$effective_proportion[cell_types, 1L]
  }

  lock_diagnostics <- do.call(
    rbind,
    lapply(all_results, `[[`, "lock_diagnostic")
  )
  rownames(lock_diagnostics) <- NULL

  elapsed_seconds <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
  )

  result <- list(
    call = call,
    method = "CDSeq5",
    method_description = paste(
      "two-stage CDSeq with near-fixed, strongly anchored plug-in RNA proportions",
      "using the unchanged compiled Gibbs sampler"
    ),
    exact_fixed_proportions = FALSE,
    proportion_lock_method = "finite Dirichlet concentration centered on plug-in RNA proportions",
    base_sampler = "CDSeq2::gibbsSampler",
    proportion_source = proportion_source,
    stage1_proportion_component = proportion_component,
    beta_source = beta_found$path,
    reference_precision_by_cell_type = rowSums(beta),
    fixed_proportions = t(fixed_matrix),
    fixed_proportions_before_floor = t(fixed_prepared$normalized_raw),
    rna_proportions = t(fixed_matrix),
    rna_proportion_role = "plug-in Stage-2 input, not a newly estimated Stage-2 endpoint",
    allocation_proportions = allocation_proportions,
    effective_proportions = effective_proportions,
    ctsGEP = ctsGEP,
    ctsGEP_normalized = ctsGEP_normalized,
    ctsGEP_normalized_posthoc = ctsGEP_normalized,
    ctsGEP_normalized_draw_mean = ctsGEP_normalized_draw_mean,
    profile_total_posterior_mean = profile_total_mean,
    profile_total_posterior_sd = profile_total_sd,
    normalization_estimator =
      "normalize the posterior mean of retained allocated-count draws",
    lock_diagnostics = lock_diagnostics,
    lock_ratio_by_sample = lock_ratio_by_sample,
    lock_concentration_by_sample = kappa_by_sample,
    proportion_floor = proportion_floor,
    max_abs_proportion_floor_change = fixed_prepared$max_abs_floor_change,
    cdseq_result = all_results,
    parameters = list(
      beta = beta,
      alpha_by_sample = lapply(all_results, `[[`, "alpha"),
      cell_types = cell_types,
      genes = genes,
      mcmc_iterations = mcmc_iterations,
      save_last_n_mc = save_last_n_mc,
      rng = rng,
      seed = seed,
      cpu_number = cpu_number,
      lock_ratio_by_sample = lock_ratio_by_sample,
      lock_concentration_by_sample = kappa_by_sample
    ),
    elapsed_seconds = elapsed_seconds,
    stage1_fit = if (isTRUE(keep_stage1_fit)) stage1_fit else NULL
  )

  class(result) <- c("cdseq5_fit", "list")
  result
}

#' @export
print.cdseq5_fit <- function(x, ...) {
  cat("CDSeq5 fit\n")
  cat("  Strategy: near-fixed, strongly anchored plug-in RNA proportions\n")
  cat("  Exact fixed proportions: no\n")
  cat("  Proportion source: ", x$proportion_source, "\n", sep = "")
  cat(
    "  Median lock ratio: ",
    format(stats::median(x$lock_ratio_by_sample), digits = 5),
    "\n",
    sep = ""
  )
  cat(
    "  Median theoretical assignment-feedback weight: ",
    format(
      stats::median(
        x$lock_diagnostics$theoretical_assignment_feedback_weight
      ),
      digits = 5
    ),
    "\n",
    sep = ""
  )
  cat(
    "  ctsGEP dimensions: ", paste(dim(x$ctsGEP), collapse = " x "),
    " (gene x cell type x sample)\n",
    sep = ""
  )
  invisible(x)
}
