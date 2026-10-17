# ==============================================================================
# ROSMAP primary benchmark workflow
#
# Run from R with:
#   source("~/data-lab5/Project_CDSeq2/code/rosmap_primary_benchmark_workflow.R")
#
# This script intentionally performs the complete workflow when sourced.  To
# load only the functions (for example, from the sensitivity workflow), use:
#   options(rosmap.primary.run_on_source = FALSE)
#   source("rosmap_primary_benchmark_workflow.R")
#
# Primary design
#   1. Construct a balanced pathological-AD/non-AD subject pool.
#   2. Recluster each component with Harmony correction by subject.
#   3. Select donor-shared clusters without using condition labels.
#   4. Require complete truth for all seven major cell types.
#   5. Select a balanced, quality-constrained 20-subject reference.
#   6. Cap reference cells within donor and cell type.
#   7. Use every remaining eligible subject for expression recovery.
#   8. Select a balanced/covariate-matched held-out ctsDE subset.
#   9. Use RNA/read proportions, not cell proportions, as proportion truth.
#  10. Build one canonical reference split; repeat splits in the sensitivity workflow.
# ==============================================================================

rosmap_primary_config <- function() {
  data_dir <- path.expand("~/data-lab5/Project_CDSeq2/data/ROSMAP_AD")

  list(
    # Input files --------------------------------------------------------------
    input_files = c(
      Ast      = file.path(data_dir, "Astrocytes.rds"),
      Excit_1  = file.path(data_dir, "Excitatory_neurons_set1.rds"),
      Excit_2  = file.path(data_dir, "Excitatory_neurons_set2.rds"),
      Excit_3  = file.path(data_dir, "Excitatory_neurons_set3.rds"),
      OPC      = file.path(data_dir, "OPCs.rds"),
      Imm      = file.path(data_dir, "Immune_cells.rds"),
      Oligo    = file.path(data_dir, "Oligodendrocytes.rds"),
      Inhib    = file.path(data_dir, "Inhibitory_neurons.rds"),
      Vascular = file.path(data_dir, "Vasculature_cells.rds")
    ),
    metadata_file = file.path(data_dir, "individual_metadata_deidentified.tsv"),
    output_dir = file.path(
      data_dir,
      "ref_heldout_benchmark_manuscript_primary"
    ),

    # Metadata -----------------------------------------------------------------
    subject_col = "subject",
    metadata_subject_col = "subject",
    diagnosis_col = "Pathologic_diagnosis_of_AD",
    high_resolution_col = "cell_type_high_resolution",
    assay = "RNA",
    covariate_cols = c("msex", "age_death", "pmi"),

    # The benchmark is performed only at these seven major classes. Fine labels
    # are retained for auditing, not used as benchmark truth classes.
    major_celltypes = c(
      "Astro", "Excit", "Inhib", "Oligo", "OPC", "Immune", "Vascular"
    ),
    component_to_major = c(
      Ast = "Astro",
      Excit_1 = "Excit",
      Excit_2 = "Excit",
      Excit_3 = "Excit",
      OPC = "OPC",
      Imm = "Immune",
      Oligo = "Oligo",
      Inhib = "Inhib",
      Vascular = "Vascular"
    ),

    # Balanced starting pool ---------------------------------------------------
    pool_n_per_condition = 100L,
    pool_seed = 123L,

    # Recluster each of the nine input objects after selecting the 100 + 100
    # starting pool. Harmony uses subject only; condition is never supplied.
    recluster_enabled = TRUE,
    recluster_npcs = 50L,
    recluster_dims = 1:20,
    recluster_resolutions = c(0.05, 0.01),
    recluster_run_umap = FALSE,
    recluster_workers = 20L,
    reuse_cluster_checkpoint = TRUE,
    save_cluster_checkpoint = TRUE,
    reuse_legacy_clustered_objects = TRUE,
    legacy_clustered_file = path.expand(
      "~/data-lab5/Project_CDSeq2/data/seurat_list_clustered.rds"
    ),
    cluster_checkpoint_file = file.path(
      data_dir,
      "seurat_list_clustered_manuscript_pool.rds"
    ),

    # Donor-shared cluster selection. The primary workflow follows the latest
    # pipeline exactly: every one of the nine component objects is filtered.
    cluster_filter_enabled = TRUE,
    cluster_filter_units = c(
      "Ast", "Excit_1", "Excit_2", "Excit_3", "OPC",
      "Oligo", "Imm", "Inhib", "Vascular"
    ),
    cluster_cols = "RNA_snn_res.0.01",
    min_cells_per_subject_cluster = c(
      Ast = 50L,
      Excit_1 = 50L,
      Excit_2 = 50L,
      Excit_3 = 50L,
      Inhib = 50L,
      OPC = 50L,
      Oligo = 50L,
      Imm = 50L,
      Vascular = 20L
    ),
    top_n_clusters = 1L,
    # FALSE uses the cluster thresholds to rank/select robust clusters, but
    # does not require each retained donor to meet the threshold separately in
    # all nine component objects. Donor eligibility is imposed once below at
    # the seven-major-cell-type truth level. This avoids stacked filters.
    strict_cluster_subject_intersection = FALSE,

    # Complete-truth and donor-quality requirements ---------------------------
    # Complete truth means that every target has observed truth in every major
    # class. Higher per-class thresholds are imposed when drawing the reference,
    # not globally on targets; retaining low-abundance targets is necessary for
    # the planned abundance-stratified performance analysis.
    min_truth_cells_per_ct = c(
      Astro = 1L,
      Excit = 1L,
      Inhib = 1L,
      Oligo = 1L,
      OPC = 1L,
      Immune = 1L,
      Vascular = 1L
    ),
    min_total_cells_eligible = 500L,

    # Reference design ---------------------------------------------------------
    n_reference_subjects = 20L,
    reference_sampling = "balanced", # "balanced" or "control_only"
    min_cells_per_ind_ref = 20L,
    max_cells_per_ind_ref = 500L,
    min_reference_subject_fraction_per_celltype = 0.50,
    min_ref_cells_per_celltype = c(
      Astro = 1000L,
      Excit = 2500L,
      Inhib = 1500L,
      Oligo = 1000L,
      OPC = 500L,
      Immune = 500L,
      Vascular = 250L
    ),
    max_reference_draw_attempts = 10000L,
    # The primary workflow creates one canonical split. Repeated split analysis
    # belongs in rosmap_sensitivity_benchmark_workflow.R and is intentionally
    # not repeated here.
    reference_seeds = 123L,

    # Held-out ctsDE subset. NULL uses the largest possible balanced subset.
    ctsde_sampling = "covariate_matched", # or "balanced_random"
    ctsde_n_per_condition = NULL,

    # Gene panel. Markers are found using reference cells only, after the donor
    # split, so held-out condition information cannot leak into gene selection.
    min_ref_gene_total = 100,
    min_bulk_gene_total = 100,
    # Reference-only FindAllMarkers is computationally expensive and was not
    # part of the original benchmark-construction script. The prevalence panel
    # is used by default; enable this only for a dedicated marker sensitivity.
    marker_augmentation = FALSE,
    marker_min_pct = 0.10,
    marker_logfc_threshold = 0.25,
    marker_padj_threshold = 0.01,
    marker_top_n_per_celltype = 200L,

    # Counts are sufficient for CDSeq2. Avoiding a second logcounts assay and
    # using uncompressed RDS writes substantially reduces runtime and peak memory.
    add_logcounts_to_reference = FALSE,
    rds_compress = FALSE,

    # Reproducibility and behavior --------------------------------------------
    strict = TRUE,
    verbose = TRUE
  )
}


# ==============================================================================
# General utilities
# ==============================================================================

assert_packages <- function(cfg = NULL) {
  packages <- c(
    "Seurat", "SeuratObject", "SingleCellExperiment", "Matrix",
    "data.table", "harmony"
  )
  if (!is.null(cfg) && isTRUE(cfg$add_logcounts_to_reference)) {
    packages <- c(packages, "scater")
  }
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Install the following packages before running the workflow: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

expand_named_param <- function(x, required_names, label) {
  if (is.null(x) || length(x) == 0L) {
    stop(label, " must be supplied.", call. = FALSE)
  }
  if (length(x) == 1L && is.null(names(x))) {
    return(stats::setNames(rep(x, length(required_names)), required_names))
  }
  check_named_values(x, required_names, label)
  x[required_names]
}

msg <- function(..., cfg = NULL) {
  if (is.null(cfg) || isTRUE(cfg$verbose)) {
    message(...)
  }
}

check_named_values <- function(x, required_names, label) {
  missing <- setdiff(required_names, names(x))
  if (length(missing) > 0L) {
    stop(label, " is missing names: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (ncol(x) == 0L) {
    x <- data.frame(note = character(), stringsAsFactors = FALSE)
  }
  data.table::fwrite(data.table::as.data.table(x), path, sep = "\t", na = "NA")
  invisible(path)
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

collapse_ids <- function(x) {
  paste(sort(unique(as.character(x))), collapse = ";")
}

condition_counts <- function(subjects, subject_metadata, stage, seed = NA_integer_) {
  subjects <- unique(as.character(subjects))
  d <- subject_metadata[subject_metadata$subject %in% subjects, , drop = FALSE]
  tab <- as.data.frame(
    table(factor(d$condition, levels = c("AD", "Control"))),
    stringsAsFactors = FALSE
  )
  names(tab) <- c("condition", "n_subjects")
  tab$stage <- stage
  tab$seed <- seed
  tab[, c("seed", "stage", "condition", "n_subjects"), drop = FALSE]
}

quantile_bin <- function(x, labels = c("low", "middle", "high")) {
  x <- as.numeric(x)
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) {
    return(factor(out, levels = labels, ordered = TRUE))
  }
  q <- unique(as.numeric(stats::quantile(
    x[ok], probs = seq(0, 1, length.out = length(labels) + 1L),
    na.rm = TRUE, names = FALSE, type = 8
  )))
  if (length(q) < 2L) {
    out[ok] <- labels[ceiling(length(labels) / 2)]
  } else {
    used_labels <- labels[seq_len(length(q) - 1L)]
    out[ok] <- as.character(cut(
      x[ok], breaks = q, include.lowest = TRUE, labels = used_labels,
      ordered_result = TRUE
    ))
  }
  factor(out, levels = labels, ordered = TRUE)
}


# ==============================================================================
# Input loading and count extraction
# ==============================================================================

read_subject_metadata <- function(cfg) {
  if (!file.exists(cfg$metadata_file)) {
    stop("Metadata file not found: ", cfg$metadata_file, call. = FALSE)
  }

  meta <- data.table::fread(cfg$metadata_file, data.table = FALSE)
  required <- c(cfg$metadata_subject_col, cfg$diagnosis_col)
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0L) {
    stop("Metadata is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  meta$subject <- as.character(meta[[cfg$metadata_subject_col]])
  dx_ad <- tolower(trimws(as.character(meta[[cfg$diagnosis_col]])))
  if (!all(dx_ad %in% c("yes", "no"))) {
    bad <- sort(unique(dx_ad[!dx_ad %in% c("yes", "no")]))
    stop(
      "Unexpected values in ", cfg$diagnosis_col, ": ",
      paste(bad, collapse = ", "), ". Expected only yes/no.",
      call. = FALSE
    )
  }
  meta$condition <- ifelse(dx_ad == "yes", "AD", "Control")

  if (anyDuplicated(meta$subject)) {
    stop("Metadata contains duplicated subject IDs.", call. = FALSE)
  }
  rownames(meta) <- meta$subject
  meta
}

extract_count_matrix <- function(object, assay = "RNA") {
  if (!assay %in% names(object@assays)) {
    stop("Assay '", assay, "' is not present in a Seurat object.", call. = FALSE)
  }

  # SeuratObject 4.x uses slot = "counts". Calling the SeuratObject 5.x
  # layer interface on a 4.x object emits "arguments are not used: layer".
  # Keep the two APIs in separate branches so older package installations are
  # warning-free and do not pay for unsuccessful layer discovery.
  seurat_object_version <- utils::packageVersion("SeuratObject")
  if (seurat_object_version < "5.0.0") {
    counts <- suppressWarnings(SeuratObject::GetAssayData(
      object = object,
      assay = assay,
      slot = "counts"
    ))
    if (nrow(counts) == 0L || ncol(counts) == 0L) {
      stop("No non-empty raw-count matrix was found in assay '", assay, "'.", call. = FALSE)
    }
    return(methods::as(counts, "dgCMatrix"))
  }

  # In SeuratObject 5.x, a single counts layer can be retrieved directly. A
  # split assay may have several counts.* layers, in which case GetAssayData()
  # fails and the layer-aware assembly below is required.
  counts <- tryCatch(
    SeuratObject::GetAssayData(object = object, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (!is.null(counts) && nrow(counts) > 0L && ncol(counts) > 0L) {
    return(methods::as(counts, "dgCMatrix"))
  }

  assay_object <- object[[assay]]
  layers <- tryCatch(
    SeuratObject::Layers(assay_object, search = "^counts($|\\.)"),
    error = function(e) character()
  )

  if (length(layers) > 0L) {
    matrices <- lapply(layers, function(layer_name) {
      SeuratObject::LayerData(object, assay = assay, layer = layer_name)
    })
    matrices <- matrices[vapply(matrices, ncol, integer(1)) > 0L]
    if (length(matrices) == 0L) {
      stop("Count layers were found but contain no cells.", call. = FALSE)
    }
    if (length(matrices) == 1L) {
      return(methods::as(matrices[[1L]], "dgCMatrix"))
    }
    common_genes <- Reduce(intersect, lapply(matrices, rownames))
    if (length(common_genes) == 0L) {
      stop("Count layers have no shared genes.", call. = FALSE)
    }
    out <- do.call(cbind, lapply(matrices, function(m) m[common_genes, , drop = FALSE]))
    out <- out[, colnames(object), drop = FALSE]
    return(methods::as(out, "dgCMatrix"))
  }

  stop("No non-empty raw-count layer was found in assay '", assay, "'.", call. = FALSE)
}

load_rosmap_components <- function(cfg) {
  missing_files <- cfg$input_files[!file.exists(cfg$input_files)]
  if (length(missing_files) > 0L) {
    stop(
      "Input RDS files not found:\n", paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }

  check_named_values(cfg$component_to_major, names(cfg$input_files), "component_to_major")

  components <- lapply(names(cfg$input_files), function(component) {
    object <- readRDS(cfg$input_files[[component]])
    if (!inherits(object, "Seurat")) {
      stop(cfg$input_files[[component]], " is not a Seurat object.", call. = FALSE)
    }
    md <- object[[]]
    if (!cfg$subject_col %in% names(md)) {
      stop(component, " is missing subject column '", cfg$subject_col, "'.", call. = FALSE)
    }

    object$subject <- as.character(object[[]][[cfg$subject_col]])
    object$cellType <- unname(cfg$component_to_major[[component]])
    object$source_component <- component
    if (cfg$high_resolution_col %in% names(object[[]])) {
      object$cell_type_high_resolution_audit <- as.character(
        object[[]][[cfg$high_resolution_col]]
      )
    } else {
      object$cell_type_high_resolution_audit <- NA_character_
    }
    object
  })
  names(components) <- names(cfg$input_files)
  components
}

prefix_component_cells <- function(components) {
  out <- lapply(names(components), function(component) {
    object <- components[[component]]
    prefix <- paste0(component, "__")
    if (!all(startsWith(colnames(object), prefix))) {
      object <- SeuratObject::RenameCells(
        object,
        new.names = paste(component, colnames(object), sep = "__")
      )
    }
    object
  })
  names(out) <- names(components)
  out
}

load_legacy_clustered_components <- function(cfg) {
  if (!isTRUE(cfg$reuse_legacy_clustered_objects)) return(NULL)

  legacy_file <- path.expand(cfg$legacy_clustered_file)
  if (!file.exists(legacy_file)) return(NULL)

  components <- tryCatch(readRDS(legacy_file), error = function(e) NULL)
  expected <- names(cfg$input_files)
  valid <- is.list(components) && length(components) == length(expected)
  if (valid && is.null(names(components))) names(components) <- expected
  valid <- valid && all(expected %in% names(components))
  if (!valid) {
    warning(
      "Ignoring legacy clustered file because it is not the expected named nine-object list: ",
      legacy_file,
      call. = FALSE
    )
    return(NULL)
  }
  components <- components[expected]

  valid_objects <- vapply(components, inherits, logical(1), what = "Seurat")
  if (!all(valid_objects)) {
    warning(
      "Ignoring legacy clustered file because these entries are not Seurat objects: ",
      paste(expected[!valid_objects], collapse = ", "),
      call. = FALSE
    )
    return(NULL)
  }

  components <- lapply(expected, function(component) {
    object <- components[[component]]
    md <- object[[]]
    if (!cfg$subject_col %in% names(md)) {
      stop(
        "The legacy clustered object ", component, " is missing subject column '",
        cfg$subject_col, "'.",
        call. = FALSE
      )
    }
    if (!any(cfg$cluster_cols %in% names(md))) {
      stop(
        "The legacy clustered object ", component,
        " has none of the requested cluster columns: ",
        paste(cfg$cluster_cols, collapse = ", "),
        call. = FALSE
      )
    }

    object$subject <- as.character(md[[cfg$subject_col]])
    object$cellType <- unname(cfg$component_to_major[[component]])
    object$source_component <- component
    if (cfg$high_resolution_col %in% names(md)) {
      object$cell_type_high_resolution_audit <- as.character(
        md[[cfg$high_resolution_col]]
      )
    } else {
      object$cell_type_high_resolution_audit <- NA_character_
    }
    object
  })
  names(components) <- expected
  components
}

subset_component_subjects <- function(object, subjects) {
  object_subjects <- as.character(object$subject)
  if (all(object_subjects %in% subjects)) return(object)
  keep <- rownames(object[[]])[object_subjects %in% subjects]
  object[, keep]
}

select_balanced_pool <- function(observed_subjects, subject_metadata, cfg) {
  observed_subjects <- intersect(unique(as.character(observed_subjects)), subject_metadata$subject)
  candidates <- subject_metadata[subject_metadata$subject %in% observed_subjects, , drop = FALSE]
  by_condition <- split(candidates$subject, candidates$condition)
  if (!all(c("AD", "Control") %in% names(by_condition))) {
    stop("Both AD and Control subjects are required for the balanced starting pool.", call. = FALSE)
  }

  n_each <- min(
    as.integer(cfg$pool_n_per_condition),
    length(by_condition$AD),
    length(by_condition$Control)
  )
  if (n_each < cfg$pool_n_per_condition) {
    warning(
      "The requested ", cfg$pool_n_per_condition,
      " subjects per condition were not available; using ", n_each,
      " per condition.",
      call. = FALSE
    )
  }
  if (n_each < 1L) {
    stop("No balanced starting pool can be constructed.", call. = FALSE)
  }

  set.seed(cfg$pool_seed)
  sort(c(
    sample(by_condition$AD, n_each, replace = FALSE),
    sample(by_condition$Control, n_each, replace = FALSE)
  ))
}

recluster_one_component <- function(object, component, cfg) {
  set.seed(cfg$pool_seed)
  object <- Seurat::NormalizeData(object = object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object = object, verbose = FALSE)
  object <- Seurat::ScaleData(object = object, verbose = FALSE)

  npcs <- min(
    as.integer(cfg$recluster_npcs),
    ncol(object) - 1L,
    length(SeuratObject::VariableFeatures(object)) - 1L
  )
  if (!is.finite(npcs) || npcs < 2L) {
    stop("Too few cells or variable genes to recluster ", component, ".", call. = FALSE)
  }
  object <- Seurat::RunPCA(object = object, npcs = npcs, verbose = FALSE)

  dims <- cfg$recluster_dims[cfg$recluster_dims <= npcs]
  if (length(dims) < 2L) {
    stop("Too few PCs are available for Harmony in ", component, ".", call. = FALSE)
  }
  object <- harmony::RunHarmony(
    object = object,
    group.by.vars = "subject",
    reduction.use = "pca",
    verbose = FALSE
  )
  object <- Seurat::FindNeighbors(
    object = object,
    reduction = "harmony",
    dims = dims,
    verbose = FALSE
  )
  for (resolution in cfg$recluster_resolutions) {
    object <- Seurat::FindClusters(
      object = object,
      resolution = resolution,
      random.seed = cfg$pool_seed,
      verbose = FALSE
    )
  }
  if (isTRUE(cfg$recluster_run_umap)) {
    object <- Seurat::RunUMAP(
      object = object,
      reduction = "harmony",
      dims = dims,
      seed.use = cfg$pool_seed,
      verbose = FALSE
    )
  }
  object
}

recluster_rosmap_components <- function(
    components,
    subject_pool,
    cfg,
    observed_subjects = NULL) {
  if (!isTRUE(cfg$recluster_enabled)) {
    return(components)
  }

  checkpoint_file <- path.expand(cfg$cluster_checkpoint_file)
  if (isTRUE(cfg$reuse_cluster_checkpoint) && file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    valid <- is.list(checkpoint) &&
      !is.null(checkpoint$components) &&
      identical(sort(checkpoint$subject_pool), sort(subject_pool)) &&
      identical(names(checkpoint$components), names(components)) &&
      identical(checkpoint$resolutions, cfg$recluster_resolutions)
    if (valid) {
      msg("Using compatible reclustering checkpoint: ", checkpoint_file, cfg = cfg)
      return(checkpoint$components)
    }
    warning(
      "Ignoring an incompatible reclustering checkpoint and rebuilding it: ",
      checkpoint_file,
      call. = FALSE
    )
  }

  workers <- min(as.integer(cfg$recluster_workers), length(components))
  if (.Platform$OS.type == "windows") workers <- 1L
  component_names <- names(components)
  clustered <- parallel::mclapply(
    component_names,
    function(component) {
      recluster_one_component(components[[component]], component, cfg)
    },
    mc.cores = workers,
    mc.preschedule = TRUE
  )
  names(clustered) <- component_names

  failed <- vapply(clustered, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop(
      "Reclustering failed for: ", paste(component_names[failed], collapse = ", "),
      call. = FALSE
    )
  }

  if (isTRUE(cfg$save_cluster_checkpoint)) {
    dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      list(
        components = clustered,
        subject_pool = sort(subject_pool),
        observed_subjects = sort(unique(observed_subjects)),
        resolutions = cfg$recluster_resolutions,
        harmony_group = "subject",
        seed = cfg$pool_seed
      ),
      checkpoint_file,
      compress = cfg$rds_compress
    )
  }
  clustered
}

load_cluster_checkpoint_fast <- function(cfg, subject_metadata) {
  if (!isTRUE(cfg$recluster_enabled) || !isTRUE(cfg$reuse_cluster_checkpoint)) {
    return(NULL)
  }
  checkpoint_file <- path.expand(cfg$cluster_checkpoint_file)
  if (!file.exists(checkpoint_file)) return(NULL)

  checkpoint <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
  if (is.null(checkpoint) || is.null(checkpoint$components) || is.null(checkpoint$subject_pool)) {
    return(NULL)
  }
  valid <- identical(names(checkpoint$components), names(cfg$input_files)) &&
    identical(checkpoint$resolutions, cfg$recluster_resolutions) &&
    identical(as.integer(checkpoint$seed), as.integer(cfg$pool_seed))
  if (!valid) return(NULL)

  pool_metadata <- subject_metadata[
    match(checkpoint$subject_pool, subject_metadata$subject),
    ,
    drop = FALSE
  ]
  if (anyNA(pool_metadata$subject)) return(NULL)
  condition_n <- table(factor(pool_metadata$condition, levels = c("AD", "Control")))
  expected_n <- as.integer(cfg$pool_n_per_condition)
  if (!all(condition_n == expected_n)) return(NULL)

  checkpoint
}


# ==============================================================================
# Donor-shared cluster selection (condition-blind)
# ==============================================================================

choose_shared_cluster_column <- function(objects, candidates, celltype) {
  available_in_all <- candidates[vapply(candidates, function(candidate) {
    all(vapply(objects, function(object) candidate %in% names(object[[]]), logical(1)))
  }, logical(1))]

  if (length(available_in_all) == 0L) {
    stop(
      "No candidate cluster column is present in every ", celltype,
      " component. Checked: ", paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  available_in_all[[1L]]
}

cluster_selection_one_celltype <- function(
    objects,
    celltype,
    cluster_candidates,
    subject_pool,
    min_cells_per_subject_cluster,
    top_n_clusters) {
  cluster_col <- choose_shared_cluster_column(objects, cluster_candidates, celltype)

  md <- do.call(rbind, lapply(names(objects), function(component) {
    x <- objects[[component]][[]]
    data.frame(
      cell_id = rownames(x),
      component = component,
      subject = as.character(x$subject),
      cluster = as.character(x[[cluster_col]]),
      stringsAsFactors = FALSE
    )
  }))
  md <- md[md$subject %in% subject_pool & !is.na(md$cluster), , drop = FALSE]
  if (nrow(md) == 0L) {
    stop("No cells remain for cluster selection in ", celltype, ".", call. = FALSE)
  }

  counts <- as.data.frame(table(md$cluster, md$subject), stringsAsFactors = FALSE)
  names(counts) <- c("cluster", "subject", "n_cells")
  counts <- counts[counts$n_cells > 0L, , drop = FALSE]
  split_counts <- split(counts, counts$cluster)

  summary_rows <- lapply(names(split_counts), function(cluster_id) {
    x <- split_counts[[cluster_id]]
    passing <- x$n_cells >= min_cells_per_subject_cluster
    data.frame(
      cellType = celltype,
      cluster_col = cluster_col,
      cluster = cluster_id,
      min_cells_per_subject_cluster = min_cells_per_subject_cluster,
      n_subjects_observed = nrow(x),
      n_subjects_passing = sum(passing),
      total_cells = sum(x$n_cells),
      median_cells_passing = if (any(passing)) stats::median(x$n_cells[passing]) else 0,
      stringsAsFactors = FALSE
    )
  })
  cluster_summary <- do.call(rbind, summary_rows)
  cluster_summary <- cluster_summary[order(
    -cluster_summary$n_subjects_passing,
    -cluster_summary$total_cells,
    -cluster_summary$median_cells_passing,
    cluster_summary$cluster
  ), , drop = FALSE]

  passing_summary <- cluster_summary[
    cluster_summary$n_subjects_passing > 0L,
    ,
    drop = FALSE
  ]
  if (nrow(passing_summary) == 0L) {
    stop(
      "No subject passes the cluster threshold for any cluster in ",
      celltype, ".",
      call. = FALSE
    )
  }
  n_select <- min(as.integer(top_n_clusters), nrow(passing_summary))
  selected_clusters <- passing_summary$cluster[seq_len(n_select)]
  selected_counts <- counts[counts$cluster %in% selected_clusters, , drop = FALSE]

  eligible_by_cluster <- lapply(selected_clusters, function(cluster_id) {
    x <- selected_counts[selected_counts$cluster == cluster_id, , drop = FALSE]
    sort(x$subject[x$n_cells >= min_cells_per_subject_cluster])
  })
  names(eligible_by_cluster) <- selected_clusters
  eligible_subjects <- if (length(eligible_by_cluster) == 1L) {
    eligible_by_cluster[[1L]]
  } else {
    Reduce(intersect, eligible_by_cluster)
  }

  verification <- merge(
    expand.grid(
      subject = sort(unique(subject_pool)),
      cluster = selected_clusters,
      stringsAsFactors = FALSE
    ),
    selected_counts,
    by = c("subject", "cluster"),
    all.x = TRUE
  )
  verification$n_cells[is.na(verification$n_cells)] <- 0L
  verification$passes <- verification$n_cells >= min_cells_per_subject_cluster
  verification$cellType <- celltype
  verification$cluster_col <- cluster_col

  list(
    celltype = celltype,
    cluster_col = cluster_col,
    selected_clusters = selected_clusters,
    cluster_summary = cluster_summary,
    eligible_subjects = eligible_subjects,
    eligible_by_cluster = eligible_by_cluster,
    verification = verification
  )
}

apply_shared_cluster_filter <- function(components, subject_pool, cfg) {
  if (!isTRUE(cfg$cluster_filter_enabled)) {
    filtered <- lapply(components, subset_component_subjects, subjects = subject_pool)
    return(list(
      components = filtered,
      shared_subjects = sort(unique(subject_pool)),
      threshold_eligible_subjects = sort(unique(subject_pool)),
      selected_clusters = data.frame(
        cellType = character(), cluster_col = character(),
        rank = integer(), cluster = character()
      ),
      cluster_summaries = data.frame(
        cellType = character(), cluster_col = character(),
        cluster = character(), n_subjects_passing = integer(),
        total_cells = numeric()
      ),
      verification = data.frame(
        cellType = character(), subject = character(), cluster = character(),
        n_cells = integer(), passes = logical()
      ),
      settings = list(enabled = FALSE)
    ))
  }

  filter_units <- cfg$cluster_filter_units
  missing_units <- setdiff(filter_units, names(components))
  if (length(missing_units) > 0L) {
    stop(
      "cluster_filter_units contains unknown components: ",
      paste(missing_units, collapse = ", "),
      call. = FALSE
    )
  }
  min_cells_by_unit <- expand_named_param(
    cfg$min_cells_per_subject_cluster,
    filter_units,
    "min_cells_per_subject_cluster"
  )
  top_n_by_unit <- expand_named_param(
    cfg$top_n_clusters,
    filter_units,
    "top_n_clusters"
  )

  selections <- lapply(filter_units, function(filter_unit) {
    cluster_selection_one_celltype(
      objects = components[filter_unit],
      celltype = filter_unit,
      cluster_candidates = cfg$cluster_cols,
      subject_pool = subject_pool,
      min_cells_per_subject_cluster = min_cells_by_unit[[filter_unit]],
      top_n_clusters = top_n_by_unit[[filter_unit]]
    )
  })
  names(selections) <- filter_units

  threshold_eligible_subjects <- Reduce(
    intersect,
    lapply(selections, `[[`, "eligible_subjects")
  )
  threshold_eligible_subjects <- sort(unique(threshold_eligible_subjects))

  # The per-component thresholds identify clusters that are broadly shared
  # across donors. Requiring the same donor to pass those thresholds in every
  # one of the nine objects, and then requiring seven-class truth again, stacks
  # two eligibility filters. By default, retain all pool donors' cells from the
  # selected clusters and impose donor eligibility once at the major-class
  # truth stage. The strict legacy behavior remains available as a sensitivity.
  if (isTRUE(cfg$strict_cluster_subject_intersection)) {
    shared_subjects <- threshold_eligible_subjects
  } else {
    shared_subjects <- sort(unique(subject_pool))
  }
  if (length(shared_subjects) == 0L) {
    stop("No subjects satisfy the selected clusters across all filtered cell types.", call. = FALSE)
  }

  filtered <- lapply(names(components), function(component) {
    object <- components[[component]]
    md <- object[[]]
    keep <- as.character(md$subject) %in% shared_subjects
    if (component %in% filter_units) {
      selection <- selections[[component]]
      keep <- keep & as.character(md[[selection$cluster_col]]) %in% selection$selected_clusters
    }
    object[, rownames(md)[keep]]
  })
  names(filtered) <- names(components)

  selected_clusters <- do.call(rbind, lapply(selections, function(x) {
    data.frame(
      cellType = x$celltype,
      cluster_col = x$cluster_col,
      rank = seq_along(x$selected_clusters),
      cluster = x$selected_clusters,
      stringsAsFactors = FALSE
    )
  }))

  list(
    components = filtered,
    shared_subjects = shared_subjects,
    threshold_eligible_subjects = threshold_eligible_subjects,
    selected_clusters = selected_clusters,
    cluster_summaries = do.call(rbind, lapply(selections, `[[`, "cluster_summary")),
    verification = do.call(rbind, lapply(selections, `[[`, "verification")),
    eligible_subjects_by_celltype = lapply(selections, `[[`, "eligible_subjects"),
    eligible_subjects_by_celltype_cluster = lapply(selections, `[[`, "eligible_by_cluster"),
    settings = list(
      enabled = TRUE,
      strict_cluster_subject_intersection = cfg$strict_cluster_subject_intersection,
      cluster_cols = cfg$cluster_cols,
      min_cells_per_subject_cluster = min_cells_by_unit,
      top_n_clusters = top_n_by_unit
    )
  )
}


# ==============================================================================
# Major-cell-type matrices and truth eligibility
# ==============================================================================

combine_components_by_major <- function(components, cfg) {
  out <- lapply(cfg$major_celltypes, function(celltype) {
    component_names <- names(cfg$component_to_major)[cfg$component_to_major == celltype]
    objects <- components[component_names]
    matrices <- lapply(component_names, function(component) {
      counts <- extract_count_matrix(objects[[component]], assay = cfg$assay)
      prefix <- paste0(component, "__")
      if (!all(startsWith(colnames(counts), prefix))) {
        colnames(counts) <- paste(component, colnames(counts), sep = "__")
      }
      counts
    })
    names(matrices) <- component_names
    common_genes <- Reduce(intersect, lapply(matrices, rownames))
    if (length(common_genes) == 0L) {
      stop("No shared genes among components of ", celltype, ".", call. = FALSE)
    }
    counts <- do.call(cbind, lapply(matrices, function(x) x[common_genes, , drop = FALSE]))
    metadata <- lapply(component_names, function(component) {
      object <- objects[[component]]
      md <- object[[]][colnames(object), , drop = FALSE]
      prefix <- paste0(component, "__")
      cell_ids <- colnames(object)
      if (!all(startsWith(cell_ids, prefix))) {
        cell_ids <- paste(component, cell_ids, sep = "__")
      }
      rownames(md) <- cell_ids
      md
    })
    md <- do.call(rbind, metadata)
    md <- md[colnames(counts), , drop = FALSE]
    list(counts = methods::as(counts, "dgCMatrix"), metadata = md)
  })
  names(out) <- cfg$major_celltypes

  shared_genes <- Reduce(intersect, lapply(out, function(x) rownames(x$counts)))
  if (length(shared_genes) == 0L) {
    stop("The seven major cell types have no shared genes.", call. = FALSE)
  }
  lapply(out, function(x) {
    x$counts <- x$counts[shared_genes, , drop = FALSE]
    x
  })
}

aggregate_counts_by_subject <- function(counts, subjects, target_subjects) {
  subjects <- as.character(subjects)
  target_subjects <- as.character(target_subjects)
  keep <- subjects %in% target_subjects
  counts <- counts[, keep, drop = FALSE]
  subjects <- subjects[keep]

  membership <- Matrix::sparseMatrix(
    i = seq_along(subjects),
    j = match(subjects, target_subjects),
    x = 1,
    dims = c(length(subjects), length(target_subjects)),
    dimnames = list(colnames(counts), target_subjects)
  )
  aggregated <- counts %*% membership
  n_cells <- tabulate(match(subjects, target_subjects), nbins = length(target_subjects))
  names(n_cells) <- target_subjects
  list(counts = methods::as(aggregated, "dgCMatrix"), n_cells = n_cells)
}

subject_celltype_availability <- function(major_data, subjects, cfg) {
  rows <- lapply(cfg$major_celltypes, function(celltype) {
    x <- major_data[[celltype]]
    cell_subjects <- as.character(x$metadata$subject)
    keep <- cell_subjects %in% subjects
    kept_subjects <- cell_subjects[keep]
    n_cells <- tabulate(match(kept_subjects, subjects), nbins = length(subjects))

    # Only per-subject RNA totals are needed here. The former implementation
    # built a full genes-by-subject matrix for every cell type and immediately
    # reduced it with colSums(), which was a large avoidable multiplication.
    cell_rna <- as.numeric(Matrix::colSums(x$counts))[keep]
    n_counts <- numeric(length(subjects))
    if (length(cell_rna) > 0L) {
      grouped <- rowsum(
        matrix(cell_rna, ncol = 1L),
        group = kept_subjects,
        reorder = FALSE
      )
      n_counts[match(rownames(grouped), subjects)] <- as.numeric(grouped[, 1L])
    }
    data.frame(
      subject = subjects,
      cellType = celltype,
      n_cells = as.integer(n_cells),
      n_counts = n_counts,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

complete_truth_subjects <- function(availability, cfg) {
  check_named_values(cfg$min_truth_cells_per_ct, cfg$major_celltypes, "min_truth_cells_per_ct")
  availability$min_required <- cfg$min_truth_cells_per_ct[availability$cellType]
  availability$passes_truth <- availability$n_cells >= availability$min_required
  by_subject <- split(availability, availability$subject)
  subjects <- names(by_subject)[vapply(by_subject, function(x) {
    all(cfg$major_celltypes %in% x$cellType) && all(x$passes_truth)
  }, logical(1))]
  sort(subjects)
}

high_resolution_composition <- function(components) {
  md <- do.call(rbind, lapply(components, function(x) {
    z <- x[[]]
    data.frame(
      subject = as.character(z$subject),
      cellType = as.character(z$cellType),
      high_resolution_label = as.character(z$cell_type_high_resolution_audit),
      source_component = as.character(z$source_component),
      stringsAsFactors = FALSE
    )
  }))
  md$high_resolution_label[is.na(md$high_resolution_label) | md$high_resolution_label == ""] <- "<missing>"
  out <- as.data.frame(table(
    md$cellType, md$source_component, md$high_resolution_label
  ), stringsAsFactors = FALSE)
  names(out) <- c("cellType", "source_component", "high_resolution_label", "n_cells")
  out[out$n_cells > 0L, , drop = FALSE]
}


# ==============================================================================
# Reference and ctsDE subject selection
# ==============================================================================

reference_coverage <- function(subjects, availability, cfg) {
  d <- availability[availability$subject %in% subjects, , drop = FALSE]
  min_by_ct <- expand_named_param(
    cfg$min_cells_per_ind_ref,
    cfg$major_celltypes,
    "min_cells_per_ind_ref"
  )
  max_by_ct <- expand_named_param(
    cfg$max_cells_per_ind_ref,
    cfg$major_celltypes,
    "max_cells_per_ind_ref"
  )
  required_cells <- expand_named_param(
    cfg$min_ref_cells_per_celltype,
    cfg$major_celltypes,
    "min_ref_cells_per_celltype"
  )
  d$min_per_subject <- min_by_ct[d$cellType]
  d$max_per_subject <- max_by_ct[d$cellType]
  d$n_contrib <- ifelse(
    d$n_cells >= d$min_per_subject,
    pmin(d$n_cells, d$max_per_subject),
    0
  )
  out <- do.call(rbind, lapply(cfg$major_celltypes, function(celltype) {
    x <- d[d$cellType == celltype, , drop = FALSE]
    data.frame(
      cellType = celltype,
      n_contributing_subjects = sum(x$n_contrib > 0),
      total_reference_cells_expected = sum(x$n_contrib),
      stringsAsFactors = FALSE
    )
  }))
  out$required_subjects <- ceiling(
    cfg$min_reference_subject_fraction_per_celltype * length(subjects)
  )
  out$required_reference_cells <- as.numeric(required_cells[out$cellType])
  out$passes_subjects <- out$n_contributing_subjects >= out$required_subjects
  out$passes_cells <- out$total_reference_cells_expected >= out$required_reference_cells
  out$passes <- out$passes_subjects & out$passes_cells
  out
}

draw_reference_candidate <- function(eligible_subjects, subject_metadata, cfg) {
  d <- subject_metadata[match(eligible_subjects, subject_metadata$subject), , drop = FALSE]
  if (anyNA(d$subject)) {
    stop("Some eligible subjects are absent from subject metadata.", call. = FALSE)
  }

  n_ref <- as.integer(cfg$n_reference_subjects)
  if (identical(cfg$reference_sampling, "balanced")) {
    if (n_ref %% 2L != 0L) {
      stop("Balanced reference sampling requires an even n_reference_subjects.", call. = FALSE)
    }
    n_each <- n_ref %/% 2L
    ad <- d$subject[d$condition == "AD"]
    control <- d$subject[d$condition == "Control"]
    if (length(ad) < n_each || length(control) < n_each) {
      stop(
        "Not enough eligible subjects for a balanced ", n_ref,
        "-subject reference.",
        call. = FALSE
      )
    }
    return(sort(c(sample(ad, n_each), sample(control, n_each))))
  }

  if (identical(cfg$reference_sampling, "control_only")) {
    control <- d$subject[d$condition == "Control"]
    if (length(control) < n_ref) {
      stop("Not enough eligible Control subjects for the reference.", call. = FALSE)
    }
    return(sort(sample(control, n_ref)))
  }

  stop("Unknown reference_sampling: ", cfg$reference_sampling, call. = FALSE)
}

select_quality_constrained_reference <- function(
    eligible_subjects,
    subject_metadata,
    availability,
    cfg,
    seed) {
  set.seed(seed)
  best <- NULL
  best_score <- -Inf

  for (attempt in seq_len(cfg$max_reference_draw_attempts)) {
    candidate <- draw_reference_candidate(eligible_subjects, subject_metadata, cfg)
    coverage <- reference_coverage(candidate, availability, cfg)
    score <- min(coverage$n_contributing_subjects - coverage$required_subjects)
    if (score > best_score) {
      best <- list(subjects = candidate, coverage = coverage, attempt = attempt)
      best_score <- score
    }
    if (all(coverage$passes) && all(cfg$major_celltypes %in% coverage$cellType)) {
      return(list(subjects = candidate, coverage = coverage, attempt = attempt))
    }
  }

  stop(
    "No reference draw satisfied cell-type donor coverage after ",
    cfg$max_reference_draw_attempts, " attempts. Best minimum margin was ",
    best_score, ". Review the saved availability table or relax only a prespecified threshold.",
    call. = FALSE
  )
}

numeric_matching_matrix <- function(d, covariate_cols) {
  covariate_cols <- intersect(covariate_cols, names(d))
  if (length(covariate_cols) == 0L) {
    return(NULL)
  }

  pieces <- lapply(covariate_cols, function(v) {
    x <- d[[v]]
    if (is.numeric(x) || is.integer(x)) {
      x <- as.numeric(x)
      x[!is.finite(x)] <- stats::median(x[is.finite(x)], na.rm = TRUE)
      if (!any(is.finite(x))) x[] <- 0
      matrix(x, ncol = 1L, dimnames = list(NULL, v))
    } else {
      x <- addNA(factor(x), ifany = TRUE)
      mm <- stats::model.matrix(~ x - 1)
      colnames(mm) <- paste0(v, "__", make.names(levels(x)))
      mm
    }
  })
  mm <- do.call(cbind, pieces)
  sds <- apply(mm, 2L, stats::sd, na.rm = TRUE)
  mm <- mm[, is.finite(sds) & sds > 0, drop = FALSE]
  if (ncol(mm) == 0L) return(NULL)
  scale(mm)
}

select_ctsde_subjects <- function(target_subjects, subject_metadata, cfg, seed) {
  d <- subject_metadata[match(target_subjects, subject_metadata$subject), , drop = FALSE]
  ad <- which(d$condition == "AD")
  control <- which(d$condition == "Control")
  n_each <- min(length(ad), length(control))
  if (!is.null(cfg$ctsde_n_per_condition)) {
    n_each <- min(n_each, as.integer(cfg$ctsde_n_per_condition))
  }
  if (n_each < 2L) {
    stop("At least two held-out subjects per condition are required for ctsDE.", call. = FALSE)
  }

  set.seed(seed + 100000L)
  if (identical(cfg$ctsde_sampling, "balanced_random")) {
    selected <- c(sample(ad, n_each), sample(control, n_each))
    return(sort(d$subject[selected]))
  }
  if (!identical(cfg$ctsde_sampling, "covariate_matched")) {
    stop("Unknown ctsde_sampling: ", cfg$ctsde_sampling, call. = FALSE)
  }

  mm <- numeric_matching_matrix(d, cfg$covariate_cols)
  if (is.null(mm)) {
    warning(
      "No usable matching covariates were found; using balanced random ctsDE sampling.",
      call. = FALSE
    )
    selected <- c(sample(ad, n_each), sample(control, n_each))
    return(sort(d$subject[selected]))
  }

  # Greedy nearest-neighbor matching without replacement. Repeating the AD
  # order by seed makes split sensitivity explicit and reproducible.
  ad_order <- sample(ad, length(ad), replace = FALSE)
  available_control <- control
  pairs <- vector("list", min(length(ad_order), length(control)))
  k <- 0L
  for (i in ad_order) {
    if (length(available_control) == 0L || k >= n_each) break
    delta <- sweep(mm[available_control, , drop = FALSE], 2L, mm[i, ], "-")
    distance <- rowSums(delta^2)
    j <- available_control[which.min(distance)]
    k <- k + 1L
    pairs[[k]] <- c(i, j)
    available_control <- setdiff(available_control, j)
  }
  selected <- unlist(pairs[seq_len(k)], use.names = FALSE)
  sort(d$subject[selected])
}

covariate_balance_table <- function(subjects, subject_metadata, covariate_cols, set_name) {
  d <- subject_metadata[subject_metadata$subject %in% subjects, , drop = FALSE]
  covariate_cols <- intersect(covariate_cols, names(d))
  rows <- list()
  index <- 0L

  for (v in covariate_cols) {
    x <- d[[v]]
    if (is.numeric(x) || is.integer(x)) {
      x <- as.numeric(x)
      ad <- x[d$condition == "AD"]
      control <- x[d$condition == "Control"]
      pooled_sd <- sqrt((stats::var(ad, na.rm = TRUE) + stats::var(control, na.rm = TRUE)) / 2)
      smd <- if (is.finite(pooled_sd) && pooled_sd > 0) {
        (mean(ad, na.rm = TRUE) - mean(control, na.rm = TRUE)) / pooled_sd
      } else {
        NA_real_
      }
      index <- index + 1L
      rows[[index]] <- data.frame(
        set = set_name,
        covariate = v,
        level = NA_character_,
        ad_value = mean(ad, na.rm = TRUE),
        control_value = mean(control, na.rm = TRUE),
        standardized_difference = smd,
        stringsAsFactors = FALSE
      )
    } else {
      levels_v <- sort(unique(as.character(x[!is.na(x)])))
      for (level in levels_v) {
        p_ad <- mean(as.character(x[d$condition == "AD"]) == level, na.rm = TRUE)
        p_control <- mean(as.character(x[d$condition == "Control"]) == level, na.rm = TRUE)
        pooled <- (p_ad + p_control) / 2
        denom <- sqrt(pooled * (1 - pooled))
        index <- index + 1L
        rows[[index]] <- data.frame(
          set = set_name,
          covariate = v,
          level = level,
          ad_value = p_ad,
          control_value = p_control,
          standardized_difference = if (denom > 0) (p_ad - p_control) / denom else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      set = character(), covariate = character(), level = character(),
      ad_value = numeric(), control_value = numeric(),
      standardized_difference = numeric()
    ))
  }
  do.call(rbind, rows)
}


# ==============================================================================
# Reference construction, marker selection, and held-out truth
# ==============================================================================

sample_reference_cells <- function(major_data, reference_subjects, cfg, seed) {
  matrices <- list()
  metadata <- list()
  contribution <- list()
  k <- 0L
  min_by_ct <- expand_named_param(
    cfg$min_cells_per_ind_ref,
    cfg$major_celltypes,
    "min_cells_per_ind_ref"
  )
  max_by_ct <- expand_named_param(
    cfg$max_cells_per_ind_ref,
    cfg$major_celltypes,
    "max_cells_per_ind_ref"
  )

  set.seed(seed + 200000L)
  for (celltype in cfg$major_celltypes) {
    x <- major_data[[celltype]]
    md <- x$metadata
    for (subject in reference_subjects) {
      cell_ids <- rownames(md)[as.character(md$subject) == subject]
      n_available <- length(cell_ids)
      if (n_available < min_by_ct[[celltype]]) {
        n_retained <- 0L
        selected <- character()
      } else {
        n_retained <- min(n_available, as.integer(max_by_ct[[celltype]]))
        selected <- if (n_available > n_retained) {
          sample(cell_ids, n_retained, replace = FALSE)
        } else {
          cell_ids
        }
      }

      k <- k + 1L
      contribution[[k]] <- data.frame(
        subject = subject,
        cellType = celltype,
        n_available = n_available,
        n_retained = n_retained,
        stringsAsFactors = FALSE
      )
      if (n_retained > 0L) {
        matrices[[length(matrices) + 1L]] <- x$counts[, selected, drop = FALSE]
        metadata[[length(metadata) + 1L]] <- md[selected, , drop = FALSE]
      }
    }
  }

  counts <- do.call(cbind, matrices)
  coldata <- do.call(rbind, metadata)
  coldata <- coldata[colnames(counts), , drop = FALSE]
  list(
    counts = methods::as(counts, "dgCMatrix"),
    coldata = coldata,
    contribution = do.call(rbind, contribution)
  )
}

find_reference_markers <- function(reference_counts, reference_coldata, cfg) {
  if (!isTRUE(cfg$marker_augmentation)) {
    return(data.frame())
  }

  object <- Seurat::CreateSeuratObject(
    counts = reference_counts,
    meta.data = reference_coldata
  )
  SeuratObject::Idents(object) <- factor(
    object$cellType,
    levels = cfg$major_celltypes
  )
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  markers <- Seurat::FindAllMarkers(
    object = object,
    only.pos = TRUE,
    min.pct = cfg$marker_min_pct,
    logfc.threshold = cfg$marker_logfc_threshold,
    return.thresh = cfg$marker_padj_threshold,
    verbose = FALSE
  )
  markers <- as.data.frame(markers)
  if (nrow(markers) == 0L) return(markers)

  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(markers))
  if (length(fc_col) == 0L || !"gene" %in% names(markers)) {
    stop("FindAllMarkers returned an unexpected column structure.", call. = FALSE)
  }
  fc_col <- fc_col[[1L]]
  markers <- markers[order(markers$cluster, -markers[[fc_col]]), , drop = FALSE]
  markers <- do.call(rbind, lapply(split(markers, markers$cluster), function(x) {
    utils::head(x, cfg$marker_top_n_per_celltype)
  }))
  rownames(markers) <- NULL
  markers
}

build_heldout_truth <- function(major_data, target_subjects, reference_sample, cfg) {
  target_by_ct <- lapply(cfg$major_celltypes, function(celltype) {
    x <- major_data[[celltype]]
    aggregate_counts_by_subject(
      counts = x$counts,
      subjects = x$metadata$subject,
      target_subjects = target_subjects
    )
  })
  names(target_by_ct) <- cfg$major_celltypes

  bulk_counts_all <- Reduce(`+`, lapply(target_by_ct, `[[`, "counts"))
  ref_total <- Matrix::rowSums(reference_sample$counts)
  bulk_total <- Matrix::rowSums(bulk_counts_all)
  prevalence_genes <- names(ref_total)[
    ref_total >= cfg$min_ref_gene_total & bulk_total >= cfg$min_bulk_gene_total
  ]

  markers <- find_reference_markers(
    reference_counts = reference_sample$counts,
    reference_coldata = reference_sample$coldata,
    cfg = cfg
  )
  marker_genes <- if (nrow(markers) > 0L) unique(as.character(markers$gene)) else character()
  marker_genes <- intersect(marker_genes, names(ref_total)[ref_total > 0 & bulk_total > 0])
  keep_genes <- rownames(reference_sample$counts)[
    rownames(reference_sample$counts) %in% union(prevalence_genes, marker_genes)
  ]
  if (length(keep_genes) == 0L) {
    stop("No genes passed the reference/bulk filters or marker augmentation.", call. = FALSE)
  }

  reference_counts <- reference_sample$counts[keep_genes, , drop = FALSE]
  target_by_ct <- lapply(target_by_ct, function(x) {
    x$counts <- x$counts[keep_genes, , drop = FALSE]
    x
  })
  bulk_counts <- Reduce(`+`, lapply(target_by_ct, `[[`, "counts"))

  truth_cell_counts <- do.call(cbind, lapply(target_by_ct, function(x) x$n_cells))
  rownames(truth_cell_counts) <- target_subjects
  colnames(truth_cell_counts) <- cfg$major_celltypes

  truth_rna_counts <- do.call(cbind, lapply(target_by_ct, function(x) {
    as.numeric(Matrix::colSums(x$counts))
  }))
  rownames(truth_rna_counts) <- target_subjects
  colnames(truth_rna_counts) <- cfg$major_celltypes
  truth_rna_proportions <- truth_rna_counts / rowSums(truth_rna_counts)

  truth_expression_counts <- setNames(vector("list", length(target_subjects)), target_subjects)
  truth_expression_metadata <- vector("list", length(target_subjects) * length(cfg$major_celltypes))
  k <- 0L
  for (j in seq_along(target_subjects)) {
    subject <- target_subjects[[j]]
    subject_matrix <- do.call(cbind, lapply(target_by_ct, function(x) {
      x$counts[, j, drop = FALSE]
    }))
    colnames(subject_matrix) <- cfg$major_celltypes
    truth_expression_counts[[subject]] <- methods::as(subject_matrix, "dgCMatrix")

    for (celltype in cfg$major_celltypes) {
      k <- k + 1L
      truth_expression_metadata[[k]] <- data.frame(
        col_id = paste(subject, celltype, sep = "||"),
        subject = subject,
        cellType = celltype,
        n_cells = truth_cell_counts[subject, celltype],
        n_counts = truth_rna_counts[subject, celltype],
        stringsAsFactors = FALSE
      )
    }
  }
  truth_expression_metadata <- do.call(rbind, truth_expression_metadata)

  gene_filter_table <- data.frame(
    gene = rownames(reference_sample$counts),
    reference_total = as.numeric(ref_total),
    heldout_bulk_total = as.numeric(bulk_total),
    passes_prevalence = rownames(reference_sample$counts) %in% prevalence_genes,
    selected_reference_marker = rownames(reference_sample$counts) %in% marker_genes,
    retained = rownames(reference_sample$counts) %in% keep_genes,
    stringsAsFactors = FALSE
  )

  list(
    reference_counts = reference_counts,
    pseudobulk_counts = methods::as(bulk_counts, "dgCMatrix"),
    truth_expression_counts = truth_expression_counts,
    truth_expression_metadata = truth_expression_metadata,
    truth_cell_counts = truth_cell_counts,
    truth_rna_counts = truth_rna_counts,
    truth_rna_proportions = truth_rna_proportions,
    markers = markers,
    gene_filter_table = gene_filter_table
  )
}

make_subject_output_metadata <- function(
    expression_target_subjects,
    ctsde_subjects,
    reference_subjects,
    subject_metadata,
    availability,
    reference_contribution,
    truth,
    cfg) {
  all_subjects <- c(reference_subjects, expression_target_subjects)
  d <- subject_metadata[
    match(all_subjects, subject_metadata$subject),
    ,
    drop = FALSE
  ]
  d$split <- ifelse(
    d$subject %in% reference_subjects,
    "reference_subjects",
    ifelse(d$subject %in% ctsde_subjects, "heldout_ctsDE", "heldout_expression_only")
  )
  d$is_ctsde_target <- d$subject %in% ctsde_subjects
  d$pseudobulk_n_counts <- NA_real_
  pb_depth <- Matrix::colSums(truth$pseudobulk_counts)
  d$pseudobulk_n_counts[match(names(pb_depth), d$subject)] <- as.numeric(pb_depth)

  available_rna_total <- stats::aggregate(n_counts ~ subject, availability, sum)

  for (celltype in cfg$major_celltypes) {
    suffix <- safe_filename(celltype)
    av <- availability[availability$cellType == celltype, , drop = FALSE]
    idx <- match(d$subject, av$subject)
    d[[paste0("n_available_", suffix)]] <- av$n_cells[idx]
    d[[paste0("rna_counts_available_", suffix)]] <- av$n_counts[idx]
    denom <- available_rna_total$n_counts[match(d$subject, available_rna_total$subject)]
    d[[paste0("rna_prop_available_", suffix)]] <- av$n_counts[idx] / denom

    rc <- reference_contribution[reference_contribution$cellType == celltype, , drop = FALSE]
    reference_column <- paste0("n_reference_", suffix)
    d[[reference_column]] <- 0L
    ref_idx <- match(d$subject, rc$subject)
    has_ref <- !is.na(ref_idx)
    d[[reference_column]][has_ref] <- rc$n_retained[ref_idx[has_ref]]
  }
  d
}

make_evaluation_strata <- function(truth, target_metadata, cfg) {
  x <- truth$truth_expression_metadata
  x$condition <- target_metadata$condition[match(x$subject, target_metadata$subject)]
  x$is_ctsde_target <- target_metadata$is_ctsde_target[match(x$subject, target_metadata$subject)]
  x$target_rna_proportion <- truth$truth_rna_proportions[
    cbind(match(x$subject, rownames(truth$truth_rna_proportions)), match(x$cellType, colnames(truth$truth_rna_proportions)))
  ]
  depth <- Matrix::colSums(truth$pseudobulk_counts)
  x$pseudobulk_depth <- as.numeric(depth[match(x$subject, names(depth))])
  x$rna_abundance_bin <- quantile_bin(x$target_rna_proportion)
  subject_depth <- data.frame(
    subject = names(depth),
    pseudobulk_depth = as.numeric(depth),
    stringsAsFactors = FALSE
  )
  subject_depth$pseudobulk_depth_bin <- quantile_bin(subject_depth$pseudobulk_depth)
  x$pseudobulk_depth_bin <- subject_depth$pseudobulk_depth_bin[
    match(x$subject, subject_depth$subject)
  ]
  x
}


# ==============================================================================
# Preparation, one-seed run, and complete primary workflow
# ==============================================================================

prepare_rosmap_benchmark <- function(cfg) {
  assert_packages(cfg)
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

  msg("1/9 Loading metadata and checking existing clustered objects.", cfg = cfg)
  subject_metadata <- read_subject_metadata(cfg)
  legacy_components <- load_legacy_clustered_components(cfg)
  used_legacy_components <- !is.null(legacy_components)

  if (used_legacy_components) {
    # This is the exact checkpoint created by the original pipeline. Reusing it
    # avoids the dominant NormalizeData/ScaleData/PCA/Harmony work entirely.
    observed_subjects <- unique(unlist(lapply(
      legacy_components,
      function(x) as.character(x$subject)
    )))
    observed_subjects <- intersect(observed_subjects, subject_metadata$subject)
    msg("2/9 Reusing the balanced subject pool in the original clustered file.", cfg = cfg)
    subject_pool <- select_balanced_pool(observed_subjects, subject_metadata, cfg)
    components <- lapply(
      legacy_components,
      subset_component_subjects,
      subjects = subject_pool
    )
    msg("3/9 Reusing completed Harmony clustering from: ",
        cfg$legacy_clustered_file, cfg = cfg)
    rm(legacy_components)
    invisible(gc())
  } else {
    checkpoint <- load_cluster_checkpoint_fast(cfg, subject_metadata)
  }

  if (!used_legacy_components && !is.null(checkpoint)) {
    msg("2/9 Reusing the balanced subject pool stored in the checkpoint.", cfg = cfg)
    subject_pool <- sort(as.character(checkpoint$subject_pool))
    components <- checkpoint$components
    observed_subjects <- checkpoint$observed_subjects
    if (is.null(observed_subjects) || length(observed_subjects) == 0L) {
      observed_subjects <- subject_pool
    }
    msg("3/9 Reusing completed Harmony reclustering from the checkpoint.", cfg = cfg)
    rm(checkpoint)
  } else if (!used_legacy_components) {
    components <- load_rosmap_components(cfg)
    observed_subjects <- unique(unlist(lapply(
      components,
      function(x) as.character(x$subject)
    )))
    observed_subjects <- intersect(observed_subjects, subject_metadata$subject)

    msg("2/9 Constructing the balanced pathological-AD/non-AD starting pool.", cfg = cfg)
    subject_pool <- select_balanced_pool(observed_subjects, subject_metadata, cfg)
    components <- lapply(components, subset_component_subjects, subjects = subject_pool)
    invisible(gc())

    msg("3/9 Re-clustering each component with Harmony correction by subject.", cfg = cfg)
    components <- recluster_rosmap_components(
      components,
      subject_pool,
      cfg,
      observed_subjects = observed_subjects
    )
  }

  msg("4/9 Selecting donor-shared clusters without using condition.", cfg = cfg)
  filtered <- apply_shared_cluster_filter(components, subject_pool, cfg)

  msg("5/9 Combining components into seven major cell types.", cfg = cfg)
  major_data <- combine_components_by_major(filtered$components, cfg)
  availability <- subject_celltype_availability(
    major_data,
    subjects = filtered$shared_subjects,
    cfg = cfg
  )
  complete_subjects <- complete_truth_subjects(availability, cfg)
  filtered_subjects <- filtered$shared_subjects
  cluster_threshold_subjects <- filtered$threshold_eligible_subjects
  total_cells <- stats::aggregate(n_cells ~ subject, availability, sum)
  names(total_cells)[2L] <- "n_available_total"
  eligible_subjects <- intersect(
    complete_subjects,
    total_cells$subject[total_cells$n_available_total >= cfg$min_total_cells_eligible]
  )
  eligible_subjects <- sort(eligible_subjects)
  if (length(eligible_subjects) <= cfg$n_reference_subjects) {
    stop(
      "Too few eligible subjects remain for a ", cfg$n_reference_subjects,
      "-subject reference: ", length(filtered_subjects),
      " donors retained selected-cluster cells; ", length(complete_subjects),
      " have complete seven-class truth; ", length(eligible_subjects),
      " also meet min_total_cells_eligible = ", cfg$min_total_cells_eligible,
      ". See the cluster verification and availability settings; do not apply ",
      "strict_cluster_subject_intersection in the primary workflow.",
      call. = FALSE
    )
  }

  subject_metadata <- merge(
    subject_metadata,
    total_cells,
    by = "subject",
    all.x = TRUE,
    sort = FALSE
  )
  rownames(subject_metadata) <- subject_metadata$subject

  common_attrition <- do.call(rbind, list(
    condition_counts(observed_subjects, subject_metadata, "observed_with_metadata"),
    condition_counts(subject_pool, subject_metadata, "balanced_starting_pool"),
    condition_counts(subject_pool, subject_metadata, "post_harmony_reclustering"),
    condition_counts(filtered_subjects, subject_metadata, "post_cluster_filter"),
    condition_counts(
      cluster_threshold_subjects,
      subject_metadata,
      "diagnostic_strict_9_component_intersection"
    ),
    condition_counts(complete_subjects, subject_metadata, "complete_truth_all_7"),
    condition_counts(eligible_subjects, subject_metadata, "reference_quality_eligible")
  ))

  audit_dir <- file.path(cfg$output_dir, "common_audit")
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(common_attrition, file.path(audit_dir, "subject_attrition_common.tsv"))
  write_tsv(availability, file.path(audit_dir, "cells_and_rna_by_subject_celltype.tsv"))
  write_tsv(filtered$selected_clusters, file.path(audit_dir, "selected_clusters.tsv"))
  write_tsv(filtered$cluster_summaries, file.path(audit_dir, "cluster_summaries.tsv"))
  write_tsv(filtered$verification, file.path(audit_dir, "cluster_verification.tsv"))
  write_tsv(
    high_resolution_composition(filtered$components),
    file.path(audit_dir, "high_resolution_label_composition.tsv")
  )
  # Do not serialize the large Seurat objects a second time. The reclustering
  # checkpoint already contains them; this file stores diagnostics only.
  cluster_filter_diagnostics <- filtered
  cluster_filter_diagnostics$components <- NULL
  saveRDS(
    cluster_filter_diagnostics,
    file.path(audit_dir, "cluster_filter_diagnostics.rds"),
    compress = cfg$rds_compress
  )
  rm(filtered, components)
  invisible(gc())

  list(
    cfg = cfg,
    subject_metadata = subject_metadata,
    major_data = major_data,
    availability = availability,
    observed_subjects = observed_subjects,
    subject_pool = subject_pool,
    filtered_subjects = filtered_subjects,
    complete_subjects = complete_subjects,
    eligible_subjects = eligible_subjects,
    common_attrition = common_attrition,
    cluster_filter = cluster_filter_diagnostics
  )
}

run_rosmap_seed <- function(prepared, seed, run_label = "primary") {
  cfg <- prepared$cfg
  seed <- as.integer(seed)
  seed_dir <- file.path(cfg$output_dir, paste0("reference_seed_", seed))
  audit_dir <- file.path(seed_dir, "audit")
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

  msg("6/9 Selecting the quality-constrained reference for seed ", seed, ".", cfg = cfg)
  reference <- select_quality_constrained_reference(
    eligible_subjects = prepared$eligible_subjects,
    subject_metadata = prepared$subject_metadata,
    availability = prepared$availability,
    cfg = cfg,
    seed = seed
  )
  expression_targets <- sort(setdiff(prepared$eligible_subjects, reference$subjects))
  ctsde_subjects <- select_ctsde_subjects(
    target_subjects = expression_targets,
    subject_metadata = prepared$subject_metadata,
    cfg = cfg,
    seed = seed
  )

  msg("7/9 Capping reference cells within donor and cell type.", cfg = cfg)
  reference_sample <- sample_reference_cells(
    major_data = prepared$major_data,
    reference_subjects = reference$subjects,
    cfg = cfg,
    seed = seed
  )

  msg("8/9 Building pseudobulk, expression truth, and RNA-proportion truth.", cfg = cfg)
  truth <- build_heldout_truth(
    major_data = prepared$major_data,
    target_subjects = expression_targets,
    reference_sample = reference_sample,
    cfg = cfg
  )
  truth$truth_expression_metadata$condition <- prepared$subject_metadata$condition[
    match(
      truth$truth_expression_metadata$subject,
      prepared$subject_metadata$subject
    )
  ]

  reference_sample$coldata$condition <- factor(
    prepared$subject_metadata$condition[
      match(reference_sample$coldata$subject, prepared$subject_metadata$subject)
    ],
    levels = c("Control", "AD")
  )
  reference_sample$coldata$split <- "reference_subjects"
  reference_sample$coldata$cellType <- factor(
    reference_sample$coldata$cellType,
    levels = cfg$major_celltypes
  )
  reference_sample$coldata$cellType_major <- reference_sample$coldata$cellType
  reference_sample$coldata$cellType_fine <- as.character(
    reference_sample$coldata$cell_type_high_resolution_audit
  )
  reference_sample$coldata$cell_type_high_resolution <- reference_sample$coldata$cellType_fine

  reference_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = truth$reference_counts),
    colData = S4Vectors::DataFrame(reference_sample$coldata)
  )
  if (isTRUE(cfg$add_logcounts_to_reference)) {
    reference_sce <- scater::logNormCounts(reference_sce)
  }
  SingleCellExperiment::metadata(reference_sce)$reference_subjects <- reference$subjects
  SingleCellExperiment::metadata(reference_sce)$reference_seed <- seed
  SingleCellExperiment::metadata(reference_sce)$major_celltypes <- cfg$major_celltypes

  reference_contribution <- merge(
    reference_sample$contribution,
    prepared$subject_metadata[, c("subject", "condition"), drop = FALSE],
    by = "subject",
    all.x = TRUE
  )
  benchmark_subject_metadata <- make_subject_output_metadata(
    expression_target_subjects = expression_targets,
    ctsde_subjects = ctsde_subjects,
    reference_subjects = reference$subjects,
    subject_metadata = prepared$subject_metadata,
    availability = prepared$availability,
    reference_contribution = reference_contribution,
    truth = truth,
    cfg = cfg
  )
  target_metadata <- benchmark_subject_metadata[
    benchmark_subject_metadata$subject %in% expression_targets,
    ,
    drop = FALSE
  ]
  evaluation_strata <- make_evaluation_strata(truth, target_metadata, cfg)

  benchmark <- list(
    pseudobulk_counts = truth$pseudobulk_counts,
    truth_expression_counts = truth$truth_expression_counts,
    truth_expression_metadata = truth$truth_expression_metadata,
    truth_cell_counts = truth$truth_cell_counts,
    truth_rna_counts = truth$truth_rna_counts,
    truth_rna_proportions = truth$truth_rna_proportions,
    target_metadata = target_metadata,
    subject_metadata = benchmark_subject_metadata,
    expression_target_subjects = expression_targets,
    ctsde_target_subjects = ctsde_subjects,
    reference_subjects = reference$subjects,
    reference_seed = seed,
    major_celltypes = cfg$major_celltypes,
    gene_panel = rownames(truth$pseudobulk_counts),
    settings = cfg
  )

  saveRDS(
    reference_sce,
    file.path(seed_dir, "rosmap_reference_sce.rds"),
    compress = cfg$rds_compress
  )
  saveRDS(
    benchmark,
    file.path(seed_dir, "rosmap_heldout_benchmark_minimal.rds"),
    compress = cfg$rds_compress
  )
  write_tsv(
    benchmark_subject_metadata,
    file.path(seed_dir, "rosmap_benchmark_subject_metadata.tsv")
  )

  attrition <- rbind(
    transform(prepared$common_attrition, seed = seed),
    condition_counts(reference$subjects, prepared$subject_metadata, "reference", seed),
    condition_counts(expression_targets, prepared$subject_metadata, "expression_targets_all", seed),
    condition_counts(ctsde_subjects, prepared$subject_metadata, "ctsDE_targets_balanced_matched", seed)
  )
  subject_ids <- rbind(
    data.frame(seed = seed, split = "reference", subject = reference$subjects),
    data.frame(seed = seed, split = "expression_target", subject = expression_targets),
    data.frame(seed = seed, split = "ctsDE_target", subject = ctsde_subjects)
  )
  reference_contribution$n_reference_counts <- 0
  for (i in seq_len(nrow(reference_contribution))) {
    sid <- reference_contribution$subject[[i]]
    ct <- reference_contribution$cellType[[i]]
    cols <- colnames(truth$reference_counts)[
      reference_sample$coldata[colnames(truth$reference_counts), "subject"] == sid &
        reference_sample$coldata[colnames(truth$reference_counts), "cellType"] == ct
    ]
    if (length(cols) > 0L) {
      reference_contribution$n_reference_counts[[i]] <- sum(truth$reference_counts[, cols, drop = FALSE])
    }
  }

  balance <- rbind(
    covariate_balance_table(
      expression_targets,
      prepared$subject_metadata,
      cfg$covariate_cols,
      "all_expression_targets"
    ),
    covariate_balance_table(
      ctsde_subjects,
      prepared$subject_metadata,
      cfg$covariate_cols,
      "ctsDE_targets"
    )
  )

  write_tsv(attrition, file.path(audit_dir, "subject_attrition_by_condition.tsv"))
  write_tsv(subject_ids, file.path(audit_dir, "exact_subject_ids.tsv"))
  write_tsv(reference$coverage, file.path(audit_dir, "reference_celltype_coverage.tsv"))
  write_tsv(reference_contribution, file.path(audit_dir, "reference_donor_contributions.tsv"))
  write_tsv(balance, file.path(audit_dir, "target_covariate_balance.tsv"))
  write_tsv(evaluation_strata, file.path(audit_dir, "evaluation_strata.tsv"))
  write_tsv(truth$gene_filter_table, file.path(audit_dir, "gene_filtering.tsv"))
  write_tsv(truth$markers, file.path(audit_dir, "reference_only_markers.tsv"))
  saveRDS(cfg, file.path(audit_dir, "settings.rds"), compress = cfg$rds_compress)
  writeLines(capture.output(utils::sessionInfo()), file.path(audit_dir, "sessionInfo.txt"))

  data.frame(
    run_label = run_label,
    reference_seed = seed,
    reference_sampling = cfg$reference_sampling,
    n_reference = length(reference$subjects),
    n_expression_targets = length(expression_targets),
    n_ctsde_targets = length(ctsde_subjects),
    n_genes = nrow(truth$pseudobulk_counts),
    reference_draw_attempt = reference$attempt,
    reference_subjects = collapse_ids(reference$subjects),
    expression_target_subjects = collapse_ids(expression_targets),
    ctsde_target_subjects = collapse_ids(ctsde_subjects),
    output_dir = seed_dir,
    stringsAsFactors = FALSE
  )
}

run_rosmap_benchmark_workflow <- function(cfg = rosmap_primary_config(), run_label = "primary") {
  prepared <- prepare_rosmap_benchmark(cfg)
  msg("9/9 Constructing the requested reference/target split(s).", cfg = cfg)

  summaries <- lapply(cfg$reference_seeds, function(seed) {
    run_rosmap_seed(prepared, seed = seed, run_label = run_label)
  })
  manifest <- do.call(rbind, summaries)
  write_tsv(manifest, file.path(cfg$output_dir, "run_manifest.tsv"))
  saveRDS(
    cfg,
    file.path(cfg$output_dir, "workflow_settings.rds"),
    compress = cfg$rds_compress
  )
  msg("ROSMAP workflow complete: ", cfg$output_dir, cfg = cfg)

  invisible(list(
    config = cfg,
    manifest = manifest,
    eligible_subjects = prepared$eligible_subjects,
    prepared = prepared
  ))
}


# Run automatically when sourced, unless explicitly disabled by the caller.
if (isTRUE(getOption("rosmap.primary.run_on_source", TRUE))) {
  rosmap_primary_results <- run_rosmap_benchmark_workflow()
}
