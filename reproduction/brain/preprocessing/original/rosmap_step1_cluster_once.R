# ==============================================================================
# ROSMAP Step 1: create the reusable clustered-data checkpoint
#
# Run once with:
#   source("~/data-lab5/Project_CDSeq2/code/rosmap_step1_cluster_once.R")
#
# This expensive step:
#   1. selects a reproducible 100-AD + 100-Control starting pool;
#   2. loads and subsets each of the nine ROSMAP Seurat objects;
#   3. runs normalization, variable-feature selection, scaling, PCA, Harmony,
#      neighbor finding, and clustering at resolutions 0.05 and 0.01;
#   4. saves each clustered component separately for fast downstream loading.
#
# It intentionally does NOT choose or remove clusters. Cluster retention is a
# Step 2 decision so the primary, sensitivity, and package-example workflows
# can use the same expensive checkpoint.
#
# To load functions without running:
#   options(rosmap.step1.run_on_source = FALSE)
#   source("rosmap_step1_cluster_once.R")
# ==============================================================================

rosmap_step1_config <- function() {
  data_dir <- path.expand("~/data-lab5/Project_CDSeq2/data/ROSMAP_AD")
  output_dir <- file.path(data_dir, "clustered_once_balanced_200")

  input_files <- c(
    Ast      = file.path(data_dir, "Astrocytes.rds"),
    Excit_1  = file.path(data_dir, "Excitatory_neurons_set1.rds"),
    Excit_2  = file.path(data_dir, "Excitatory_neurons_set2.rds"),
    Excit_3  = file.path(data_dir, "Excitatory_neurons_set3.rds"),
    OPC      = file.path(data_dir, "OPCs.rds"),
    Imm      = file.path(data_dir, "Immune_cells.rds"),
    Oligo    = file.path(data_dir, "Oligodendrocytes.rds"),
    Inhib    = file.path(data_dir, "Inhibitory_neurons.rds"),
    Vascular = file.path(data_dir, "Vasculature_cells.rds")
  )

  output_files <- c(
    Ast      = file.path(output_dir, "Astrocytes.rds"),
    Excit_1  = file.path(output_dir, "Excitatory_neurons_set1.rds"),
    Excit_2  = file.path(output_dir, "Excitatory_neurons_set2.rds"),
    Excit_3  = file.path(output_dir, "Excitatory_neurons_set3.rds"),
    OPC      = file.path(output_dir, "OPCs.rds"),
    Imm      = file.path(output_dir, "Immune_cells.rds"),
    Oligo    = file.path(output_dir, "Oligodendrocytes.rds"),
    Inhib    = file.path(output_dir, "Inhibitory_neurons.rds"),
    Vascular = file.path(output_dir, "Vasculature_cells.rds")
  )

  list(
    input_files = input_files,
    output_files = output_files,
    metadata_file = file.path(data_dir, "individual_metadata_deidentified.tsv"),
    output_dir = output_dir,
    subject_col = "subject",
    metadata_subject_col = "subject",
    diagnosis_col = "Pathologic_diagnosis_of_AD",
    high_resolution_col = "cell_type_high_resolution",
    assay = "RNA",
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

    # Starting-pool design
    pool_n_per_condition = 100L,
    pool_seed = 123L,

    # Expensive clustering settings. Both resolutions are retained because
    # Step 2 sensitivity analyses use each of them.
    npcs = 50L,
    harmony_dims = 1:20,
    resolutions = c(0.05, 0.01),
    workers = 20L,
    limit_blas_threads_per_worker = TRUE,
    run_umap = FALSE,

    # FALSE gives faster downstream reads but larger files.
    rds_compress = FALSE,
    force_recluster = FALSE,
    verbose = TRUE
  )
}


# ==============================================================================
# Utilities
# ==============================================================================

step1_message <- function(..., cfg) {
  if (isTRUE(cfg$verbose)) message(...)
}

assert_step1_packages <- function() {
  packages <- c("Seurat", "SeuratObject", "harmony", "data.table")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Install the following packages before running Step 1: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

validate_step1_config <- function(cfg) {
  if (is.null(names(cfg$input_files)) || is.null(names(cfg$output_files))) {
    stop("input_files and output_files must be named vectors.", call. = FALSE)
  }
  if (!identical(names(cfg$input_files), names(cfg$output_files))) {
    stop("input_files and output_files must have identical names and order.", call. = FALSE)
  }
  missing_mapping <- setdiff(names(cfg$input_files), names(cfg$component_to_major))
  if (length(missing_mapping) > 0L) {
    stop(
      "component_to_major is missing: ", paste(missing_mapping, collapse = ", "),
      call. = FALSE
    )
  }
  missing_inputs <- cfg$input_files[!file.exists(cfg$input_files)]
  if (length(missing_inputs) > 0L) {
    stop(
      "Step 1 input files were not found:\n",
      paste(missing_inputs, collapse = "\n"),
      call. = FALSE
    )
  }
  if (!file.exists(cfg$metadata_file)) {
    stop("Metadata file was not found: ", cfg$metadata_file, call. = FALSE)
  }
  invisible(TRUE)
}

read_step1_metadata <- function(cfg) {
  meta <- data.table::fread(cfg$metadata_file, data.table = FALSE)
  required <- c(cfg$metadata_subject_col, cfg$diagnosis_col)
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0L) {
    stop("Metadata is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  meta$subject <- as.character(meta[[cfg$metadata_subject_col]])
  diagnosis <- tolower(trimws(as.character(meta[[cfg$diagnosis_col]])))
  if (!all(diagnosis %in% c("yes", "no"))) {
    stop("The pathological-AD diagnosis column must contain only yes/no.", call. = FALSE)
  }
  meta$condition <- ifelse(diagnosis == "yes", "AD", "Control")
  if (anyDuplicated(meta$subject)) {
    stop("Metadata contains duplicated subject IDs.", call. = FALSE)
  }
  meta
}

select_step1_subject_pool <- function(subject_metadata, cfg) {
  ad <- subject_metadata$subject[subject_metadata$condition == "AD"]
  control <- subject_metadata$subject[subject_metadata$condition == "Control"]
  n_each <- as.integer(cfg$pool_n_per_condition)
  if (length(ad) < n_each || length(control) < n_each) {
    stop(
      "Cannot select ", n_each, " subjects per condition: available AD = ",
      length(ad), ", Control = ", length(control), ".",
      call. = FALSE
    )
  }

  # This matches the original script: sample AD first and Control second from
  # the metadata using one shared seed.
  set.seed(cfg$pool_seed)
  sort(c(
    sample(ad, n_each, replace = FALSE),
    sample(control, n_each, replace = FALSE)
  ))
}

load_and_subset_step1_component <- function(component, subject_pool, cfg) {
  object <- readRDS(cfg$input_files[[component]])
  if (!inherits(object, "Seurat")) {
    stop(cfg$input_files[[component]], " is not a Seurat object.", call. = FALSE)
  }
  md <- object[[]]
  if (!cfg$subject_col %in% names(md)) {
    stop(component, " is missing subject column '", cfg$subject_col, "'.", call. = FALSE)
  }

  object$subject <- as.character(md[[cfg$subject_col]])
  keep <- rownames(object[[]])[as.character(object$subject) %in% subject_pool]
  if (length(keep) == 0L) {
    stop("No selected-pool cells were found in ", component, ".", call. = FALSE)
  }
  object <- object[, keep]
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
}

cluster_step1_component <- function(object, component, cfg) {
  set.seed(cfg$pool_seed)
  object <- Seurat::NormalizeData(object = object, assay = cfg$assay, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(
    object = object,
    assay = cfg$assay,
    verbose = FALSE
  )
  object <- Seurat::ScaleData(object = object, assay = cfg$assay, verbose = FALSE)

  npcs <- min(
    as.integer(cfg$npcs),
    ncol(object) - 1L,
    length(SeuratObject::VariableFeatures(object)) - 1L
  )
  if (!is.finite(npcs) || npcs < 2L) {
    stop("Too few cells or variable genes to cluster ", component, ".", call. = FALSE)
  }
  object <- Seurat::RunPCA(
    object = object,
    assay = cfg$assay,
    npcs = npcs,
    verbose = FALSE
  )

  dims <- cfg$harmony_dims[cfg$harmony_dims <= npcs]
  if (length(dims) < 2L) {
    stop("Too few PCs are available for ", component, ".", call. = FALSE)
  }
  # Harmony uses subject only. Diagnosis is not supplied to clustering.
  object <- harmony::RunHarmony(
    object = object,
    group.by.vars = cfg$subject_col,
    reduction.use = "pca",
    verbose = FALSE
  )
  object <- Seurat::FindNeighbors(
    object = object,
    reduction = "harmony",
    dims = dims,
    verbose = FALSE
  )
  for (resolution in cfg$resolutions) {
    object <- Seurat::FindClusters(
      object = object,
      resolution = resolution,
      random.seed = cfg$pool_seed,
      verbose = FALSE
    )
  }
  if (isTRUE(cfg$run_umap)) {
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

step1_component_summary <- function(components, subject_pool, cfg) {
  do.call(rbind, lapply(names(components), function(component) {
    object <- components[[component]]
    md <- object[[]]
    resolution_columns <- paste0("RNA_snn_res.", cfg$resolutions)
    missing_columns <- setdiff(resolution_columns, names(md))
    if (length(missing_columns) > 0L) {
      stop(
        component, " is missing expected cluster columns: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    data.frame(
      component = component,
      major_cell_type = unname(cfg$component_to_major[[component]]),
      n_cells = ncol(object),
      n_genes = nrow(object),
      n_subjects_observed = length(unique(as.character(object$subject))),
      n_pool_subjects = length(subject_pool),
      n_clusters_res_0_05 = length(unique(as.character(md$RNA_snn_res.0.05))),
      n_clusters_res_0_01 = length(unique(as.character(md$RNA_snn_res.0.01))),
      stringsAsFactors = FALSE
    )
  }))
}


# ==============================================================================
# Complete Step 1 workflow
# ==============================================================================

run_rosmap_step1_cluster_once <- function(cfg = rosmap_step1_config()) {
  assert_step1_packages()
  validate_step1_config(cfg)

  # With forked component-level parallelism, multithreaded BLAS inside every
  # worker can oversubscribe the machine. Request one BLAS/OpenMP thread per
  # worker; at most nine component workers are launched below.
  thread_variables <- c(
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
  )
  old_thread_values <- Sys.getenv(thread_variables, unset = NA_character_)
  if (isTRUE(cfg$limit_blas_threads_per_worker)) {
    new_thread_values <- stats::setNames(
      as.list(rep("1", length(thread_variables))),
      thread_variables
    )
    do.call(Sys.setenv, new_thread_values)
    on.exit({
      to_unset <- thread_variables[is.na(old_thread_values)]
      if (length(to_unset) > 0L) Sys.unsetenv(to_unset)
      to_restore <- !is.na(old_thread_values)
      if (any(to_restore)) {
        restore_values <- stats::setNames(
          as.list(old_thread_values[to_restore]),
          thread_variables[to_restore]
        )
        do.call(Sys.setenv, restore_values)
      }
    }, add = TRUE)
  }

  complete_checkpoint <- all(file.exists(cfg$output_files))
  partial_checkpoint <- any(file.exists(cfg$output_files)) && !complete_checkpoint
  if (complete_checkpoint && !isTRUE(cfg$force_recluster)) {
    step1_message(
      "Step 1 checkpoint already exists. No clustering was rerun: ",
      cfg$output_dir,
      cfg = cfg
    )
    return(invisible(list(
      skipped = TRUE,
      output_files = cfg$output_files,
      output_dir = cfg$output_dir
    )))
  }
  if (partial_checkpoint && !isTRUE(cfg$force_recluster)) {
    stop(
      "A partial Step 1 checkpoint exists in ", cfg$output_dir,
      ". Inspect it, then rerun with cfg$force_recluster <- TRUE if replacement is intended.",
      call. = FALSE
    )
  }

  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  step1_message("1/5 Reading metadata and selecting 100 AD + 100 Control subjects.", cfg = cfg)
  subject_metadata <- read_step1_metadata(cfg)
  subject_pool <- select_step1_subject_pool(subject_metadata, cfg)

  step1_message("2/5 Loading and subsetting the nine Seurat objects.", cfg = cfg)
  components <- lapply(names(cfg$input_files), function(component) {
    load_and_subset_step1_component(component, subject_pool, cfg)
  })
  names(components) <- names(cfg$input_files)
  invisible(gc())

  step1_message(
    "3/5 Running PCA, Harmony, neighbors, and both clustering resolutions.",
    cfg = cfg
  )
  workers <- min(as.integer(cfg$workers), length(components))
  if (.Platform$OS.type == "windows") workers <- 1L
  clustered <- parallel::mclapply(
    names(components),
    function(component) {
      cluster_step1_component(components[[component]], component, cfg)
    },
    mc.cores = workers,
    mc.preschedule = TRUE
  )
  names(clustered) <- names(components)
  failed <- vapply(clustered, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop(
      "Step 1 clustering failed for: ",
      paste(names(clustered)[failed], collapse = ", "),
      call. = FALSE
    )
  }
  rm(components)
  invisible(gc())

  step1_message("4/5 Validating and saving the reusable component files.", cfg = cfg)
  component_summary <- step1_component_summary(clustered, subject_pool, cfg)
  for (component in names(clustered)) {
    saveRDS(
      clustered[[component]],
      cfg$output_files[[component]],
      compress = cfg$rds_compress
    )
  }

  pool_table <- subject_metadata[
    match(subject_pool, subject_metadata$subject),
    c("subject", "condition"),
    drop = FALSE
  ]
  pool_table <- pool_table[order(pool_table$condition, pool_table$subject), , drop = FALSE]
  data.table::fwrite(
    data.table::as.data.table(pool_table),
    file.path(cfg$output_dir, "step1_subject_pool.tsv"),
    sep = "\t"
  )
  data.table::fwrite(
    data.table::as.data.table(component_summary),
    file.path(cfg$output_dir, "step1_component_summary.tsv"),
    sep = "\t"
  )
  saveRDS(
    list(
      settings = cfg,
      subject_pool = subject_pool,
      output_files = cfg$output_files,
      component_summary = component_summary
    ),
    file.path(cfg$output_dir, "step1_checkpoint_manifest.rds")
  )
  writeLines(
    capture.output(utils::sessionInfo()),
    file.path(cfg$output_dir, "step1_sessionInfo.txt")
  )

  step1_message("5/5 Step 1 complete: ", cfg$output_dir, cfg = cfg)
  invisible(list(
    skipped = FALSE,
    subject_pool = subject_pool,
    component_summary = component_summary,
    output_files = cfg$output_files,
    output_dir = cfg$output_dir
  ))
}


if (isTRUE(getOption("rosmap.step1.run_on_source", TRUE))) {
  rosmap_step1_results <- run_rosmap_step1_cluster_once()
}
