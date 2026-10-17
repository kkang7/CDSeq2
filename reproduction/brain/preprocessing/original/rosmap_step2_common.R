# ==============================================================================
# ROSMAP Step 2 common loader and configuration
#
# This file is sourced by the three lightweight Step 2 entry scripts:
#   rosmap_step2_primary.R
#   rosmap_step2_sensitivity.R
#   rosmap_step2_cdseq2_example.R
#
# Keep it beside rosmap_primary_benchmark_workflow.R. It loads that file as the
# downstream benchmark engine but explicitly disables Harmony/reclustering.
# ==============================================================================

rosmap_current_source_file <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile)) {
      return(normalizePath(ofile, mustWork = FALSE))
    }
  }
  stop("Could not determine the path of the sourced ROSMAP script.", call. = FALSE)
}

.rosmap_step2_common_file <- rosmap_current_source_file()
.rosmap_step2_code_dir <- dirname(.rosmap_step2_common_file)
.rosmap_step2_engine_file <- file.path(
  .rosmap_step2_code_dir,
  "rosmap_primary_benchmark_workflow.R"
)
if (!file.exists(.rosmap_step2_engine_file)) {
  stop(
    "The Step 2 benchmark engine was not found beside this file: ",
    .rosmap_step2_engine_file,
    call. = FALSE
  )
}

.old_primary_run_option <- getOption("rosmap.primary.run_on_source")
options(rosmap.primary.run_on_source = FALSE)
source(.rosmap_step2_engine_file, local = environment())
if (is.null(.old_primary_run_option)) {
  options(rosmap.primary.run_on_source = NULL)
} else {
  options(rosmap.primary.run_on_source = .old_primary_run_option)
}


rosmap_step1_component_files <- function() {
  data_dir <- path.expand(
    "~/data-lab5/Project_CDSeq2/data/ROSMAP_AD/clustered_once_balanced_200"
  )
  c(
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
}

rosmap_step2_base_config <- function(active_components = NULL) {
  cfg <- rosmap_primary_config()
  step1_files <- rosmap_step1_component_files()
  if (!is.null(active_components)) {
    unknown <- setdiff(active_components, names(step1_files))
    if (length(unknown) > 0L) {
      stop(
        "Unknown active Step 1 components: ", paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    step1_files <- step1_files[active_components]
  }

  cfg$input_files <- step1_files
  cfg$cluster_filter_units <- intersect(cfg$cluster_filter_units, names(step1_files))

  # The defining Step 2 property: never run normalization/PCA/Harmony/clustering.
  cfg$recluster_enabled <- FALSE
  cfg$reuse_cluster_checkpoint <- FALSE
  cfg$save_cluster_checkpoint <- FALSE
  cfg$reuse_legacy_clustered_objects <- FALSE
  cfg$legacy_clustered_file <- NULL
  cfg$cluster_checkpoint_file <- NULL
  cfg$step1_checkpoint_dir <- dirname(step1_files[[1L]])
  cfg
}

assert_rosmap_step1_checkpoint <- function(cfg) {
  missing <- cfg$input_files[!file.exists(cfg$input_files)]
  if (length(missing) > 0L) {
    stop(
      "The reusable Step 1 clustered checkpoint is incomplete. Missing files:\n",
      paste(missing, collapse = "\n"),
      "\nRun rosmap_step1_cluster_once.R before this Step 2 workflow.",
      call. = FALSE
    )
  }

  manifest <- file.path(dirname(cfg$input_files[[1L]]), "step1_checkpoint_manifest.rds")
  if (!file.exists(manifest)) {
    warning(
      "The clustered component files exist, but the Step 1 manifest is missing: ",
      manifest,
      call. = FALSE
    )
  } else {
    checkpoint <- tryCatch(readRDS(manifest), error = function(e) NULL)
    if (is.null(checkpoint) || is.null(checkpoint$settings)) {
      stop("The Step 1 checkpoint manifest cannot be read or is incomplete.", call. = FALSE)
    }
    required_resolutions <- c(0.05, 0.01)
    saved_resolutions <- as.numeric(checkpoint$settings$resolutions)
    if (!all(required_resolutions %in% saved_resolutions)) {
      stop(
        "The Step 1 checkpoint does not contain both required clustering ",
        "resolutions 0.05 and 0.01.",
        call. = FALSE
      )
    }
    if (!identical(
      as.integer(checkpoint$settings$pool_seed),
      as.integer(cfg$pool_seed)
    )) {
      stop("The Step 1 pool seed does not match the Step 2 configuration.", call. = FALSE)
    }
    if (!identical(
      as.integer(checkpoint$settings$pool_n_per_condition),
      as.integer(cfg$pool_n_per_condition)
    )) {
      stop("The Step 1 pool size does not match the Step 2 configuration.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

run_rosmap_step2 <- function(cfg, run_label) {
  assert_rosmap_step1_checkpoint(cfg)
  if (isTRUE(cfg$recluster_enabled)) {
    stop("Step 2 must have recluster_enabled = FALSE.", call. = FALSE)
  }
  run_rosmap_benchmark_workflow(cfg = cfg, run_label = run_label)
}
