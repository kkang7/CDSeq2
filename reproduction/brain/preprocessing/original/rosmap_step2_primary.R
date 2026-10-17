# ==============================================================================
# ROSMAP Step 2A: primary manuscript benchmark
#
# Prerequisite (run once):
#   source("~/data-lab5/Project_CDSeq2/code/rosmap_step1_cluster_once.R")
#
# Run this lightweight downstream workflow with:
#   source("~/data-lab5/Project_CDSeq2/code/rosmap_step2_primary.R")
# ==============================================================================

rosmap_step2_primary_source_file <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile)) return(normalizePath(ofile, mustWork = FALSE))
  }
  stop("Could not determine the Step 2 primary script path.", call. = FALSE)
}

.rosmap_step2_primary_file <- rosmap_step2_primary_source_file()
.rosmap_step2_common_file <- file.path(
  dirname(.rosmap_step2_primary_file),
  "rosmap_step2_common.R"
)
if (!file.exists(.rosmap_step2_common_file)) {
  stop("rosmap_step2_common.R was not found beside this script.", call. = FALSE)
}
source(.rosmap_step2_common_file, local = environment())


rosmap_step2_primary_config <- function() {
  cfg <- rosmap_step2_base_config()
  cfg$output_dir <- path.expand(
    paste0(
      "~/data-lab5/Project_CDSeq2/data/ROSMAP_AD/",
      "ref_heldout_benchmark_manuscript_primary"
    )
  )

  # Primary prespecification
  cfg$cluster_filter_enabled <- TRUE
  cfg$cluster_cols <- "RNA_snn_res.0.01"
  cfg$top_n_clusters <- 1L
  cfg$strict_cluster_subject_intersection <- FALSE
  cfg$reference_sampling <- "balanced"
  cfg$reference_seeds <- 123L
  cfg$marker_augmentation <- FALSE
  cfg$add_logcounts_to_reference <- FALSE
  cfg$rds_compress <- FALSE
  cfg
}

run_rosmap_step2_primary <- function(cfg = rosmap_step2_primary_config()) {
  run_rosmap_step2(cfg = cfg, run_label = "primary")
}


if (isTRUE(getOption("rosmap.step2.primary.run_on_source", TRUE))) {
  rosmap_step2_primary_results <- run_rosmap_step2_primary()
}
