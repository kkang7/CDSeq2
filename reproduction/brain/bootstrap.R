argv <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
bundle <- dirname(normalizePath(script, mustWork=TRUE))
config_file <- if(length(argv)) argv[1] else file.path(bundle,"config.R")
source(config_file)
stopifnot(is.numeric(workers),length(workers)==1L,is.finite(workers),workers>=1,workers==as.integer(workers))
m <- read.csv(file.path(bundle,"engine_checksums.csv"),stringsAsFactors=FALSE)
stopifnot(identical(unname(tools::md5sum(file.path(bundle,m$file))),m$md5))
output_root <- path.expand(output_root)
direct <- file.path(output_root,"direct")
gam6 <- file.path(direct,"gam6_no_composition_loo_r2")
continuous <- file.path(direct,"continuous_associations_v1")
eqtl <- file.path(direct,"eqtl_three_types_v1")
bp <- file.path(direct,"eqtl_bayesprism_comparison_v1")
engine <- function(name) source(file.path(bundle,"engines",name),local=.GlobalEnv)
initialize <- function() {
  signature <- list(training=normalizePath(path.expand(training_file),mustWork=TRUE),
    reference=normalizePath(path.expand(reference_file),mustWork=TRUE),
    transfer=normalizePath(path.expand(transfer_dir),mustWork=TRUE),engines=m)
  marker <- file.path(output_root,"BUNDLE_RUN.rds")
  if(file.exists(marker)) stopifnot(identical(readRDS(marker),signature)) else {
    if(dir.exists(output_root) && length(list.files(output_root,all.files=TRUE,no..=TRUE)))
      stop("Choose a new empty output_root; refusing to use an unrelated existing run.")
    dir.create(output_root,recursive=TRUE,showWarnings=FALSE)
    saveRDS(signature,marker)
  }
}
direct_action <- function(action) {
 initialize(); engine("direct.R")
 run_direct_aggregation(action=action,workers=workers,training_file=training_file,
  reference_file=reference_file,transfer_dir=transfer_dir,output_dir=direct,
  folds=5L,seed=20261009L,min_profiles=100L)
}
