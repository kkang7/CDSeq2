script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"bootstrap.R"))
initialize(); engine("continuous.R")
# This audited engine requires a new output directory; never delete partial results automatically.
run_continuous_associations(run_dir=direct,gam6_dir=gam6,output_dir=continuous,
 programs=3L,program_size=10L,signal_fractions=c(.25,.5,.75),replicates=3L,
 null_replicates=200L,retentions=seq(.1,1,.1),seed=20261010L,
 covariates=NULL,bh_threshold=.05)
