script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"bootstrap.R"))
initialize(); engine("step2_three_type_simulation.R")
action <- if(length(argv)>=2L) argv[2] else "run"
run_eqtl_three_types(action=action,workers=workers,run_dir=direct,
 gam6_file=file.path(gam6,"GAM6_model.rds"),output_dir=eqtl,
 target_types=c("Excit","Inhib","Oligo"),targets_per_quintile=10L,
 shared_per_pair=8L,offtype_per_type=4L,null_genes=40L,
 signal_ratios=c(0,.5,.8),mafs=.3,replicates=3L,seed=20261015L,
 expressed_fraction=.8,max_abs_beta=8,max_source_logCPM_sd=2,
 max_source_panel_share=.02,exclude_files=file.path(bundle,"excluded_previous_genes.csv"))
