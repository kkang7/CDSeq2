# Run from the repository root: Rscript --vanilla examples/run_tutorial.R
args <- commandArgs(TRUE)
file_arg <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
here <- dirname(normalizePath(file_arg))
out <- if(length(args))args[1]else file.path(here,'output')
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop('Choose a new empty output directory; completed examples are not overwritten.')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(requireNamespace('CDSeq2',quietly=TRUE))
for(f in c('01_simulate.R','02_deconvolve.R','03_prioritize.R','04_associations.R'))source(file.path(here,f))
started <- proc.time()[['elapsed']]
data <- simulate_tutorial()
message('Fitting 20 synthetic training donors')
training <- fit_tutorial(data$training,data,101)
message('Fitting 20 independent synthetic validation donors')
validation <- fit_tutorial(data$validation,data,202)
message('Training and validating the six-feature prioritization model')
priority <- prioritize_tutorial(data,training,validation)
CDSeq2::save_variation_scorer(priority$scorer,file.path(out,'toy_scorer.rds'))
reloaded <- CDSeq2::load_variation_scorer(file.path(out,'toy_scorer.rds'))
CDSeq2::check_variation_scorer(reloaded)
saveRDS(list(data=data,training=training,validation=validation,priority=priority),file.path(out,'tutorial_objects.rds'))
write.csv(priority$selected,file.path(out,'profile_scores.csv'),row.names=FALSE)
write.csv(priority$evaluation$donor_checks,file.path(out,'donor_separation.csv'),row.names=FALSE)
for(nm in c('phenotype','genotype'))write.csv(association_tutorial(validation$counts,priority$selected,
  data$validation[[nm]]),file.path(out,paste0(nm,'_associations.csv')),row.names=FALSE)
pr <- priority$evaluation$predictions
key <- function(x)paste(x$gene,x$cell_type,sep='|')
pr$retained <- priority$selected$retained[match(key(pr),key(priority$selected))]
recovery <- do.call(rbind,lapply(unique(pr$cell_type),function(ct){
  d<-pr[pr$cell_type==ct,];data.frame(cell_type=ct,eligible=sum(is.finite(d$variation_score)),
    retained=sum(d$retained),fraction=mean(d$retained),
    median_all=median(d$observed_r,na.rm=TRUE),median_selected=median(d$observed_r[d$retained],na.rm=TRUE))}))
write.csv(recovery,file.path(out,'recovery_summary.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
writeLines(sprintf('Elapsed seconds: %.1f',proc.time()[['elapsed']]-started),file.path(out,'run_summary.txt'))
print(recovery)
message('Tutorial completed: ',normalizePath(out))
