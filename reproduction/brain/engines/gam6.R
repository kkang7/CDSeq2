# GAM6 ablation v1.0.0: remove ONLY composition_loo_r2 from original direct GAM7.
# Source defines functions only. No CDSeq2 fitting or package modification.
# source(file.choose())
# run_gam6_ablation()
# Default input: ~/Downloads/CDSeq2_direct_aggregation_native_v1
# Output: <input>/gam6_no_composition_loo_r2
# Uses the exact saved training rows, five gene folds, transformations, k=5,
# shrinkage cubic smooths, Gaussian REML, Fisher response and range-flag policy.
# Smoothing parameters are re-estimated on training rows, as required for GAM6.
# Validation is already examined: this is exploratory model development.
# Existing GAM7 files are read-only. Checkpoints allow resuming this ablation.
# Native R/mgcv execution has not been tested in the assistant environment.

.g6_script <- local({
  p<-unlist(lapply(sys.frames(),function(z)z$ofile),use.names=FALSE)
  if(length(p))normalizePath(tail(p,1),mustWork=TRUE)else NA_character_
})
.g6_assert <- function(ok,...)if(!isTRUE(ok))stop(...,call.=FALSE)
.g6_key <- function(d)paste(d$gene,d$cell_type,sep='\r')

# Launch the scoring experiment in a fresh Rscript session to avoid the previous
# RStudio-versus-worker BLAS/source-package mismatch. GAM fits are serial.
run_gam6_ablation <- function(run_dir='~/Downloads/CDSeq2_direct_aggregation_native_v1',
                              output_dir=NULL) {
  .g6_assert(!is.na(.g6_script)&&file.exists(.g6_script),'Source this saved script first.')
  run_dir<-normalizePath(path.expand(run_dir),mustWork=TRUE)
  if(is.null(output_dir))output_dir<-file.path(run_dir,'gam6_no_composition_loo_r2')
  output_dir<-path.expand(output_dir)
  job<-tempfile(fileext='.rds');launcher<-tempfile(fileext='.R')
  on.exit(unlink(c(job,launcher)),add=TRUE)
  saveRDS(list(run_dir=run_dir,output_dir=output_dir),job)
  writeLines(c(paste0('source(',encodeString(.g6_script,quote='"'),')'),
    paste0('do.call(.g6_run,readRDS(',encodeString(job,quote='"'),'))')),launcher)
  status<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(launcher)))
  .g6_assert(status==0L,'GAM6 ablation stopped; see the diagnostic above. Completed ablation fits are preserved.')
  invisible(output_dir)
}

.g6_run <- function(run_dir,output_dir) {
  .g6_assert(requireNamespace('mgcv',quietly=TRUE),'Install mgcv.')
  workflow<-file.path(run_dir,'workflow.R')
  .g6_assert(file.exists(workflow),'Original workflow.R is missing from run_dir.')
  # Original source only defines functions. Isolated environments prevent changing
  # any original global feature definitions or the saved GAM7 model.
  e7<-new.env(parent=globalenv());sys.source(workflow,e7)
  load<-e7$.da_load;save<-e7$.da_save;hash<-e7$.da_hash
  experiment<-load(file.path(run_dir,'experiment.rds'))
  .g6_assert(identical(unname(tools::md5sum(workflow)),experiment$config$script_md5),
    'Original workflow bytes differ from the completed experiment.')
  a<-file.path(run_dir,'analysis')
  tr<-load(file.path(a,'training_features_targets.rds'))
  va<-load(file.path(a,'validation_features_targets.rds'))
  old<-load(file.path(a,'variation_scorer.rds'))
  saved_validation<-load(file.path(a,'validation_results.rds'))$predictions
  fold_map<-load(file.path(a,'fold_by_gene.rds'))
  .g6_assert(identical(old$features,e7$.da_features)&&identical(old$method,'GAM7_no_theta'),
    'Expected the original direct GAM7 model.')
  runtime<-list(R=R.version.string,platform=R.version$platform,
    mgcv=as.character(utils::packageVersion('mgcv')),BLAS=unname(extSoftVersion()['BLAS']))
  previous<-old$provenance$runtime
  expected<-list(R=previous$R,platform=previous$platform,
    mgcv=unname(previous$versions['mgcv']),BLAS=previous$BLAS)
  different<-names(runtime)[!vapply(names(runtime),function(n)identical(runtime[[n]],expected[[n]]),TRUE)]
  .g6_assert(!length(different),'Use the same native R/mgcv runtime as the direct run. Differences: ',paste(different,collapse=', '))
  .g6_assert(!anyDuplicated(.g6_key(tr))&&!anyDuplicated(.g6_key(va))&&
    setequal(unique(tr$gene),names(fold_map))&&identical(fold_map,old$development$fold_by_gene),
    'Training identifiers or saved gene folds differ.')
  .g6_assert(setequal(unique(tr$gene),unique(va$gene))&&
    setequal(unique(tr$cell_type),unique(va$cell_type)),'Training/validation panel differs.')
  tr$fold<-unname(fold_map[tr$gene]);folds<-sort(unique(tr$fold));types<-names(old$models)
  .g6_assert(length(folds)==5L&&!anyNA(tr$fold),'Expected the five saved gene folds.')
  # Confirm original model prediction replay before any alternative fits.
  p7<-e7$predict_direct_variation(va,old)
  j<-match(.g6_key(p7),.g6_key(saved_validation))
  .g6_assert(!anyNA(j),'Saved validation keys are incomplete.')
  b<-saved_validation$predicted_r[j]
  .g6_assert(identical(is.na(p7$predicted_r),is.na(b))&&
    max(c(0,abs(p7$predicted_r-b)),na.rm=TRUE)<1e-8,'Original GAM7 validation replay failed.')
  e6<-new.env(parent=globalenv());sys.source(workflow,e6)
  e6$.da_features<-setdiff(e7$.da_features,'composition_loo_r2')
  .g6_assert(identical(e6$.da_logs,e7$.da_logs)&&length(e6$.da_features)==6L,'Unexpected transformation change.')
  # Fixed training eligibility prevents removal of one predictor from changing
  # the training population as well as the model specification.
  eligible<-is.finite(tr$observed_r)&rowSums(!is.finite(as.matrix(tr[e7$.da_features])))==0L
  minimum<-experiment$config$min_profiles
  .g6_assert(length(minimum)==1L&&is.finite(minimum)&&minimum>=50,'Missing original minimum support setting.')
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  output_dir<-normalizePath(output_dir,mustWork=TRUE)
  .g6_assert(output_dir!=normalizePath(run_dir,mustWork=TRUE)&&output_dir!=normalizePath(a,mustWork=TRUE),
    'Use a separate ablation output directory.')
  lock<-file.path(output_dir,'RUNNING.lock')
  .g6_assert(dir.create(lock,showWarnings=FALSE),'Ablation lock exists. Confirm no old process is running before removing a stale lock: ',lock)
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  writeLines(paste('PID',Sys.getpid()),file.path(lock,'owner.txt'))
  dir.create(file.path(output_dir,'fits'),showWarnings=FALSE)
  spec<-list(version='GAM6_ablation_1.0.0',source_model=old$model_id,
    original_workflow_md5=experiment$config$script_md5,
    ablation_script_md5=unname(tools::md5sum(.g6_script)),training=hash(tr),
    validation=hash(va),folds=hash(fold_map),eligibility=hash(eligible),runtime=runtime,
    features=e6$.da_features,min_profiles=minimum)
  id<-hash(spec);manifest<-file.path(output_dir,'experiment.rds')
  if(file.exists(manifest)).g6_assert(identical(load(manifest)$spec,spec),'Ablation context differs; choose a new output_dir.')
  else save(list(id=id,spec=spec),manifest)
  models<-list();oof<-list();audit<-list()
  message('Original GAM7 replay passed. Fitting 35 CV + 7 full GAM6 models; no CDSeq2 fitting.')
  for(ct in types)for(f in c(folds,NA_integer_)) {
    ti<-which(tr$cell_type==ct&eligible&(is.na(f)|tr$fold!=f))
    vi<-if(is.na(f))integer()else which(tr$cell_type==ct&tr$fold==f)
    .g6_assert(!length(intersect(tr$gene[ti],tr$gene[vi])),'Gene-fold leakage.')
    label<-if(is.na(f))'full'else paste0('fold',f)
    path<-file.path(output_dir,'fits',paste0(ct,'_',label,'.rds'))
    expected_n<-old$development$fit_audit
    audit_label<-if(is.na(f))'full'else as.character(f)
    # Original direct workflow labels its audit folds fold1, ..., fold5.
    row<-expected_n[expected_n$cell_type==ct&expected_n$fold==label,,drop=FALSE]
    .g6_assert(nrow(row)==1L&&row$training_rows==length(ti)&&row$training_genes==length(unique(tr$gene[ti])),
      'Training population differs from original GAM7 for ',ct,' ',label)
    if(file.exists(path)) {
      z<-load(path);.g6_assert(identical(z$id,id)&&identical(z$cell_type,ct)&&identical(z$label,label),'Wrong GAM6 cache.')
      m<-z$model
    } else {
      message(ct,' | ',label)
      m<-e6$.da_fit_gam(tr[ti,,drop=FALSE],minimum)
      save(list(id=id,cell_type=ct,label=label,model=m),path)
    }
    audit[[length(audit)+1L]]<-data.frame(cell_type=ct,fold=label,
      training_rows=m$n_training,training_genes=m$n_distinct_genes)
    if(is.na(f))models[[ct]]<-m else {
      p<-e6$predict_direct_variation(tr[vi,,drop=FALSE],list(features=e6$.da_features,models=setNames(list(m),ct)))
      oof[[length(oof)+1L]]<-e6$.da_join_predictions(tr[vi,,drop=FALSE],p)
    }
  }
  oof<-do.call(rbind,oof);rownames(oof)<-NULL
  model<-list(schema='CDSeq2_GAM6_no_composition_loo_ablation_v1',method='GAM6_no_composition_loo_r2',
    features=e6$.da_features,log_features=e6$.da_logs,models=models,gene_panel=old$gene_panel,
    context=old$context,donor_ids=old$donor_ids,reference_donor_ids=old$reference_donor_ids,
    oracle=old$oracle,source_model_id=old$model_id,model_id=paste0('GAM6_',id),
    status='exploratory_after_validation_review',runtime=runtime,fold_by_gene=fold_map,
    fit_audit=do.call(rbind,audit))
  p6<-e6$predict_direct_variation(va,model)
  model$sentinel<-list(features=va[seq_len(min(20L,nrow(va))),,drop=FALSE],predictions=p6[seq_len(min(20L,nrow(p6))),,drop=FALSE])
  save(model,file.path(output_dir,'GAM6_model.rds'))
  reloaded<-load(file.path(output_dir,'GAM6_model.rds'))
  .g6_assert(identical(p6,e6$predict_direct_variation(va,reloaded)),'GAM6 save/load prediction replay failed.')
  validation6<-e6$.da_join_predictions(va,p6)
  # Save full available predictions, then enforce a common finite-score universe
  # for GAM7/GAM6/baseline retention comparisons. Truth does not define eligibility.
  utils::write.csv(validation6,file.path(output_dir,'validation_GAM6_predictions.csv'),row.names=FALSE,na='')
  utils::write.csv(oof,file.path(output_dir,'training_OOF_GAM6_predictions.csv'),row.names=FALSE,na='')
  summary<-curves<-list()
  for(stage in c('training_OOF','validation')) {
    seven<-if(stage=='validation')e7$.da_join_predictions(va,p7)else old$development$predictions
    six<-if(stage=='validation')validation6 else oof
    j<-match(.g6_key(seven),.g6_key(six));.g6_assert(!anyNA(j),'Comparison keys differ.');six<-six[j,,drop=FALSE]
    common<-is.finite(seven$predicted_r)&is.finite(six$predicted_r)&
      is.finite(seven$contribution_mean)&is.finite(seven$expected_log_mean_count)
    for(name in c('GAM7','GAM6')) {
      d<-if(name=='GAM7')seven else six
      d$predicted_r[!common]<-NA_real_;d$variation_score[!common]<-NA_real_
      met<-e7$.da_metrics(d)
      met$summary$stage<-stage;met$summary$model<-name
      met$retention$stage<-stage;met$retention$scorer_universe<-name
      met$retention$method[met$retention$method=='GAM7']<-name
      summary[[length(summary)+1L]]<-met$summary
      # Support-baseline curves are identical on this shared ranking universe.
      if(name=='GAM6')met$retention<-met$retention[met$retention$method=='GAM6',,drop=FALSE]
      curves[[length(curves)+1L]]<-met$retention
    }
  }
  summary<-do.call(rbind,summary);curves<-do.call(rbind,curves)
  utils::write.csv(summary,file.path(output_dir,'comparison_summary.csv'),row.names=FALSE,na='')
  utils::write.csv(curves,file.path(output_dir,'comparison_retention_curves.csv'),row.names=FALSE,na='')
  utils::write.csv(curves[abs(curves$retention-.2)<1e-10,],file.path(output_dir,'comparison_top20.csv'),row.names=FALSE,na='')
  utils::write.csv(model$fit_audit,file.path(output_dir,'GAM6_fit_audit.csv'),row.names=FALSE)
  save(list(model_id=model$model_id,validation=validation6,training_OOF=oof,summary=summary,curves=curves),
    file.path(output_dir,'ablation_results.rds'))
  capture.output(sessionInfo(),file=file.path(output_dir,'sessionInfo.txt'))
  writeLines(c('Completed GAM6 feature-removal experiment.',
    'Original helper code, fixed training eligibility, saved gene folds and GAM fitting specification reused.',
    'Only composition_loo_r2 removed; REML smoothing re-estimated on training folds/full training cohort.',
    'No deconvolution, input regeneration, validation fitting, clipping, threshold tuning, or original model replacement.',
    'Exploratory: validation donors were already reviewed before this feature ablation.',
    'GAM6_model.rds deliberately uses a separate schema; do not pass it to the GAM7-only package loader.',
    'Primary comparisons are comparison_summary.csv, comparison_top20.csv and comparison_retention_curves.csv.'),
    file.path(output_dir,'SUMMARY.txt'))
  writeLines(as.character(Sys.time()),file.path(output_dir,'COMPLETED'))
  message('COMPLETED: ',output_dir)
  print(summary[summary$stage=='validation',c('cell_type','model','MAE','rank_Spearman','range_flagged')],row.names=FALSE)
  invisible(output_dir)
}
