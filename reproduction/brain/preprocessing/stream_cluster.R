# Same original per-component functions and seeds; bounded component concurrency.
run_rosmap_step1_stream <- function(cfg) {
  assert_step1_packages();validate_step1_config(cfg)
  Sys.setenv(OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
  dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
  metadata <- read_step1_metadata(cfg); pool <- select_step1_subject_pool(metadata,cfg)
  signature<-list(settings=cfg,pool=pool,source_sizes=file.info(cfg$input_files)$size)
  state<-file.path(cfg$output_dir,'stream_settings.rds')
  if(file.exists(state))stopifnot(identical(readRDS(state),signature)) else saveRDS(signature,state)
  summaries<-list()
  for(component in names(cfg$input_files)) {
    output<-cfg$output_files[[component]];summary_file<-paste0(output,'.summary.rds')
    if(file.exists(output)&&file.exists(summary_file)) {
      message('Reusing completed component: ',component)
      summaries[[component]]<-readRDS(summary_file);next
    }
    if(file.exists(output))stop('Incomplete component output: ',output)
    message(format(Sys.time()),' Loading ',component)
    object<-load_and_subset_step1_component(component,pool,cfg)
    message(format(Sys.time()),' Clustering ',component,': ',ncol(object),' nuclei')
    object<-cluster_step1_component(object,component,cfg)
    summary<-step1_component_summary(setNames(list(object),component),pool,cfg)
    saveRDS(object,paste0(output,'.partial'),compress=cfg$rds_compress)
    stopifnot(file.rename(paste0(output,'.partial'),output))
    saveRDS(summary,summary_file);summaries[[component]]<-summary
    rm(object,summary);gc();message(format(Sys.time()),' Completed ',component)
  }
  summary<-do.call(rbind,summaries)
  table<-metadata[match(pool,metadata$subject),c('subject','condition'),drop=FALSE]
  table<-table[order(table$condition,table$subject),,drop=FALSE]
  data.table::fwrite(table,file.path(cfg$output_dir,'step1_subject_pool.tsv'),sep='\t')
  data.table::fwrite(summary,file.path(cfg$output_dir,'step1_component_summary.tsv'),sep='\t')
  saveRDS(list(settings=cfg,subject_pool=pool,output_files=cfg$output_files,component_summary=summary),
    file.path(cfg$output_dir,'step1_checkpoint_manifest.rds'))
  capture.output(sessionInfo(),file=file.path(cfg$output_dir,'step1_sessionInfo.txt'))
}
