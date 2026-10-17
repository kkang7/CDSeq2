# Run from the repository root; source data remain external to the repository.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4L || !args[1]%in%c('verify-downloads','cluster','build'))
  stop('Usage: Rscript reproduction/brain/preprocessing/run.R verify-downloads|cluster|build DATA_DIR WORK_DIR OUTPUT_DIR')
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
here <- dirname(normalizePath(self)); source(file.path(here,'contract.R'))
data_dir <- normalizePath(path.expand(args[2]),mustWork=TRUE)
work <- normalizePath(path.expand(args[3]),mustWork=FALSE)
out <- normalizePath(path.expand(args[4]),mustWork=FALSE)
if(args[1]=='verify-downloads') {
  manifest<-read.csv(file.path(here,'original_download_md5.csv'),stringsAsFactors=FALSE)
  paths<-file.path(data_dir,manifest$file)
  stopifnot(all(file.exists(paths)))
  message('Reading all source files to verify MD5; this scans approximately 174 GiB.')
  observed<-unname(tools::md5sum(paths))
  if(!identical(observed,manifest$md5))stop('Source bytes differ: ',paste(manifest$file[observed!=manifest$md5],collapse=', '))
  message('All nine source RDS files match the historical identities.');quit(save='no')
}
original<-read.csv(file.path(here,'original_code_md5.csv'),stringsAsFactors=FALSE)
stopifnot(identical(unname(tools::md5sum(file.path(here,'original',original$file))),original$md5))

if(work==data_dir || out==data_dir || work==out) stop('Use distinct data, work and output directories.')
options(rosmap.step1.run_on_source=FALSE,rosmap.step2.primary.run_on_source=FALSE,
        rosmap.primary.run_on_source=FALSE)
source(file.path(here,'original/rosmap_step1_cluster_once.R'))
cfg1 <- rosmap_step1_config()
cfg1$input_files <- setNames(file.path(data_dir,basename(cfg1$input_files)),names(cfg1$input_files))
cfg1$output_dir <- file.path(work,'clustered_once_balanced_200')
cfg1$output_files <- setNames(file.path(cfg1$output_dir,basename(cfg1$output_files)),names(cfg1$output_files))
cfg1$metadata_file <- file.path(data_dir,'individual_metadata_deidentified.tsv')
cfg1$workers <- 1L # Bound memory: each source object can occupy tens of GB.
if(args[1]=='cluster') {
  run_rosmap_step1_cluster_once(cfg1)
  quit(save='no')
}
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop('OUTPUT_DIR must be empty; existing results are never overwritten.')
source(file.path(here,'original/rosmap_step2_primary.R'))
cfg <- rosmap_step2_primary_config()
cfg$input_files <- cfg1$output_files
cfg$step1_checkpoint_dir <- cfg1$output_dir
cfg$metadata_file <- cfg1$metadata_file
cfg$recluster_enabled <- FALSE; cfg$reuse_cluster_checkpoint <- FALSE
cfg$save_cluster_checkpoint <- FALSE; cfg$reuse_legacy_clustered_objects <- FALSE
cfg$output_dir <- file.path(work,'primary_preparation')
assert_rosmap_step1_checkpoint(cfg)
pr <- prepare_rosmap_benchmark(cfg)
d <- read.csv(file.path(here,'donors.csv'),stringsAsFactors=FALSE)
genes <- readLines(file.path(here,'genes.txt')); types <- cfg$major_celltypes
train <- d$donor[d$role=='training']; refs <- d$donor[d$role=='reference']
eval <- d[d$role=='evaluation',,drop=FALSE]
stopifnot(!anyDuplicated(d$donor),length(train)==40L,length(refs)==20L,nrow(eval)==40L,
          all(c(train,refs)%in%pr$eligible_subjects))
byct <- setNames(lapply(types,function(k) {
  x <- pr$major_data[[k]]
  aggregate_counts_by_subject(x$counts[genes,,drop=FALSE],x$metadata$subject,train)
}),types)
truth <- setNames(lapply(seq_along(train),function(i) {
  m <- do.call(cbind,lapply(byct,function(x)x$counts[,i,drop=FALSE])); colnames(m)<-types;m
}),train)
a <- array(0,c(length(genes),length(types),length(train)),list(genes,types,train))
for(i in seq_along(train)) a[,,i] <- as.matrix(truth[[i]])
cells <- do.call(cbind,lapply(byct,`[[`,'n_cells'));dimnames(cells)<-list(train,types)
md <- as.data.frame(pr$subject_metadata);md<-md[match(train,md$subject),,drop=FALSE]
md$sample<-md$donor<-train
sampled <- sample_reference_cells(pr$major_data,refs,cfg,20260919L)
rc <- sampled$counts[genes,,drop=FALSE];rmd<-as.data.frame(sampled$coldata)
rmd<-rmd[colnames(rc),,drop=FALSE]
rmd$condition<-as.character(pr$subject_metadata$condition[match(rmd$subject,pr$subject_metadata$subject)])
keep <- Matrix::colSums(rc)>0;rc<-rc[,keep,drop=FALSE];rmd<-rmd[keep,,drop=FALSE]
sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=rc),colData=rmd)
b <- list(gene_panel=genes,major_celltypes=types,pseudobulk_counts=apply(a,c(1,3),sum),
  truth_expression_counts=truth,truth_cell_counts=cells,target_metadata=md,reference_subjects=refs)
pooled <- sapply(types,function(k)as.numeric(Matrix::rowSums(rc[,rmd$cellType==k,drop=FALSE])))
dimnames(pooled)<-list(genes,types)
rm(pr,sampled,rc,rmd,byct,a);gc()
# Evaluation cohort uses ALL source nuclei for the frozen donors, without the
# training cohort's largest-cluster restriction. Do not substitute clustered files.
canon <- function(x) {
  x<-toupper(trimws(as.character(x)));ok<-!is.na(x)&grepl('^(ROSMAP[-_]?)?[0-9]+$',x)
  z<-sub('^ROSMAP[-_]?','',x[ok]);z<-sub('^0+(?=.)','',z,perl=TRUE)
  x[ok]<-paste0('ROSMAP-',z);if(anyNA(x)||any(!nzchar(x)))stop('Missing donor IDs');x
}
vcounts<-array(0,c(length(genes),length(types),nrow(eval)),list(genes,types,eval$sample))
ncells<-matrix(0,length(types),nrow(eval),dimnames=list(types,eval$sample))
for(component in names(cfg1$input_files)) {
  message('Aggregating evaluation nuclei: ',component)
  obj<-readRDS(cfg1$input_files[[component]])
  metadata<-obj[[]]; counts<-extract_count_matrix(obj,'RNA')
  stopifnot(identical(rownames(metadata),colnames(counts)),all(genes%in%rownames(counts)))
  donor<-canon(metadata$subject);k<-match(cfg1$component_to_major[[component]],types)
  for(i in seq_len(nrow(eval))) {
    take<-which(donor==eval$donor[i]);ncells[k,i]<-ncells[k,i]+length(take)
    if(length(take))vcounts[,k,i]<-vcounts[,k,i]+as.numeric(Matrix::rowSums(counts[genes,take,drop=FALSE]))
  }
  rm(obj,metadata,counts);gc()
}
v<-list(counts=vcounts,metadata=data.frame(sample=eval$sample,donor=eval$donor),
        n_cells=ncells,pooled_reference=pooled,reference_donors=refs)
req<-list(reference_mean=sweep(pooled,2,colSums(pooled),'/'),reference_donors=refs,
          fit=dget(file.path(here,'fit.R')))
# Fail closed if source revisions or clustering changes alter the paper inputs.
rosmap_verify_objects(b,sce,v,req,file.path(here,'expected_content.csv'))
dir.create(file.path(out,'TRANSFER_TO_MAC_FRESH_V1/inputs'),recursive=TRUE)
dir.create(file.path(out,'TRANSFER_TO_MAC_FRESH_V1/FRESH_REQUEST'),recursive=TRUE)
paths<-c('rosmap_heldout_benchmark_minimal.rds','rosmap_reference_sce.rds',
 'TRANSFER_TO_MAC_FRESH_V1/inputs/new_donor_scaffold.rds','TRANSFER_TO_MAC_FRESH_V1/FRESH_REQUEST/request.rds')
objects<-list(b,sce,v,req)
for(i in seq_along(paths))saveRDS(objects[[i]],file.path(out,paths[i]),version=2)
write.csv(data.frame(file=paths,md5=unname(tools::md5sum(file.path(out,paths)))),
          file.path(out,'generated_files.csv'),row.names=FALSE)
writeLines(c(paste0('training_file <- ',encodeString(file.path(out,paths[1]),quote='"')),
 paste0('reference_file <- ',encodeString(file.path(out,paths[2]),quote='"')),
 paste0('transfer_dir <- ',encodeString(file.path(out,'TRANSFER_TO_MAC_FRESH_V1'),quote='"')),
 paste0('output_root <- ',encodeString(paste0(out,'_analysis'),quote='"')),'workers <- 1L'),file.path(out,'brain_config.R'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
writeLines('Numerical identity verified against manuscript inputs.',file.path(out,'VERIFIED.txt'))
message('Verified inputs saved in ',out,'. Run preflight with ',file.path(out,'brain_config.R'))
