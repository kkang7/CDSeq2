# CDSeq2 direct aggregation experiment v1.0.0
# Source defines functions only; it never starts fits automatically.
# source("~/Downloads/CDSeq2_direct_aggregation.R")
# run_direct_aggregation(action="prepare")
# run_direct_aggregation(action="probe", workers=1)
# run_direct_aggregation(action="run", workers=4)
# Repeat the identical call to resume. Use a fresh R session, not devtools::load_all().
# Requires installed CDSeq2, Matrix, SingleCellExperiment, SummarizedExperiment, mgcv.
# See companion README for paths, interpretation and verification limitations.

.da_script <- local({
  p <- unlist(lapply(sys.frames(), function(z) z$ofile), use.names=FALSE)
  if(length(p)) normalizePath(tail(p,1),mustWork=TRUE) else NA_character_
})
.da_assert <- function(ok,...) if(!isTRUE(ok)) stop(...,call.=FALSE)
.da_md5 <- function(p) unname(tools::md5sum(p))
.da_hash <- function(x) {p<-tempfile();on.exit(unlink(p));saveRDS(x,p,version=2);.da_md5(p)}
.da_save <- function(x,p) {
  dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile(tmpdir=dirname(p));on.exit(unlink(tmp));saveRDS(x,tmp,version=2)
  .da_assert(file.rename(tmp,p),"Cannot save ",p)
  writeLines(.da_md5(p),paste0(p,".md5"));invisible(p)
}
.da_load <- function(p) {
  .da_assert(file.exists(p)&&file.exists(paste0(p,".md5"))&&
    identical(.da_md5(p),readLines(paste0(p,".md5"),warn=FALSE)),"Missing/changed checkpoint: ",p)
  readRDS(p)
}
.da_csv <- function(d,p) utils::write.csv(d,p,row.names=FALSE,na="")
.da_log <- function(root,...) {
  z<-paste(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),paste0(...,collapse=""))
  message(z);cat(z,"\n",file=file.path(root,"progress.log"),append=TRUE)
}
.da_seed <- function(seed,fun) {
  had<-exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  if(had) old<-get(".Random.seed",.GlobalEnv)
  kinds<-RNGkind()
  on.exit({do.call(RNGkind,as.list(kinds));if(had)assign(".Random.seed",old,.GlobalEnv)
    else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)})
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed);fun()
}
.da_cor <- function(a,b,method="pearson") {
  if(length(a)<3L||any(!is.finite(c(a,b)))||sd(a)<1e-12||sd(b)<1e-12)return(NA_real_)
  unname(stats::cor(a,b,method=method))
}
.da_ids <- function(x) !is.null(x)&&length(x)>0L&&!anyNA(x)&&all(nzchar(x))&&!anyDuplicated(x)
.da_array <- function(a) {
  .da_assert(is.array(a)&&length(dim(a))==3L&&is.numeric(a)&&all(is.finite(a))&&
    all(a>=0)&&all(a==floor(a))&&all(vapply(dimnames(a),.da_ids,TRUE)),"Invalid observed count array")
  .da_assert(all(apply(a,c(2,3),sum)>0),"A donor/type library is zero")
  a
}
.da_features <- c("contribution_mean","expected_log_mean_count","ref_specificity",
  "bulk_logCPM_sd","bulk_max_donor_share","composition_loo_r2","bulk_theta_r")
.da_logs <- c("contribution_mean","ref_specificity","bulk_logCPM_sd","bulk_max_donor_share")

# Same seven predictor definitions as variation interface 2.3.0.
# bulk_theta_r uses CPM, not logCPM. RNA proportions are estimated by CDSeq2.
.da_observables <- function(bulk,mu,theta) {
  genes<-rownames(bulk);types<-colnames(mu);samples<-colnames(bulk)
  .da_assert(setequal(rownames(mu),genes)&&setequal(rownames(theta),samples)&&
    setequal(colnames(theta),types),"Observable axes differ")
  mu<-mu[genes,types,drop=FALSE];theta<-theta[samples,types,drop=FALSE]
  .da_assert(all(colSums(bulk)>0)&&all(is.finite(theta))&&all(theta>=0)&&
    max(abs(rowSums(theta)-1))<1e-6&&max(abs(colSums(mu)-1))<1e-6,"Invalid normalized inputs")
  B<-sweep(bulk,2,colSums(bulk),"/")*1e6;logB<-log2(B+.5)
  X<-cbind(1,theta[,-ncol(theta),drop=FALSE]);q<-qr(X)
  residuals<-t(as.matrix(stats::lm.fit(X,t(logB))$residuals))
  tss<-rowSums((logB-rowMeans(logB))^2)
  h<-rowSums(qr.Q(q)[,seq_len(q$rank),drop=FALSE]^2)
  loo<-if(all(h<1-1e-8))1-rowSums(sweep(residuals,2,1-h,"/")^2)/pmax(tss,1e-12) else rep(NA_real_,length(genes))
  contribution<-allocated<-matrix(0,nrow(mu),ncol(mu),dimnames=dimnames(mu))
  for(i in seq_along(samples)) {
    x<-sweep(mu,2,theta[i,],"*");w<-x/pmax(rowSums(x),1e-300)
    contribution<-contribution+w;allocated<-allocated+w*bulk[,i]
  }
  contribution<-contribution/length(samples);allocated<-allocated/length(samples)
  d<-do.call(rbind,lapply(seq_along(types),function(k)data.frame(
    gene=genes,cell_type=types[k],cohort_n=length(samples),
    contribution_mean=contribution[,k],expected_log_mean_count=log2(allocated[,k]+.5),
    ref_specificity=mu[,k]/pmax(rowSums(mu),1e-300),bulk_logCPM_sd=apply(logB,1,sd),
    bulk_max_donor_share=apply(bulk,1,max)/pmax(rowSums(bulk),1),composition_loo_r2=loo,
    bulk_theta_r=apply(B,1,function(x).da_cor(x,theta[,k])),theta_mean=mean(theta[,k]),row.names=NULL)))
  attr(d,"modeled_gene_panel")<-genes;attr(d,"feature_schema")<-"CDSeq2_GAM7_features_v1";d
}
.da_targets <- function(truth,estimated) {
  .da_assert(identical(dimnames(truth),dimnames(estimated)),"Recovery arrays disagree")
  out<-list()
  for(k in seq_len(dim(truth)[2])) {
    a<-matrix(estimated[,k,],nrow=dim(truth)[1]);b<-matrix(truth[,k,],nrow=dim(truth)[1])
    ok<-all(colSums(a)>0&colSums(b)>0)
    r<-sp<-rep(NA_real_,nrow(a))
    if(ok) {
      a<-sweep(a,2,colSums(a),"/");b<-sweep(b,2,colSums(b),"/")
      r<-vapply(seq_len(nrow(a)),function(g).da_cor(a[g,],b[g,]),0)
      sp<-vapply(seq_len(nrow(a)),function(g).da_cor(a[g,],b[g,],"spearman"),0)
    }
    out[[k]]<-data.frame(gene=dimnames(truth)[[1]],cell_type=dimnames(truth)[[2]][k],
      observed_r=r,observed_spearman=sp,complete_type_libraries=ok)
  }
  do.call(rbind,out)
}
.da_x <- function(d) {
  x<-as.matrix(d[.da_features]);storage.mode(x)<-"double"
  for(nm in .da_logs)x[,nm]<-log2(x[,nm]+1e-6)
  x
}
.da_fit_gam <- function(d,min_profiles) {
  ok<-is.finite(d$observed_r)&rowSums(!is.finite(as.matrix(d[.da_features])))==0L
  d<-d[ok,,drop=FALSE]
  .da_assert(length(unique(d$gene))>=min_profiles,"Insufficient distinct eligible genes for GAM; no fallback fitted")
  x<-.da_x(d);reg<-as.data.frame(x);reg$.y<-atanh(pmax(-.995,pmin(.995,d$observed_r)))
  nu<-vapply(reg[colnames(x)],function(z)length(unique(z)),1L);active<-names(nu)[nu>1L]
  terms<-vapply(active,function(nm)if(nu[[nm]]<5L)nm else paste0('s(',nm,', k=5, bs="cs")'),"")
  s<-mgcv::s
  f<-stats::as.formula(paste(".y ~",if(length(terms))paste(terms,collapse=" + ")else"1"),env=environment())
  fit<-mgcv::gam(f,data=reg,family=stats::gaussian(),method="REML",control=mgcv::gam.control(nthreads=1L))
  .da_assert(isTRUE(fit$converged),"GAM failed to converge")
  if(!is.null(fit$outer.info$conv)).da_assert(identical(fit$outer.info$conv,"full convergence"),"Smoothing optimization failed")
  for(nm in c("formula","terms","pterms","pred.formula"))if(!is.null(fit[[nm]]))environment(fit[[nm]])<-asNamespace("mgcv")
  list(kind="gam",fit=fit,terms=terms,dropped=names(nu)[nu<=1L],lower=apply(x,2,min),upper=apply(x,2,max),
    n_training=nrow(d),n_distinct_genes=length(unique(d$gene)))
}

# Also available for future cohorts with matching feature definitions/context.
# Retains and flags finite out-of-range predictions, as in the prior experiment.
predict_direct_variation <- function(features,scorer) {
  .da_assert(identical(scorer$features,.da_features),"Unsupported scorer features")
  .da_assert(all(c("gene","cell_type",.da_features)%in%names(features)),"Missing features")
  out<-features[c("gene","cell_type")];n<-nrow(out)
  out$predicted_r<-out$variation_score<-rep(NA_real_,n)
  out$outside_training_range<-rep(NA,n);out$range_features<-rep("",n);out$status<-rep("undefined_feature",n)
  for(ct in unique(out$cell_type)) {
    ii<-which(out$cell_type==ct);m<-scorer$models[[ct]]
    if(is.null(m)){out$status[ii]<-"unsupported_cell_type";next}
    good<-rowSums(!is.finite(as.matrix(features[ii,.da_features,drop=FALSE])))==0L
    out$range_features[ii[!good]]<-"undefined_feature";jj<-ii[good];if(!length(jj))next
    x<-.da_x(features[jj,,drop=FALSE]);p<-tanh(as.numeric(stats::predict(m$fit,as.data.frame(x),type="response")))
    .da_assert(length(p)==length(jj)&&all(is.finite(p)),"Invalid GAM prediction")
    b<-sweep(x,2,m$lower-1e-9,"<")|sweep(x,2,m$upper+1e-9,">")
    out$predicted_r[jj]<-p;out$variation_score[jj]<-100*pmax(0,p)
    out$outside_training_range[jj]<-rowSums(b)>0
    out$range_features[jj]<-apply(b,1,function(z)paste(colnames(x)[z],collapse=";"))
    out$status[jj]<-ifelse(rowSums(b)>0,"outside_training_range","ok")
  }
  out$contribution_mean<-features$contribution_mean;out$expected_log_mean_count<-features$expected_log_mean_count
  out
}

.da_extract <-
function (raw, bulk, types, cfg) 
{
    sample <- colnames(bulk)
    genes <- rownames(bulk)
    .da_assert(is.list(raw$cdseq_result) && sample %in% names(raw$cdseq_result), "Missing sample in cdseq_result.")
    rr <- raw$cdseq_result[[sample]]
    a <- rr$estGEP_array
    dn <- dimnames(a)
    .da_assert(length(dim(a)) == 3L && length(dn) == 3L, "estGEP_array must retain allocated-count draws.")
    gi <- which(vapply(dn, function(x) !is.null(x) && setequal(x, genes), TRUE))
    ki <- which(vapply(dn, function(x) !is.null(x) && setequal(x, types), TRUE))
    .da_assert(length(gi) == 1 && length(ki) == 1 && gi != ki, "Cannot orient retained draws by names.")
    a <- aperm(a, c(gi, ki, setdiff(1:3, c(gi, ki))))
    a <- a[genes, types, , drop = FALSE]
    .da_assert(dim(a)[3] == cfg$fit$save_last_n_mc && all(is.finite(a) & a >= 0), "Wrong retained draws or invalid allocations.")
    error <- max(abs(sweep(apply(a, c(1, 3), sum), 1, bulk[, 1], "-")))
    .da_assert(error <= max(1e-06, max(bulk) * 1e-10), "Allocated counts do not reconstruct every bulk gene.")
    theta <- rr$estProp
    if (is.null(dim(theta))) {
        .da_assert(setequal(names(theta), types), "Unnamed estProp.")
        theta <- theta[types]
    }
    else {
        theta <- as.matrix(drop(theta))
        if (setequal(rownames(theta), types)) 
            theta <- rowMeans(theta[types, , drop = FALSE])
        else if (setequal(colnames(theta), types)) 
            theta <- colMeans(theta[, types, drop = FALSE])
        else stop("Cannot orient estProp.", call. = FALSE)
    }
    .da_assert(all(is.finite(theta) & theta >= 0) && abs(sum(theta) - 1) < 1e-06, "Invalid returned estProp.")
    beta <- raw$parameters$beta
    if (is.null(beta)) 
        beta <- raw$parameters$BETA
    .da_assert(!is.null(beta), "Returned reference prior BETA is required.")
    beta <- as.matrix(beta)
    if (setequal(rownames(beta), types) && setequal(colnames(beta), genes)) 
        beta <- t(beta)
    .da_assert(setequal(rownames(beta), genes) && setequal(colnames(beta), types), "BETA axis mismatch.")
    beta <- beta[genes, types, drop = FALSE]
    .da_assert(all(is.finite(beta) & beta > 0) && max(abs(colSums(beta)/cfg$fit$reference_precision - 1)) < 1e-06, 
        "Reference prior concentration changed.")
    counts <- apply(a, c(1, 2), mean)
    dimnames(counts) <- list(genes, types)
    list(counts = counts, proportions = setNames(as.numeric(theta), types), reference_mean = sweep(beta, 2, colSums(beta), 
        "/"), reconstruction_max_error = error, retained_draws = dim(a)[3], proportion_source = "estProp", count_source = "mean(estGEP_array)")
}


.da_runtime <- function() {
  pk<-c('CDSeq2','mgcv','Matrix','SingleCellExperiment','SummarizedExperiment')
  for(p in pk).da_assert(requireNamespace(p,quietly=TRUE),'Install required package: ',p)
  home<-find.package('CDSeq2')
  files<-sort(unique(c(file.path(home,'DESCRIPTION'),
    list.files(file.path(home,'R'),recursive=TRUE,full.names=TRUE),
    list.files(file.path(home,'libs'),recursive=TRUE,full.names=TRUE))))
  stamps<-tools::md5sum(files);names(stamps)<-substring(files,nchar(home)+2L)
  list(R=R.version.string,platform=R.version$platform,
    versions=vapply(pk,function(p)as.character(utils::packageVersion(p)),''),
    cdseq_files=stamps,BLAS=unname(extSoftVersion()['BLAS']),
    cdseq_formals=formals(getExportedValue('CDSeq2','cdseq2')))
}
.da_path <- function(p) normalizePath(path.expand(p),mustWork=TRUE)
.da_find_transfer <- function(path=NULL) {
  if(!is.null(path))return(.da_path(path))
  candidates<-c('~/Downloads/TRANSFER_TO_MAC_FRESH_V1',
    '~/Downloads/preparation_6057b67c1ccd_TRANSFER_TO_MAC/TRANSFER_TO_MAC_FRESH_V1',
    '~/Downloads/preparation_6057b67c1ccd_TRANSFER_TO_MAC(1)/TRANSFER_TO_MAC_FRESH_V1')
  hit<-candidates[dir.exists(path.expand(candidates))]
  .da_assert(length(hit)>0L,'Unzip the transfer archive, then supply transfer_dir pointing to TRANSFER_TO_MAC_FRESH_V1.')
  .da_path(hit[1])
}
.da_inputs <- function(training_file,reference_file,transfer_dir) {
  scaffold_file<-file.path(transfer_dir,'inputs','new_donor_scaffold.rds')
  request_file<-file.path(transfer_dir,'FRESH_REQUEST','request.rds')
  files<-c(training=training_file,reference=reference_file,scaffold=scaffold_file,request=request_file)
  .da_assert(all(file.exists(files)),'Missing input: ',paste(files[!file.exists(files)],collapse=', '))
  expected<-c(training='1bb46db5e1b8f69088dd273deb9c3ead',reference='2831486744ba826a824145398c9a6991',
    scaffold='9c341b5907dd16338365665ef3414a27')
  b<-readRDS(training_file);v<-readRDS(scaffold_file);req<-readRDS(request_file)
  # Accept a regenerated serialization only when its numerical content matches
  # the frozen manuscript inputs. The original RDS hashes remain the fast path.
  original <- identical(unname(tools::md5sum(files[names(expected)])),unname(expected)) &&
    identical(.da_md5(request_file),'6057b67c1ccdead0306fe14c20fb36e0')
  if(!original) {
    pre <- file.path(dirname(dirname(.da_script)),'preprocessing')
    ce <- new.env(parent=.GlobalEnv);sys.source(file.path(pre,'contract.R'),ce)
    ce$rosmap_verify_objects(b,readRDS(reference_file),v,req,file.path(pre,'expected_content.csv'))
  }
  genes<-as.character(b$gene_panel);types<-as.character(b$major_celltypes)
  samples<-colnames(b$pseudobulk_counts)
  .da_assert(length(genes)==2000L&&length(types)==7L&&length(samples)==40L,'Unexpected training dimensions')
  a<-array(0,c(length(genes),length(types),length(samples)),list(genes,types,samples))
  for(i in seq_along(samples)) {
    m<-as.matrix(b$truth_expression_counts[[samples[i]]])
    .da_assert(setequal(rownames(m),genes)&&setequal(colnames(m),types),'Training component axes differ')
    a[,,i]<-m[genes,types,drop=FALSE]
  }
  a<-.da_array(a);vb<-.da_array(v$counts)
  .da_assert(setequal(dimnames(vb)[[1]],genes)&&setequal(dimnames(vb)[[2]],types)&&dim(vb)[3]==40L,
    'Validation panel/types differ')
  vb<-vb[genes,types,,drop=FALSE]
  .da_assert(identical(unname(apply(a,c(1,3),sum)),unname(as.matrix(b$pseudobulk_counts[genes,samples,drop=FALSE]))),
    'Original training bulk does not equal direct component sum')
  tm<-b$target_metadata;vm<-v$metadata
  for(md in list(tm,vm)).da_assert(all(c('sample','donor')%in%names(md))&&.da_ids(as.character(md$sample))&&
    .da_ids(as.character(md$donor)),'Require one distinct biological donor per sample')
  .da_assert(setequal(tm$sample,samples)&&setequal(vm$sample,dimnames(vb)[[3]]),'Metadata/sample axes differ')
  tm<-tm[match(samples,tm$sample),,drop=FALSE];vm<-vm[match(dimnames(vb)[[3]],vm$sample),,drop=FALSE]
  train_ids<-as.character(tm$donor);val_ids<-as.character(vm$donor);ref_ids<-as.character(b$reference_subjects)
  .da_assert(length(ref_ids)==20L&&.da_ids(ref_ids)&&setequal(ref_ids,v$reference_donors)&&
    setequal(ref_ids,req$reference_donors),'Reference donor set differs')
  .da_assert(!length(intersect(train_ids,val_ids))&&!length(intersect(train_ids,ref_ids))&&
    !length(intersect(val_ids,ref_ids)),'Biological donor leakage detected')
  .da_assert(length(unique(c(train_ids,val_ids,ref_ids)))==100L,'Expected 100 distinct biological donors')
  sce<-readRDS(reference_file)
  .da_assert(inherits(sce,'SingleCellExperiment'),'Reference must be a SingleCellExperiment')
  rc<-SummarizedExperiment::assay(sce,'counts');md<-as.data.frame(SummarizedExperiment::colData(sce))
  .da_assert(setequal(rownames(rc),genes)&&identical(rownames(md),colnames(rc))&&
    all(c('subject','cellType')%in%names(md))&&setequal(md$subject,ref_ids)&&
    setequal(md$cellType,types)&&ncol(rc)==39692L,'Reference axes/donors/nuclei differ')
  rc<-rc[genes,,drop=FALSE];ann<-as.character(md$cellType)
  pooled<-sapply(types,function(ct)as.numeric(Matrix::rowSums(rc[,ann==ct,drop=FALSE])))
  dimnames(pooled)<-list(genes,types)
  .da_assert(max(abs(pooled-v$pooled_reference[genes,types]))==0,'Pooled support reference differs')
  mu<-sweep(pooled,2,colSums(pooled),'/')
  .da_assert(max(abs(mu-req$reference_mean[genes,types]))<1e-12,'Feature reference mean differs')
  list(training=list(truth=a,bulk=apply(a,c(1,3),sum),metadata=tm,donors=train_ids),
    validation=list(truth=vb,bulk=apply(vb,c(1,3),sum),metadata=vm,donors=val_ids),
    reference=list(counts=as.matrix(rc),annotation=ann),reference_mean=mu,reference_donors=ref_ids,
    fit=req$fit,checksums=tools::md5sum(files),training_cells=b$truth_cell_counts,validation_cells=v$n_cells)
}
.da_worker <- function(job) {
  w<-get('.da_worker_state',.GlobalEnv);p<-file.path(w$root,'sample_fits',paste0(job$cohort,'_',sprintf('%03d',job$i),'.rds'))
  tryCatch({
    bulk<-w$bulks[[job$cohort]][,job$i,drop=FALSE]
    sig<-list(run_id=w$id,cohort=job$cohort,sample=colnames(bulk),seed=job$seed,bulk_hash=.da_hash(bulk))
    if(file.exists(p)) {
      z<-.da_load(p);.da_assert(identical(z$signature,sig),'Sample cache identity differs')
      return(data.frame(cohort=job$cohort,sample=colnames(bulk),status='cached',seconds=z$seconds,error=''))
    }
    args<-c(w$fit,list(bulk_data=bulk,reference_scRNA_count=w$ref$counts,
      reference_scRNA_annotation=w$ref$annotation,cpu_number=1L))
    optional<-list(seed=job$seed,n_chains=1L,return_mcmc_arrays=TRUE,return_assignments=FALSE,
      verbose=FALSE,drop_all_zero_genes=FALSE)
    fn<-getExportedValue('CDSeq2','cdseq2')
    args<-c(args,optional[intersect(names(optional),names(formals(fn)))])
    message('START ',job$cohort,' ',colnames(bulk),' | PID ',Sys.getpid())
    start<-proc.time()[3];warnings<-character()
    raw<-.da_seed(job$seed,function()withCallingHandlers(do.call(fn,args),warning=function(e){
      warnings<<-c(warnings,conditionMessage(e));invokeRestart('muffleWarning')}))
    result<-.da_extract(raw,bulk,w$types,list(fit=w$fit))
    z<-list(signature=sig,result=result,seconds=unname(proc.time()[3]-start),warnings=unique(warnings))
    .da_save(z,p);message('FINISHED ',job$cohort,' ',colnames(bulk))
    data.frame(cohort=job$cohort,sample=colnames(bulk),status='fitted',seconds=z$seconds,error='')
  },error=function(e)data.frame(cohort=job$cohort,sample=NA_character_,status='failed',seconds=NA_real_,error=conditionMessage(e)))
}
.da_run_fits <- function(root,state,jobs,workers,script) {
  cl<-parallel::makePSOCKcluster(min(workers,length(jobs)),rscript_args='--vanilla',outfile=file.path(root,'workers.log'))
  on.exit(parallel::stopCluster(cl),add=TRUE)
  parallel::clusterCall(cl,function(script,paths,state){
    .libPaths(paths)
    stopifnot(identical(unname(tools::md5sum(script)),state$script_md5))
    sys.source(script,.GlobalEnv)
    .da_assert(identical(.da_runtime(),state$runtime),'Worker R/package runtime differs')
    state$ref<-.da_load(file.path(state$root,'inputs','reference.rds'))
    assign('.da_worker_state',state,.GlobalEnv);TRUE
  },script,.libPaths(),state)
  ans<-parallel::parLapplyLB(cl,jobs,.da_worker)
  ans<-do.call(rbind,ans);.da_csv(ans,file.path(root,'latest_fit_status.csv'))
  .da_assert(!any(ans$status=='failed'),'Some CDSeq2 fits failed. See latest_fit_status.csv and workers.log; completed fits are preserved.')
  invisible(ans)
}
.da_assemble <- function(root,cohort,input,id) {
  a<-array(NA_real_,dim(input$truth),dimnames(input$truth));types<-dimnames(a)[[2]]
  theta<-matrix(NA_real_,dim(a)[3],length(types),dimnames=list(dimnames(a)[[3]],types))
  prior<-NULL
  for(i in seq_len(dim(a)[3])) {
    z<-.da_load(file.path(root,'sample_fits',paste0(cohort,'_',sprintf('%03d',i),'.rds')))
    .da_assert(identical(z$signature$run_id,id)&&identical(z$signature$sample,dimnames(a)[[3]][i])&&
      identical(z$signature$bulk_hash,.da_hash(input$bulk[,i,drop=FALSE])),'Wrong sample fit identity')
    .da_assert(identical(dimnames(z$result$counts),dimnames(a)[1:2]),'Fitted count axes differ')
    a[,,i]<-z$result$counts;theta[i,]<-z$result$proportions[types]
    if(is.null(prior))prior<-z$result$reference_mean else .da_assert(max(abs(prior-z$result$reference_mean))<1e-12,'Fitting prior differs by donor')
  }
  .da_assert(all(is.finite(a))&&all(a>=0)&&max(abs(apply(a,c(1,3),sum)-input$bulk)/pmax(input$bulk,1))<1e-6,
    'Assembled allocations do not conserve bulk')
  list(counts=a,proportions=theta,fitted_reference_prior_mean=prior)
}
.da_metrics <- function(d,retentions=seq(0,1,.05)) {
  summary<-curves<-list()
  rankers<-c(GAM7='variation_score',contribution='contribution_mean',expected_count='expected_log_mean_count')
  for(ct in unique(d$cell_type)) {
    z<-d[d$cell_type==ct,,drop=FALSE];ok<-is.finite(z$observed_r)&is.finite(z$predicted_r)
    summary[[ct]]<-data.frame(cell_type=ct,profiles=nrow(z),defined_truth=sum(is.finite(z$observed_r)),
      finite_scores=sum(is.finite(z$predicted_r)),evaluated=sum(ok),
      MAE=if(any(ok))mean(abs(z$observed_r[ok]-z$predicted_r[ok]))else NA_real_,
      RMSE=if(any(ok))sqrt(mean((z$observed_r[ok]-z$predicted_r[ok])^2))else NA_real_,
      rank_Spearman=.da_cor(z$observed_r[ok],z$predicted_r[ok],'spearman'),
      range_flagged=sum(z$outside_training_range%in%TRUE),mean_estimated_RNA_fraction=mean(z$theta_mean))
    # Eligibility uses observable ranking availability, not truth availability.
    common<-is.finite(z$variation_score)&is.finite(z$contribution_mean)&is.finite(z$expected_log_mean_count)
    for(method in names(rankers)) {
      ix<-which(common);ix<-ix[order(-z[[rankers[[method]]]][ix],z$gene[ix],method='radix')]
      for(frac in retentions) {
        jj<-head(ix,floor(length(ix)*frac+1e-9));r<-z$observed_r[jj];sp<-z$observed_spearman[jj]
        good<-is.finite(r);gs<-is.finite(sp)
        curves[[length(curves)+1L]]<-data.frame(cell_type=ct,method=method,retention=frac,
          all_profiles=nrow(z),rankable=sum(common),retained=length(jj),defined_Pearson=sum(good),
          mean_Pearson=if(any(good))mean(r[good])else NA_real_,
          median_Pearson=if(any(good))median(r[good])else NA_real_,
          fraction_Pearson_ge_0_5=if(any(good))mean(r[good]>=.5)else NA_real_,
          defined_Spearman=sum(gs),median_Spearman=if(any(gs))median(sp[gs])else NA_real_)
      }
    }
  }
  list(summary=do.call(rbind,summary),retention=do.call(rbind,curves))
}
.da_join_predictions <- function(d,p) {
  .da_assert(identical(d$gene,p$gene)&&identical(d$cell_type,p$cell_type),'Prediction row mismatch')
  cbind(d,p[setdiff(names(p),names(d))])
}
.da_train <- function(d,root,id,context,donors,refs,mu,folds,seed,min_profiles) {
  genes<-unique(d$gene);types<-unique(d$cell_type)
  fold_by_gene<-.da_seed(seed,function()setNames(sample(rep(seq_len(folds),length.out=length(genes))),sort(genes,method='radix')))
  d$fold<-unname(fold_by_gene[d$gene]);d$gene_group<-d$gene
  .da_save(fold_by_gene,file.path(root,'analysis','fold_by_gene.rds'))
  models<-list();oof<-list();audit<-list()
  for(ct in types)for(fold in c(seq_len(folds),NA_integer_)) {
    dd<-d[d$cell_type==ct,,drop=FALSE]
    train<-if(is.na(fold))dd else dd[dd$fold!=fold,,drop=FALSE]
    test<-if(is.na(fold))dd[FALSE,,drop=FALSE] else dd[dd$fold==fold,,drop=FALSE]
    .da_assert(!length(intersect(train$gene_group,test$gene_group)),'Grouped gene leakage')
    label<-if(is.na(fold))'full'else paste0('fold',fold)
    path<-file.path(root,'gam_fits',paste0(ct,'_',label,'.rds'))
    sig<-list(id=id,cell_type=ct,fold=fold,training=.da_hash(train),min_profiles=min_profiles)
    if(file.exists(path)) {z<-.da_load(path);.da_assert(identical(z$signature,sig),'GAM cache identity differs');m<-z$model}
    else {.da_log(root,'GAM ',ct,' ',label);m<-.da_fit_gam(train,min_profiles);.da_save(list(signature=sig,model=m),path)}
    audit[[length(audit)+1L]]<-data.frame(cell_type=ct,fold=label,training_rows=m$n_training,training_genes=m$n_distinct_genes)
    if(is.na(fold))models[[ct]]<-m else {
      p<-predict_direct_variation(test,list(features=.da_features,models=setNames(list(m),ct)))
      oof[[length(oof)+1L]]<-.da_join_predictions(test,p)
    }
  }
  oof<-do.call(rbind,oof);rownames(oof)<-NULL
  model<-list(schema='CDSeq2_variation_GAM7_v1',method='GAM7_no_theta',features=.da_features,
    models=models,gene_panel=genes,context=context,fit_settings=context$fit_settings,
    target='Across-donor normalized Pearson correlation',oracle='direct_observed_component_counts',
    donor_ids=donors,reference_donor_ids=refs,reference_hashes=character(),training_benchmarks=id,
    validation_status='grouped_gene_development_CV_only',model_id=paste0('direct_GAM7_',id),
    development=list(predictions=oof,metrics=.da_metrics(oof),fold_by_gene=fold_by_gene,fit_audit=do.call(rbind,audit)),
    provenance=list(experiment_version='1.0.0',feature_definition='variation_interface_2.3.0',
      created=as.character(Sys.time()),runtime=.da_runtime(),seed=seed,folds=folds))
  sent<-d[unlist(lapply(types,function(ct)head(which(d$cell_type==ct),10L))),,drop=FALSE]
  saved<-predict_direct_variation(sent,model)[c('gene','cell_type','predicted_r','variation_score','outside_training_range','range_features')]
  names(saved)[3:6]<-paste0('GAM7_no_theta_',c('predicted_r','score','outside_range','range_features'))
  model$sentinel<-list(features=sent,predictions=saved)
  class(model)<-c('cdseq2_variation_scorer','list')
  # Reference hash uses the package API serialization convention.
  tf<-tempfile();writeBin(serialize(mu,NULL,version=2),tf);model$reference_hashes<-.da_md5(tf);unlink(tf)
  .da_save(model,file.path(root,'analysis','variation_scorer.rds'))
  reloaded<-.da_load(file.path(root,'analysis','variation_scorer.rds'))
  .da_assert(identical(predict_direct_variation(sent,model),predict_direct_variation(sent,reloaded)), 'Saved scorer replay failed')
  if('check_variation_scorer'%in%getNamespaceExports('CDSeq2'))
    getExportedValue('CDSeq2','check_variation_scorer')(reloaded)
  model
}

.da_plots <- function(d,metrics,path) {
  grDevices::pdf(path,width=11,height=8.5)
  old<-par(mfrow=c(2,4),mar=c(4,4,3,1));on.exit({par(old);grDevices::dev.off()},add=TRUE)
  colors<-c(GAM7='#0072B2',contribution='#D55E00',expected_count='#009E73')
  for(ct in unique(d$cell_type)) {
    z<-metrics$retention[metrics$retention$cell_type==ct,,drop=FALSE]
    plot(NA,xlim=c(0,1),ylim=c(-1,1),xlab='Fraction retained',ylab='Mean normalized Pearson',main=ct)
    for(m in names(colors)) {q<-z[z$method==m,];lines(q$retention,q$mean_Pearson,col=colors[m],lwd=2)}
    abline(h=0,col='grey80',lty=3)
  }
  plot.new();legend('center',legend=names(colors),col=colors,lwd=2,bty='n')
  for(ct in unique(d$cell_type)) {
    z<-d[d$cell_type==ct,];plot(z$predicted_r,z$observed_r,pch=16,cex=.35,
      col=grDevices::adjustcolor('#0072B2',.35),xlim=c(-1,1),ylim=c(-1,1),
      xlab='Predicted Pearson',ylab='Observed Pearson',main=ct);abline(0,1,col='grey40',lty=2)
  }
  plot.new();text(.5,.5,'Direct aggregation\nValidation donors\nNo model tuning')
}
.da_analysis <- function(root,input,id,folds,seed,min_profiles,runtime) {
  dat<-list();estimates<-list()
  for(cohort in c('training','validation')) {
    est<-.da_assemble(root,cohort,input[[cohort]],id);estimates[[cohort]]<-est
    .da_save(est,file.path(root,'analysis',paste0(cohort,'_estimates.rds')))
    f<-.da_observables(input[[cohort]]$bulk,input$reference_mean,est$proportions)
    targets<-.da_targets(input[[cohort]]$truth,est$counts)
    .da_assert(identical(f$gene,targets$gene)&&identical(f$cell_type,targets$cell_type),'Feature/target keys differ')
    d<-cbind(f,targets[setdiff(names(targets),names(f))]);d$cohort_id<-paste0('direct_',cohort)
    dat[[cohort]]<-d;.da_save(d,file.path(root,'analysis',paste0(cohort,'_features_targets.rds')))
    .da_csv(d,file.path(root,'analysis',paste0(cohort,'_features_targets.csv')))
    truth_theta<-t(sweep(apply(input[[cohort]]$truth,c(2,3),sum),2,colSums(input[[cohort]]$bulk),'/'))
    prop<-do.call(rbind,lapply(colnames(est$proportions),function(ct)data.frame(cell_type=ct,
      mean_true_RNA_fraction=mean(truth_theta[,ct]),mean_estimated_RNA_fraction=mean(est$proportions[,ct]),
      MAE=mean(abs(truth_theta[,ct]-est$proportions[,ct])),
      Pearson=.da_cor(truth_theta[,ct],est$proportions[,ct]))))
    .da_csv(prop,file.path(root,'analysis',paste0(cohort,'_RNA_fraction_accuracy.csv')))
  }
  context<-list(application='ROSMAP direct aggregation',tissue='PFC',deconvolution_method='CDSeq2',
    fit_settings=input$fit,simulation='none; sum observed donor/type components',
    oracle='direct_observed_component_counts',range_policy='flag',
    validation_history='Donor-disjoint; these validation donors were examined in the previous resampled experiment')
  # Training function receives training rows only. Validation truth never enters fitting.
  model<-.da_train(dat$training,root,id,context,input$training$donors,input$reference_donors,
    input$reference_mean,folds,seed,min_profiles)
  validation<-.da_join_predictions(dat$validation,predict_direct_variation(dat$validation,model))
  metrics<-.da_metrics(validation)
  .da_save(list(model_id=model$model_id,predictions=validation,metrics=metrics,
    status='completed_donor_disjoint_direct_aggregation_evaluation',context=context),file.path(root,'analysis','validation_results.rds'))
  .da_csv(validation,file.path(root,'analysis','validation_predictions.csv'))
  .da_csv(metrics$summary,file.path(root,'analysis','validation_summary_by_type.csv'))
  .da_csv(metrics$retention,file.path(root,'analysis','validation_retention_curves.csv'))
  .da_csv(model$development$predictions,file.path(root,'analysis','training_OOF_predictions.csv'))
  .da_csv(model$development$metrics$summary,file.path(root,'analysis','training_OOF_summary.csv'))
  .da_csv(model$development$metrics$retention,file.path(root,'analysis','training_OOF_retention_curves.csv'))
  .da_csv(model$development$fit_audit,file.path(root,'analysis','GAM_fit_audit.csv'))
  selected<-validation;selected$retained_top20<-FALSE;selected$rank_within_type<-NA_integer_
  for(ct in unique(selected$cell_type)) {
    ii<-which(selected$cell_type==ct&is.finite(selected$variation_score))
    ii<-ii[order(-selected$variation_score[ii],selected$gene[ii],method='radix')]
    selected$rank_within_type[ii]<-seq_along(ii)
    selected$retained_top20[head(ii,floor(.2*length(ii)+1e-9))]<-TRUE
  }
  .da_csv(selected,file.path(root,'analysis','validation_top20_profiles.csv'))
  .da_plots(validation,metrics,file.path(root,'analysis','validation_figures.pdf'))
  writeLines(c('COMPLETED: direct aggregation profile recovery experiment.',
    '100 distinct donors: 20 reference, 40 training, 40 validation.',
    '80 direct-bulk CDSeq2 fits; 35 grouped-gene CV GAM fits plus 7 full GAM fits.',
    'Truth = observed components, normalized over the complete 2000-gene panel within donor/type.',
    'No multinomial resampling, perturbation, phenotype generation or DE testing.',
    'Estimated RNA fractions enter scorer features only, not pseudobulk construction.',
    'Development CV holds out genes, not donors. Final validation holds out biological donors.',
    'Validation donors were previously examined under a different count-generation experiment.',
    'Results assess recovery of observed single-nucleus components; no claim of noise-free biology or real-bulk transportability.',
    'Report each cell type, RNA support, undefined targets and extrapolation flags.',
    'Scores rank profile recovery; they are not association probabilities, p-values or FDR guarantees.',
    'The top20 cutoff is descriptive and was not optimized on validation truth.'),file.path(root,'SUMMARY.txt'))
  writeLines(as.character(Sys.time()),file.path(root,'COMPLETED'))
  .da_log(root,'COMPLETED. Results: ',file.path(root,'analysis'))
  invisible(list(run_dir=root,model_file=file.path(root,'analysis','variation_scorer.rds')))
}

#' Run the verified ROSMAP direct aggregation experiment
#' @param action prepare: audit and save inputs only; probe: fit largest training
#'   donor only; run: all remaining fits and analysis; analyze: use complete fits.
#' @param workers Fresh PSOCK fitting workers, each with cpu_number=1. Default 4.
#'   Reference is ~0.64 GB as a dense matrix per worker before fitting overhead.
#' @param training_file Original registered_base benchmark, NOT the sim_final file.
#' @param reference_file Exact verified 39,692-nucleus reference SCE.
#' @param transfer_dir Folder containing inputs/ and FRESH_REQUEST/ from transfer ZIP.
#' @param output_dir Dedicated results directory. Changed code/data/settings/runtime
#'   require a new directory. Existing incompatible results are never silently reused.
#' @param folds Grouped gene development folds; default 5, same mapping in every type.
#' @param seed Fixed gene-fold seed; caller RNG is restored.
#' @param min_profiles Minimum distinct finite-target genes per GAM fit, default 100.
#' @return Invisibly, run information. Source alone never launches this function.
run_direct_aggregation <- function(action=c('prepare','probe','run','analyze'),workers=4L,
    training_file='~/Downloads/score_validation_v1/registered_base/rosmap_heldout_benchmark_minimal.rds',
    reference_file='~/Downloads/score_validation_v1/score_validation_data/reference/rosmap_reference_sce.rds',
    transfer_dir=NULL,output_dir='~/Downloads/CDSeq2_direct_aggregation_v1',
    folds=5L,seed=20261009L,min_profiles=100L) {
  action<-match.arg(action)
  for(nm in c('workers','folds','seed','min_profiles')) {
    v<-get(nm);.da_assert(is.numeric(v)&&length(v)==1L&&is.finite(v)&&v==floor(v)&&v>=1&&v<=.Machine$integer.max,'Invalid ',nm)
  }
  .da_assert(folds>=2L&&folds<=20L&&min_profiles>=50L,'Require 2--20 folds and min_profiles >=50')
  .da_assert(!is.na(.da_script)&&file.exists(.da_script),'Source the saved script before calling this function')
  if('CDSeq2'%in%loadedNamespaces()) {
    ns<-asNamespace('CDSeq2')
    .da_assert(!isTRUE(get0('.__DEVTOOLS__',envir=ns,inherits=FALSE)),
      'Restart R and use the installed CDSeq2 build, not load_all()')
  }
  runtime<-.da_runtime()
  training_file<-.da_path(training_file);reference_file<-.da_path(reference_file);transfer_dir<-.da_find_transfer(transfer_dir)
  root<-path.expand(output_dir);dir.create(root,recursive=TRUE,showWarnings=FALSE);root<-.da_path(root)
  lock<-file.path(root,'RUNNING.lock')
  .da_assert(dir.create(lock,showWarnings=FALSE),'Run lock exists: ',lock,
    '. If a previous process crashed, verify no controller/workers are running before manually removing this directory.')
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  writeLines(c(paste('PID',Sys.getpid()),paste('HOST',Sys.info()[['nodename']])),file.path(lock,'owner.txt'))
  for(d in c('inputs','sample_fits','gam_fits','analysis'))dir.create(file.path(root,d),showWarnings=FALSE)
  .da_log(root,'Checking exact inputs, original components, reference and biological donor separation')
  input<-.da_inputs(training_file,reference_file,transfer_dir)
  required<-c('bulk_data','reference_scRNA_count','reference_scRNA_annotation','cpu_number',names(input$fit))
  .da_assert(all(required%in%names(formals(getExportedValue('CDSeq2','cdseq2')))),
    'Installed CDSeq2 interface does not support the verified fitting settings')
  config<-list(version='direct_aggregation_1.0.0',script_md5=.da_md5(.da_script),inputs=input$checksums,
    runtime=runtime,fit=input$fit,folds=as.integer(folds),seed=as.integer(seed),min_profiles=as.integer(min_profiles),
    fit_seed_base=20301009L,oracle='direct_observed_component_counts')
  id<-.da_hash(config);rp<-file.path(root,'experiment.rds')
  if(file.exists(rp)) .da_assert(identical(.da_load(rp)$config,config),
    'Existing run differs in code, inputs, settings or runtime. Choose a new output_dir.')
  else {
    existing<-list.files(file.path(root,'sample_fits'),all.files=FALSE)
    .da_assert(!length(existing),'Sample fits exist without experiment provenance; choose a new output_dir')
    .da_save(list(id=id,config=config),rp)
  }
  copied<-file.path(root,'workflow.R')
  if(!identical(.da_path(.da_script),normalizePath(copied,mustWork=FALSE)))
    .da_assert(file.copy(.da_script,copied,overwrite=TRUE),'Cannot preserve script')
  .da_save(input$reference,file.path(root,'inputs','reference.rds'));input$reference<-NULL
  .da_save(input,file.path(root,'inputs','direct_inputs.rds'))
  donors<-rbind(data.frame(role='reference',donor=input$reference_donors),
    data.frame(role='training',donor=input$training$donors),data.frame(role='validation',donor=input$validation$donors))
  .da_csv(donors,file.path(root,'donor_split_PRIVATE.csv'))
  .da_csv(data.frame(check=c('unique_reference','unique_training','unique_validation','unique_total',
    'training_validation_overlap','training_reference_overlap','validation_reference_overlap'),value=c(20,40,40,100,0,0,0)),file.path(root,'donor_checks.csv'))
  coverage<-do.call(rbind,lapply(c('training','validation'),function(cohort){
    a<-input[[cohort]]$truth;lib<-apply(a,c(2,3),sum)
    data.frame(cohort=cohort,cell_type=rownames(lib),min_component_library=apply(lib,1,min),
      median_component_library=apply(lib,1,median),mean_true_RNA_fraction=rowMeans(sweep(lib,2,colSums(lib),'/')))
  }))
  .da_csv(coverage,file.path(root,'input_RNA_support.csv'))
  capture.output(utils::sessionInfo(),file=file.path(root,'sessionInfo.txt'))
  .da_log(root,'VERIFIED: 100 biological donors; direct observed counts; run ID ',id)
  if(action=='prepare')return(invisible(list(run_dir=root,id=id,status='prepared_no_fits')))
  jobs<-list();number<-0L
  for(cohort in c('training','validation'))for(i in seq_len(ncol(input[[cohort]]$bulk))) {
    number<-number+1L;jobs[[number]]<-list(cohort=cohort,i=i,seed=as.integer(config$fit_seed_base+number))
  }
  state<-list(root=root,id=id,fit=input$fit,runtime=runtime,script_md5=config$script_md5,types=colnames(input$reference_mean),
    bulks=lapply(input[c('training','validation')],`[[`,'bulk'))
  if(action!='analyze') {
    # Probe largest training library first; its fit is reused in the full run.
    probe<-which.max(colSums(input$training$bulk))
    .da_log(root,'Checking/fitting one largest-library training donor first')
    .da_run_fits(root,state,jobs[probe],1L,copied)
    if(action=='probe')return(invisible(list(run_dir=root,id=id,status='probe_completed')))
    .da_log(root,'Fitting/resuming 80 direct-bulk samples using ',workers,' workers')
    .da_run_fits(root,state,jobs,as.integer(workers),copied)
  }
  .da_analysis(root,input,id,as.integer(folds),as.integer(seed),as.integer(min_profiles),runtime)
}
