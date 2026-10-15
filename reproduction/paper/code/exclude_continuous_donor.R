# Reanalysis of archived continuous simulations after a post hoc donor exclusion.
.ca_assert <- function(x,...) if(!isTRUE(x)) stop(...,call.=FALSE)
.ca_md5 <- function(p) unname(tools::md5sum(p))
.ca_read <- function(p) {
  .ca_assert(file.exists(p)&&file.exists(paste0(p,'.md5')), 'Missing verified input: ',p)
  .ca_assert(identical(.ca_md5(p),readLines(paste0(p,'.md5'),warn=FALSE)), 'Input checksum failed: ',p)
  readRDS(p)
}
.ca_save <- function(x,p) {
  tmp<-tempfile(tmpdir=dirname(p));on.exit(unlink(tmp));saveRDS(x,tmp,version=2)
  .ca_assert(file.rename(tmp,p),'Cannot save ',p)
  writeLines(.ca_md5(p),paste0(p,'.md5'))
}
.ca_key <- function(d) paste(d$gene,d$cell_type,sep='\r')
.ca_cor <- function(x,y) {
  ok<-is.finite(x)&is.finite(y)
  if(sum(ok)<3L||stats::sd(x[ok])<1e-12||stats::sd(y[ok])<1e-12) return(NA_real_)
  unname(stats::cor(x[ok],y[ok]))
}
.ca_mean <- function(x) if(any(is.finite(x)))mean(x[is.finite(x)])else NA_real_
.ca_div <- function(a,b) if(b>0)a/b else NA_real_
.ca_expression <- function(a) {
  .ca_assert(is.array(a)&&length(dim(a))==3L&&all(is.finite(a))&&all(a>=0), 'Invalid component counts')
  g<-dimnames(a)[[1]];ct<-dimnames(a)[[2]];s<-dimnames(a)[[3]]
  for(v in list(g,ct,s)) .ca_assert(length(v)>0&&!anyNA(v)&&!anyDuplicated(v)&&all(nzchar(v)), 'Invalid count axes')
  z<-do.call(rbind,lapply(seq_along(ct),function(k){
    m<-matrix(a[,k,],nrow=length(g));lib<-colSums(m)
    if(any(lib<=0)) return(matrix(NA_real_,length(g),length(s)))
    log2(sweep(m,2,lib,'/')*1e6+.5)
  }))
  colnames(z)<-s
  list(x=z,keys=data.frame(gene=rep(g,length(ct)),cell_type=rep(ct,each=length(g)),stringsAsFactors=FALSE))
}
.ca_residual <- function(x,q) {
  ok<-rowSums(!is.finite(x))==0L
  z<-matrix(NA_real_,nrow(x),ncol(x),dimnames=dimnames(x))
  if(any(ok))z[ok,]<-t(qr.resid(q,t(x[ok,,drop=FALSE])))
  z
}
.ca_associate <- function(x,y,df) {
  # x and y have already been residualized against identical covariates.
  ss<-rowSums(x*x); sy<-sum(y*y)
  r<-drop(x%*%y)/sqrt(ss*sy)
  r[!is.finite(r)|ss<1e-20|sy<1e-20]<-NA_real_
  r<-pmax(-1,pmin(1,r));stat<-abs(r)*sqrt(df/pmax(1-r*r,1e-15))
  p<-2*stats::pt(stat,df=df,lower.tail=FALSE)
  # Fixed full panel x type family; untestable rows count conservatively as p=1.
  q<-stats::p.adjust(replace(p,!is.finite(p),1),method='BH')
  q[!is.finite(p)]<-NA_real_
  list(r=r,p=p,q=q)
}


args<-commandArgs(TRUE)
if(length(args)<2L)stop('Usage: Rscript exclude_continuous_donor.R EXPERIMENT_ROOT OUTPUT_FOLDER [DETAIL_CACHE]')
p<-normalizePath(args[1]);o<-args[2];td<-file.path(o,'supplementary_tables')
cache<-if(length(args)>=3L)args[3]else NULL
if(!is.null(cache)){dir.create(cache,recursive=TRUE,showWarnings=FALSE);dir.create(file.path(cache,'signal_tests'),showWarnings=FALSE)}
sources<-list()
load<-function(rel){f<-file.path(p,rel);x<-.ca_read(f);sources[[length(sources)+1L]]<<-data.frame(source=rel,md5=.ca_md5(f));x}
input<-load('inputs/direct_inputs.rds');estimate<-load('analysis/validation_estimates.rds')
phenotypes<-load('continuous_associations_v1/synthetic_phenotypes.rds')
null_y<-load('continuous_associations_v1/null_phenotypes.rds');sets<-load('continuous_associations_v1/selection_sets.rds')
stopifnot(identical(dimnames(input$validation$truth),dimnames(estimate$counts)))
tt<-.ca_expression(input$validation$truth);ee<-.ca_expression(estimate$counts);keys<-tt$keys
samples<-colnames(tt$x);stopifnot(length(samples)==40,identical(rownames(null_y),samples),ncol(null_y)==200)
# Single fixed exclusion from every continuous-phenotype signal/null analysis.
# Existing simulated values and ranks are retained; no new PCs/noise/inference.
excluded<-3L;keep<-setdiff(seq_along(samples),excluded);n<-length(keep);df<-n-2L
stopifnot(which.min(scale(phenotypes[['trait_031']]$phenotype))==excluded)
# Verify alignment and the association implementation on three archived 40-donor traits.
cq0<-qr(matrix(1,40,1));tx0<-.ca_residual(tt$x,cq0);ex0<-.ca_residual(ee$x,cq0)
for(id in c('trait_001','trait_031','trait_189')){
 y<-phenotypes[[id]]$phenotype;stopifnot(identical(names(y),samples))
 a<-.ca_associate(tx0,drop(qr.resid(cq0,y)),38);b<-.ca_associate(ex0,drop(qr.resid(cq0,y)),38)
 old<-load(file.path('continuous_associations_v1/signal_tests',paste0(id,'.rds')))
 stopifnot(identical(.ca_key(old),.ca_key(keys)))
 for(nm in c('r','p','q'))for(en in c('reference','estimated')){
  v<-if(en=='reference')a[[nm]]else b[[nm]];r<-old[[paste0(en,'_',nm)]]
  stopifnot(identical(is.na(v),is.na(r)),max(abs(v-r),na.rm=TRUE)<1e-10)
 }
}
# Frozen membership uses observable archived scores; no re-ranking on outcomes.
stopifnot(length(sets)==280)
full_ix<-split(seq_len(nrow(keys)),keys$cell_type)
cq<-qr(matrix(1,n,1));tx<-.ca_residual(tt$x[,keep],cq);ex<-.ca_residual(ee$x[,keep],cq)
# Per-trait statistics use 37 residual degrees of freedom and BH over all 14,000 tests.
records<-vector('list',length(phenotypes))
for(i in seq_along(phenotypes)){
 t<-phenotypes[[i]];stopifnot(identical(names(t$phenotype),samples))
 yr<-drop(qr.resid(cq,t$phenotype[keep]));a<-.ca_associate(tx,yr,df);b<-.ca_associate(ex,yr,df)
 ac_all<-is.finite(a$q)&a$q<=.05;bc_all<-is.finite(b$q)&b$q<=.05
 agree_all<-ac_all&bc_all&sign(a$r)==sign(b$r);agree_all[is.na(agree_all)]<-FALSE
 total_ref<-vapply(full_ix,function(j)sum(ac_all[j]),integer(1))
 rows<-lapply(sets,function(z){
  j<-z$ix;defined<-is.finite(a$r[j])&is.finite(b$r[j]);ac<-ac_all[j];bc<-bc_all[j];agree<-agree_all[j]
  data.frame(trait=t$id,source_type=t$source_type,program=t$program,signal_fraction=t$signal_fraction,replicate=t$replicate,
   evaluated_type=z$cell_type,method=z$method,retention=z$retention,rankable=z$rankable,selected=length(j),paired_defined=sum(defined),
   association_r_MAE=.ca_mean(abs(a$r[j]-b$r[j])),association_r_concordance=.ca_cor(a$r[j],b$r[j]),
   reference_calls=sum(ac),estimated_calls=sum(bc),same_sign_overlap=sum(agree),
   reference_agreement_among_calls=.ca_div(sum(agree),sum(bc)),reference_recall_within_selected=.ca_div(sum(agree),sum(ac)),
   reference_recall_full_type=.ca_div(sum(agree),total_ref[[z$cell_type]]),reference_unavailable_among_calls=sum(bc&!is.finite(a$q[j])))
 })
 records[[i]]<-do.call(rbind,rows)
 if(!is.null(cache)).ca_save(cbind(keys,data.frame(reference_r=a$r,reference_p=a$p,reference_q=a$q,estimated_r=b$r,estimated_p=b$p,estimated_q=b$q)),file.path(cache,'signal_tests',paste0(t$id,'.rds')))
 if(i%%27==0)cat('Recomputed',i,'of',length(phenotypes),'signal phenotypes\n')
}
result<-do.call(rbind,records)
metrics<-c('association_r_MAE','association_r_concordance','reference_agreement_among_calls','reference_recall_within_selected','reference_recall_full_type','estimated_calls','same_sign_overlap')
summary<-aggregate(result[metrics],result[c('source_type','signal_fraction','evaluated_type','method','retention')],.ca_mean)
result$n_donors<-summary$n_donors<-n;result$excluded_donor_index<-summary$excluded_donor_index<-excluded
write.csv(result,file.path(td,'S9_continuous_replicates.csv'),row.names=FALSE,na='')
write.csv(summary,file.path(td,'S3_continuous_summary.csv'),row.names=FALSE,na='')
nulls<-vector('list',ncol(null_y))
for(bi in seq_len(ncol(null_y))){
 yr<-drop(qr.resid(cq,null_y[keep,bi]));aa<-.ca_associate(tx,yr,df);bb<-.ca_associate(ex,yr,df);rows<-list()
 for(engine in c('observed_components','CDSeq2')){
  z<-if(engine=='CDSeq2')bb else aa;call<-is.finite(z$q)&z$q<=.05
  rows[[length(rows)+1L]]<-data.frame(replicate=bi,engine=engine,cell_type='ALL',method='unfiltered',retention=1,selected=length(call),false_calls=sum(call),any_false_call=any(call))
  for(ss in sets)rows[[length(rows)+1L]]<-data.frame(replicate=bi,engine=engine,cell_type=ss$cell_type,method=ss$method,retention=ss$retention,selected=length(ss$ix),false_calls=sum(call[ss$ix]),any_false_call=any(call[ss$ix]))
 }
 nulls[[bi]]<-do.call(rbind,rows)
 if(bi%%50==0)cat('Recomputed',bi,'of',ncol(null_y),'null phenotypes\n')
}
nulls<-do.call(rbind,nulls)
ns<-aggregate(nulls[c('false_calls','any_false_call')],nulls[c('engine','cell_type','method','retention')],mean)
names(ns)[names(ns)=='false_calls']<-'mean_false_calls';names(ns)[names(ns)=='any_false_call']<-'probability_any_false_call'
ns$null_replicates<-ncol(null_y);ns$Monte_Carlo_SE<-sqrt(ns$probability_any_false_call*(1-ns$probability_any_false_call)/ncol(null_y));ns$n_donors<-n;ns$excluded_donor_index<-excluded
write.csv(ns,file.path(td,'S4_null_results.csv'),row.names=FALSE,na='')
if(!is.null(cache)){write.csv(nulls,file.path(cache,'null_replicate_results.csv'),row.names=FALSE);saveRDS(list(donor_indices=keep,expression_reference=tt$x[,keep],expression_estimated=ee$x[,keep],keys=keys,phenotypes=lapply(phenotypes,function(t)t$phenotype[keep]),null_phenotypes=null_y[keep,],selection_sets=sets),file.path(cache,'verification_inputs.rds'))}
write.csv(do.call(rbind,sources),file.path(td,'continuous_exclusion_source_manifest.csv'),row.names=FALSE)
cat('Continuous analyses use',n,'donors; excluded original donor index',excluded,'\n')
print(subset(summary,source_type=='Excit'&evaluated_type=='Excit'&method=='GAM6'&retention %in% c(.2,1)),row.names=FALSE)
print(subset(ns,cell_type=='ALL'),row.names=FALSE)
