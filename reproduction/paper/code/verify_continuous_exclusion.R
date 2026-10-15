# Independent regression check of the donor-excluded continuous analyses.
# DETAIL_CACHE is produced by exclude_continuous_donor.R and is not distributed
# because it contains donor-level expression and original sample identifiers.
a <- commandArgs(TRUE)
if(length(a)!=2L)stop('Usage: Rscript verify_continuous_exclusion.R DETAIL_CACHE OUTPUT_FOLDER')
cache<-a[1];td<-file.path(a[2],'supplementary_tables')
v<-readRDS(file.path(cache,'verification_inputs.rds'))
stopifnot(identical(v$donor_indices,setdiff(1:40,3L)))
read<-function(n)read.csv(file.path(td,n),stringsAsFactors=FALSE)
ols<-function(x,y){
  design<-cbind(1,y);ok<-rowSums(!is.finite(x))==0
  ok[ok]<-apply(x[ok,,drop=FALSE],1,var)*(length(y)-1)>=1e-20
  fit<-.lm.fit(design,t(x[ok,,drop=FALSE]));s2<-colSums(fit$residuals^2)/(length(y)-2)
  se<-sqrt(s2*solve(crossprod(design))[2,2]);b<-fit$coefficients[2,]
  p<-slope<-rep(NA_real_,nrow(x));p[ok]<-2*pt(abs(b/se),df=length(y)-2,lower.tail=FALSE);slope[ok]<-b
  q<-p.adjust(replace(p,is.na(p),1),'BH');q[is.na(p)]<-NA_real_
  list(p=p,q=q,slope=slope,call=is.finite(q)&q<=.05)
}
close<-function(x,y,tol=1e-9){stopifnot(identical(is.na(x),is.na(y)));stopifnot(max(abs(x-y),na.rm=TRUE)<tol)}
replicates<-read('S9_continuous_replicates.csv')
for(id in names(v$phenotypes)){
  y<-v$phenotypes[[id]];ref<-ols(v$expression_reference,y);est<-ols(v$expression_estimated,y)
  d<-readRDS(file.path(cache,'signal_tests',paste0(id,'.rds')))
  stopifnot(identical(d$gene,v$keys$gene),identical(d$cell_type,v$keys$cell_type))
  close(ref$p,d$reference_p);close(ref$q,d$reference_q);close(est$p,d$estimated_p);close(est$q,d$estimated_q)
  z<-replicates[replicates$trait==id,];stopifnot(nrow(z)==length(v$selection_sets))
  overlap<-ref$call&est$call&sign(ref$slope)==sign(est$slope);overlap[is.na(overlap)]<-FALSE
  for(i in seq_along(v$selection_sets)){
    s<-v$selection_sets[[i]];j<-s$ix;r<-z[i,]
    stopifnot(r$method==s$method,r$evaluated_type==s$cell_type,abs(r$retention-s$retention)<1e-12,
      r$selected==length(j),r$reference_calls==sum(ref$call[j]),r$estimated_calls==sum(est$call[j]),r$same_sign_overlap==sum(overlap[j]))
  }
}
cat('PASS: all 189 signal traits; independent least-squares p/q and all 52,920 selection-level call counts\n')
nulls<-read.csv(file.path(cache,'null_replicate_results.csv'),stringsAsFactors=FALSE)
for(i in 1:ncol(v$null_phenotypes))for(engine in c('observed_components','CDSeq2')){
  x<-if(engine=='CDSeq2')v$expression_estimated else v$expression_reference
  r<-ols(x,v$null_phenotypes[,i]);z<-nulls[nulls$replicate==i&nulls$engine==engine,]
  expected<-c(sum(r$call),vapply(v$selection_sets,function(s)sum(r$call[s$ix]),integer(1)))
  stopifnot(identical(as.integer(z$false_calls),as.integer(expected)),all(z$any_false_call==(expected>0)))
}
cat('PASS: all 200 null traits, both expression matrices, and all selection-level rejection counts\n')
mean_defined<-function(x)if(any(is.finite(x)))mean(x[is.finite(x)])else NA_real_
metrics<-c('association_r_MAE','association_r_concordance','reference_agreement_among_calls','reference_recall_within_selected','reference_recall_full_type','estimated_calls','same_sign_overlap')
group<-c('source_type','signal_fraction','evaluated_type','method','retention')
s<-aggregate(replicates[metrics],replicates[group],mean_defined);reported<-read('S3_continuous_summary.csv')
stopifnot(identical(s[group],reported[group]));for(nm in metrics)close(s[[nm]],reported[[nm]])
n<-aggregate(nulls[c('false_calls','any_false_call')],nulls[c('engine','cell_type','method','retention')],mean)
reported<-read('S4_null_results.csv');stopifnot(identical(n[1:4],reported[1:4]))
close(n$false_calls,reported$mean_false_calls);close(n$any_false_call,reported$probability_any_false_call)
for(f in c('S3_continuous_summary.csv','S4_null_results.csv','S9_continuous_replicates.csv')){
 d<-read(f);stopifnot(all(d$n_donors==39),all(d$excluded_donor_index==3))
}
cat('PASS: aggregated signal/null tables and 39-donor exclusion metadata\n')
