# New diagnostic analyses; primary inference, selections and tests are unchanged.
# Rscript prioritization_audit.R INPUT_RDS LD_RESULTS OUTPUT_DIR
a<-commandArgs(TRUE); input<-readRDS(a[1]); ld<-a[2]; out<-a[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
tests<-read.csv(file.path(ld,'test_results.csv'));expr<-read.csv(file.path(ld,'expression_recovery.csv'))
features<-input$ld_features; key<-function(d)paste(d$id,d$gene,d$cell_type,sep='|')
d<-tests[tests$method=='CDSeq2',];f<-features[match(key(d),key(features)),];e<-expr[expr$method=='CDSeq2',];e<-e[match(key(d),key(e)),]
stopifnot(!anyNA(f$gene),!anyNA(e$gene),identical(is.na(d$score),is.na(f$variation_score)),max(abs(d$score-f$variation_score),na.rm=TRUE)<1e-10)
d$contribution<-f$contribution_mean;d$expected_count<-f$expected_log_mean_count;d$recovery<-e$r_CPM
ranker<-function(x,g){i<-which(is.finite(x));r<-rep(NA_integer_,length(x));r[i[order(-x[i],g[i],method='radix')]]<-seq_along(i);r}
metrics<-function(z,keep){called<-keep & z$available & z$q<=.05;tp<-sum(called&z$injected);fp<-sum(called&!z$injected);n<-sum(z$injected);r<-z$recovery[keep&is.finite(z$recovery)];data.frame(selected=sum(keep),defined_recovery=length(r),median_recovery=if(length(r))median(r)else NA_real_,calls=tp+fp,TP=tp,FP=fp,targets=n,precision=if(tp+fp)tp/(tp+fp)else NA_real_,recall=if(n)tp/n else NA_real_,F1=if(n)2*tp/(tp+fp+n)else NA_real_)}
dec<-controls<-list();set.seed(20261020)
groups<-split(d,paste(d$id,d$cell_type));
for(z in groups){
 r<-ranker(z$score,z$gene); n<-sum(is.finite(r)); decile<-ceiling(r/n*10)
 for(b in 1:10){keep<-!is.na(decile)&decile==b;dec[[length(dec)+1]]<-cbind(z[1,c('id','setting','cell_type')],decile=b,metrics(z,keep))}
 for(method in c('GAM6','Contribution','Expected_count','Random')){
  nr<-if(method=='Random')100 else 1
  for(j in seq_len(nr)){
   x<-switch(method,GAM6=z$score,Contribution=z$contribution,Expected_count=z$expected_count,Random=runif(nrow(z)))
   rr<-ranker(x,z$gene);keep<-is.finite(rr)&rr<=floor(.4*nrow(z))
   stopifnot(sum(keep)==80)
   controls[[length(controls)+1]]<-cbind(z[1,c('id','setting','cell_type')],method=method,random_repeat=if(method=='Random')j else 0,metrics(z,keep))
  }
 }
}
dec<-do.call(rbind,dec);controls<-do.call(rbind,controls)
write.csv(dec,file.path(out,'S36_score_deciles.csv'),row.names=FALSE,na='')
write.csv(controls,file.path(out,'S37_matched_retention_controls.csv'),row.names=FALSE,na='')
# Resample whole cohorts, preserving paired comparisons and repeated genes.
# Random selections are Monte Carlo-averaged within each cohort, not independent cohorts.
cm<-aggregate(controls[c('TP','FP','targets','selected','median_recovery')],controls[c('id','setting','cell_type','method')],mean)
stat<-function(x) c(precision=if(sum(x$TP+x$FP))sum(x$TP)/sum(x$TP+x$FP)else NA_real_,recall=if(sum(x$targets))sum(x$TP)/sum(x$targets)else NA_real_,median_recovery=mean(x$median_recovery))
ci<-list();set.seed(20261021)
for(z in split(cm,paste(cm$setting,cm$cell_type))){
 ids<-sort(unique(z$id));draw<-replicate(2000,sample(ids,length(ids),TRUE));base<-z[z$method=='GAM6',];base<-base[match(ids,base$id),]
 for(method in unique(z$method)){
  x<-z[z$method==method,];x<-x[match(ids,x$id),];boot<-apply(draw,2,function(ii)stat(x[match(ii,ids),]));diff<-apply(draw,2,function(ii)stat(base[match(ii,ids),])-stat(x[match(ii,ids),]))
  for(metric in rownames(boot)){
   q<-function(v) if(any(is.finite(v)))quantile(v[is.finite(v)],c(.025,.975),names=FALSE)else c(NA,NA)
   b<-q(boot[metric,]);dif<-q(diff[metric,]);ci[[length(ci)+1]]<-data.frame(setting=z$setting[1],cell_type=z$cell_type[1],method=method,metric=metric,estimate=stat(x)[metric],lower=b[1],upper=b[2],GAM6_minus_method=stat(base)[metric]-stat(x)[metric],difference_lower=dif[1],difference_upper=dif[2],cohorts=length(ids),bootstrap_repetitions=2000)
  }
 }
}
write.csv(do.call(rbind,ci),file.path(out,'S38_control_uncertainty.csv'),row.names=FALSE,na='')
cat('Score groups and equal-budget controls complete.\n')
