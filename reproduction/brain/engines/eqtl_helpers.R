# Shared numerical/reporting helpers, v1.0.0. Sourcing never starts analysis.
.eqs_assert <- function(ok, ...) if (!isTRUE(ok)) stop(..., call.=FALSE)
.eqs_key <- function(d) paste(d$gene,d$cell_type,sep='\r')
.eqs_md5 <- function(p) unname(tools::md5sum(p))
.eqs_hash <- function(x) {
 p<-tempfile();on.exit(unlink(p));saveRDS(x,p,version=2);.eqs_md5(p)
}
.eqs_save <- function(x,p) {
 dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
 tmp<-tempfile(tmpdir=dirname(p));on.exit(unlink(tmp));saveRDS(x,tmp,version=2)
 .eqs_assert(file.rename(tmp,p),'Cannot save ',p)
 writeLines(.eqs_md5(p),paste0(p,'.md5'));invisible(p)
}
.eqs_load <- function(p) {
 .eqs_assert(file.exists(p)&&file.exists(paste0(p,'.md5'))&&
  identical(.eqs_md5(p),readLines(paste0(p,'.md5'),warn=FALSE)),'Missing/changed checkpoint: ',p)
 readRDS(p)
}
.eqs_csv <- function(d,p) utils::write.csv(d,p,row.names=FALSE,na='')
.eqs_div <- function(a,b) ifelse(b>0,a/b,NA_real_)
.eqs_seed <- function(seed,fun) {
 had<-exists('.Random.seed',.GlobalEnv,inherits=FALSE)
 if(had)old<-get('.Random.seed',.GlobalEnv)
 kinds<-RNGkind()
 on.exit({do.call(RNGkind,as.list(kinds));if(had)assign('.Random.seed',old,.GlobalEnv)
  else if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)})
 RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed);fun()
}
.eqs_ids <- function(d,cols) {
 .eqs_assert(is.data.frame(d)&&!anyDuplicated(names(d))&&all(cols%in%names(d)),'Missing/duplicated identifier columns')
 for(n in cols).eqs_assert(!anyNA(d[[n]])&&all(nzchar(as.character(d[[n]])))&&
  !any(grepl('\r',as.character(d[[n]]),fixed=TRUE)),'Invalid identifier: ',n)
 .eqs_assert(!anyDuplicated(do.call(paste,c(d[cols],sep='\r'))),'Duplicated test identifiers')
}
.eqs_groups <- function(d,cols) split(seq_len(nrow(d)),do.call(paste,c(d[cols],sep='\r')))
.eqs_measure <- function(z,ii) {
 v<-z[ii,,drop=FALSE];nt<-sum(v$injected);hit<-sum(v$correct_detection)
 calls<-sum(v$called);fp<-sum(v$called&!v$injected);wrong<-sum(v$wrong_sign_call)
 ne<-sum(v$expected_detectable);no<-sum(v$observed_detectable)
 re<-sum(v$correct_detection&v$expected_detectable);ro<-sum(v$correct_detection&v$observed_detectable)
 other<-v$null_kind=='other_type_of_injected_gene';null<-v$null_kind=='fully_null_gene'
 data.frame(n_tests=nrow(v),n_injected=nt,n_injected_all=sum(z$injected),correct_detected=hit,
  calls=calls,false_calls=fp,wrong_sign_calls=wrong,null_tests=sum(!v$injected),
  other_type_null_tests=sum(other),other_type_false_calls=sum(v$called&other),
  fully_null_tests=sum(null),fully_null_false_calls=sum(v$called&null),
  expected_detectable=ne,observed_detectable=no,recovered_expected=re,recovered_observed=ro,
  detection_rate=.eqs_div(hit,nt),recall_all_injected=.eqs_div(hit,sum(z$injected)),
  sensitivity_expected_detectable=.eqs_div(re,ne),sensitivity_observed_detectable=.eqs_div(ro,no),
  FDP=fp/max(1,calls),incorrect_call_fraction=(fp+wrong)/max(1,calls))
}
.eqs_pool <- function(d,keys) {
 sums<-c('n_tests','n_injected','n_injected_all','correct_detected','calls','false_calls','wrong_sign_calls',
  'null_tests','other_type_null_tests','other_type_false_calls','fully_null_tests','fully_null_false_calls',
  'expected_detectable','observed_detectable','recovered_expected','recovered_observed')
 out<-lapply(.eqs_groups(d,keys),function(ii){
  z<-d[ii,,drop=FALSE];a<-z[1,keys,drop=FALSE];a$n_replicates<-length(unique(z$replicate))
  for(n in sums)a[[n]]<-sum(z[[n]])
  a$detection_rate<-.eqs_div(a$correct_detected,a$n_injected)
  a$recall_all_injected<-.eqs_div(a$correct_detected,a$n_injected_all)
  a$sensitivity_expected_detectable<-.eqs_div(a$recovered_expected,a$expected_detectable)
  a$sensitivity_observed_detectable<-.eqs_div(a$recovered_observed,a$observed_detectable)
  a$pooled_FDP<-a$false_calls/max(1,a$calls)
  a$pooled_incorrect_call_fraction<-(a$false_calls+a$wrong_sign_calls)/max(1,a$calls)
  a$mean_replicate_FDP<-mean(z$FDP)
  p<-z$detection_rate[is.finite(z$detection_rate)]
  a$replicate_rate_min<-if(length(p))min(p)else NA_real_
  a$replicate_rate_max<-if(length(p))max(p)else NA_real_;a
 });out<-do.call(rbind,out);rownames(out)<-NULL;out
}
# Full-panel groups use saved quintiles. Matched groups rank only nominated tests.
# Missing scores never determine the BH family; they only affect ranking.
.eqs_selections <- function(z,scope,methods,retentions=c(.2,.4)) {
 cols<-paste(scope,methods,'score',sep='_')
 common<-rowSums(!is.finite(as.matrix(z[cols])))==0L
 out<-list();k<-0L
 for(m in methods){
  q<-z[[paste(scope,m,'quintile',sep='_')]];score<-z[[paste(scope,m,'score',sep='_')]]
  groups<-c(list(all=seq_len(nrow(z))),stats::setNames(lapply(1:5,function(b)which(q==b)),paste0('Q',1:5)),
   list(top20=which(q==5),top40=which(q>=4)))
  for(g in names(groups)) {k<-k+1L;out[[k]]<-list(method=m,policy='full_panel',group=g,ii=groups[[g]])}
  eligible<-which(common);eligible<-eligible[order(-score[eligible],z$gene[eligible],method='radix')]
  for(f in retentions){k<-k+1L;out[[k]]<-list(method=m,policy='matched_tested_budget',
   group=paste0('top',round(100*f)),ii=head(eligible,floor(length(eligible)*f+1e-9)))}
 };out
}
.eqs_report <- function(d,methods,root,retentions=c(.2,.4),patterns=FALSE) {
 keys<-c('scenario','engine','maf','signal_ratio','replicate','cell_type')
 rows<-membership<-shifts<-list();nr<-nm<-nh<-0L
 for(ix in .eqs_groups(d,keys)) {
  z<-d[ix,,drop=FALSE];sels<-list()
  for(scope in c('baseline','scenario')){
   sels[[scope]]<-.eqs_selections(z,scope,methods,retentions)
   for(s in sels[[scope]]){
    strata<-if(patterns)c('all',sort(unique(z$pattern)))else'all'
    for(pattern in strata){
     jj<-if(pattern=='all')seq_len(nrow(z))else which(z$pattern==pattern)
     pos<-match(intersect(s$ii,jj),jj)
     nr<-nr+1L;rows[[nr]]<-cbind(z[1,keys,drop=FALSE],data.frame(score_scope=scope,method=s$method,
      policy=s$policy,group=s$group,pattern=pattern),.eqs_measure(z[jj,,drop=FALSE],pos))
    }
    if(z$engine[1]=='CDSeq2'&&s$group%in%c('top20','top40')){
     nm<-nm+1L;membership[[nm]]<-cbind(z[c('scenario','gene','cell_type')],
      data.frame(score_scope=scope,method=s$method,policy=s$policy,group=s$group,retained=seq_len(nrow(z))%in%s$ii))
    }
   }
  }
  if(z$engine[1]=='CDSeq2')for(i in seq_along(sels$baseline)){
   a<-sels$baseline[[i]];b<-sels$scenario[[i]]
   if(!a$group%in%c('top20','top40'))next
   added<-setdiff(b$ii,a$ii);lost<-setdiff(a$ii,b$ii)
   nh<-nh+1L;shifts[[nh]]<-cbind(z[1,setdiff(keys,'engine'),drop=FALSE],data.frame(method=a$method,
    policy=a$policy,group=a$group,n_baseline=length(a$ii),n_scenario=length(b$ii),
    overlap=length(intersect(a$ii,b$ii)),Jaccard=.eqs_div(length(intersect(a$ii,b$ii)),length(union(a$ii,b$ii))),
    added_correct=sum(z$correct_detection[added]),lost_correct=sum(z$correct_detection[lost]),
    added_false=sum(z$called[added]&!z$injected[added]),lost_false=sum(z$called[lost]&!z$injected[lost])))
  }
 }
 byrep<-do.call(rbind,rows);rownames(byrep)<-NULL
 pooled<-.eqs_pool(byrep,c('engine','maf','signal_ratio','cell_type','score_scope','method','policy','group','pattern'))
 .eqs_csv(byrep,file.path(root,'detection_by_replicate.csv'));.eqs_csv(pooled,file.path(root,'detection_summary.csv'))
 .eqs_csv(do.call(rbind,membership),file.path(root,'selection_membership.csv'))
 .eqs_csv(do.call(rbind,shifts),file.path(root,'baseline_to_scenario_changes.csv'))
 .eqs_csv(pooled[pooled$engine=='CDSeq2'&pooled$pattern=='all'&pooled$group%in%c('all','top20','top40'),],
  file.path(root,'primary_comparison.csv'))
 invisible(list(byrep=byrep,pooled=pooled))
}
.eqs_validate_details <- function(d,methods) {
 .eqs_ids(d,c('scenario','engine','gene','cell_type'))
 req<-c('maf','signal_ratio','replicate','p','q','effect','injected','called','correct_detection','wrong_sign_call',
  'expected_detectable','observed_detectable','null_kind')
 .eqs_assert(nrow(d)>0&&all(req%in%names(d)),'Incomplete test details')
 for(n in c('injected','called','correct_detection','wrong_sign_call','expected_detectable','observed_detectable'))
  .eqs_assert(is.logical(d[[n]])&&!anyNA(d[[n]]),'Invalid logical column: ',n)
 for(n in c('p','q')).eqs_assert(is.numeric(d[[n]])&&all(is.na(d[[n]])|(is.finite(d[[n]])&d[[n]]>=0&d[[n]]<=1)),
  'Invalid p/q column: ',n)
 .eqs_assert(all(d$called==(is.finite(d$q)&d$q<=.05)),'Saved called decisions differ from original q<=0.05')
 .eqs_assert(all(!d$correct_detection|(d$called&d$injected))&&
  all(d$correct_detection+d$wrong_sign_call==as.integer(d$called&d$injected)), 'Inconsistent correct/wrong-sign flags')
 for(scope in c('baseline','scenario'))for(m in methods){
  n<-paste(scope,m,c('score','quintile'),sep='_')
  .eqs_assert(all(n%in%names(d))&&all(vapply(d[n],is.numeric,TRUE)),'Saved scores/quintiles missing: ',scope,' / ',m)
  .eqs_assert(!any(is.infinite(d[[n[1]]])),'Infinite score')
  q<-d[[n[2]]];.eqs_assert(all(is.na(q)|q%in%1:5),'Invalid saved quintile')
  .eqs_assert(all(is.na(q)|is.finite(d[[n[1]]])),'A defined quintile has an undefined score')
 }
 invisible(TRUE)
}
