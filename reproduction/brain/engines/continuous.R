# CDSeq2 continuous phenotype benchmark, version 1.0.0
# Source only defines functions. Uses completed direct-aggregation + GAM6 files.
# No CDSeq2/GAM fitting; base R only. No real age/pathology phenotype is implied.
# source(file.choose())
# run_continuous_associations()
# Defaults: 7 types x 3 programs x 3 signal fractions x 3 noise replicates;
# 200 independent Gaussian null phenotypes. Output is exploratory, not a new
# independent validation of the post-hoc GAM6 choice.

.ca_script <- local({
  p <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names=FALSE)
  if(length(p)) normalizePath(tail(p,1),mustWork=TRUE) else NA_character_
})
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

#' Benchmark continuous associations on the existing validation donors.
#' @param run_dir Completed direct aggregation root (contains inputs/, analysis/).
#' @param gam6_dir GAM6 ablation root; default subdirectory of run_dir.
#' @param output_dir New directory. Existing directories are refused to avoid
#'   overwriting evidence. Use a different name for a sensitivity analysis.
#' @param programs Number of randomly selected gene programs per type, default 3.
#' @param program_size Genes per program, default 10. Chosen using training
#'   expression variability only; independent of validation scores/recovery.
#' @param signal_fractions Nominal signal variance fractions: .25,.5,.75.
#' @param replicates Independent noise replicates per program and signal fraction.
#' @param null_replicates Independent Gaussian null traits, default 200.
#' @param retentions Fractions ranked within each type on a common observable
#'   rankability universe. Default .1,.2,...,1. Includes .2 for main summaries.
#' @param seed Reproducible seed. Caller RNG kind/state are restored.
#' @param covariates Optional finite numeric samples x covariates matrix with
#'   rownames exactly matching validation sample IDs. No intercept column.
#'   Default intercept-only; the same adjustment is used for truth and estimates.
#' @param bh_threshold Global BH reporting threshold, default .05.
#' @details Expression is log2(within-type panel CPM+.5). A program's validation
#'   PC1 is a synthetic expression-linked trait signal, not age, not an injected
#'   DEG, and not a causal effect. Phenotype = sqrt(h)*standardized PC1 +
#'   sqrt(1-h)*independent N(0,1) noise; h is nominal, not exact sample R-squared.
#'   Correlated genes in other types may truly associate; they are NOT labeled
#'   false positives. Observed-component BH agreement is NOT biological FDR.
#'   Counts are held fixed; null summaries measure repeated-trait error on these
#'   donors, not uncertainty over cohorts or full biological generalization.
#' @return Invisibly output directory. Saves traits, test tables, comparisons,
#'   null diagnostics, plots, provenance, and a private donor/sample map.
run_continuous_associations <- function(
    run_dir='~/Downloads/CDSeq2_direct_aggregation_native_v1',
    gam6_dir=NULL,output_dir=NULL,programs=3L,program_size=10L,
    signal_fractions=c(.25,.5,.75),replicates=3L,null_replicates=200L,
    retentions=seq(.1,1,.1),seed=20261010L,covariates=NULL,bh_threshold=.05) {
  for(nm in c('programs','program_size','replicates','null_replicates','seed')) {
    v<-get(nm);.ca_assert(is.numeric(v)&&length(v)==1L&&is.finite(v)&&v>=1&&v==floor(v)&&v<=.Machine$integer.max,'Invalid ',nm)
  }
  .ca_assert(is.numeric(signal_fractions)&&length(signal_fractions)>0&&!anyNA(signal_fractions)&&all(signal_fractions>0&signal_fractions<1)&&!anyDuplicated(signal_fractions),'Invalid signal fractions')
  .ca_assert(is.numeric(retentions)&&length(retentions)>0&&all(is.finite(retentions))&&all(retentions>0&retentions<=1),'Invalid retentions')
  .ca_assert(length(bh_threshold)==1L&&is.finite(bh_threshold)&&bh_threshold>0&&bh_threshold<1,'Invalid BH threshold')
  retentions<-sort(unique(retentions))
  run_dir<-normalizePath(path.expand(run_dir),mustWork=TRUE)
  if(is.null(gam6_dir))gam6_dir<-file.path(run_dir,'gam6_no_composition_loo_r2')
  gam6_dir<-normalizePath(path.expand(gam6_dir),mustWork=TRUE)
  if(is.null(output_dir))output_dir<-file.path(run_dir,'continuous_associations_v1')
  output_dir<-path.expand(output_dir)
  .ca_assert(!dir.exists(output_dir)&&!file.exists(output_dir),'Output exists; choose a NEW output_dir: ',output_dir)
  files<-c(input=file.path(run_dir,'inputs','direct_inputs.rds'),
    estimate=file.path(run_dir,'analysis','validation_estimates.rds'),
    gam7=file.path(run_dir,'analysis','validation_results.rds'),
    gam6=file.path(gam6_dir,'ablation_results.rds'))
  message('Checking counts, donor separation and saved scoring results')
  inp<-.ca_read(files['input']);est<-.ca_read(files['estimate'])
  seven<-.ca_read(files['gam7']);six<-.ca_read(files['gam6'])
  .ca_assert(!any(inp$validation$donors%in%c(inp$training$donors,inp$reference_donors)),'Validation donor overlap')
  .ca_assert(length(inp$validation$donors)==dim(inp$validation$truth)[3]&&!anyDuplicated(inp$validation$donors),'Expected one validation sample per biological donor')
  .ca_assert(identical(dimnames(est$counts),dimnames(inp$validation$truth)),'Estimate/truth axes differ')
  .ca_assert(max(abs(apply(est$counts,c(1,3),sum)-inp$validation$bulk)/pmax(inp$validation$bulk,1))<1e-6,'Allocated counts fail conservation')
  .ca_assert(max(abs(apply(inp$validation$truth,c(1,3),sum)-inp$validation$bulk))<1e-8,'Bulk is not direct aggregation')
  tr<-.ca_expression(inp$training$truth);tt<-.ca_expression(inp$validation$truth);ee<-.ca_expression(est$counts)
  keys<-tt$keys;samples<-colnames(tt$x);types<-unique(keys$cell_type);n<-length(samples)
  .ca_assert(identical(.ca_key(tr$keys),.ca_key(keys)),'Training/validation panels differ')
  align<-function(d) {
    .ca_assert(is.data.frame(d)&&!anyDuplicated(.ca_key(d))&&setequal(.ca_key(d),.ca_key(keys)),'Prediction keys differ')
    d[match(.ca_key(keys),.ca_key(d)),,drop=FALSE]
  }
  p7<-align(seven$predictions);p6<-align(six$validation)
  rankings<-data.frame(GAM6=p6$variation_score,GAM7=p7$variation_score,
    contribution=p7$contribution_mean,expected_count=p7$expected_log_mean_count)
  common<-rowSums(!is.finite(as.matrix(rankings)))==0L
  C<-matrix(1,n,1,dimnames=list(samples,'intercept'))
  if(!is.null(covariates)) {
    .ca_assert(is.matrix(covariates)&&is.numeric(covariates)&&!is.null(rownames(covariates))&&
      !anyDuplicated(rownames(covariates))&&setequal(rownames(covariates),samples)&&all(is.finite(covariates)), 'covariates must be finite numeric with unique sample rownames')
    C<-cbind(C,covariates[samples,,drop=FALSE])
  }
  cq<-qr(C);df<-n-cq$rank-1L
  .ca_assert(cq$rank==ncol(C)&&df>=3,'Covariate design singular or insufficient residual degrees of freedom')
  tx<-.ca_residual(tt$x,cq);ex<-.ca_residual(ee$x,cq)
  # Freeze ranking before any phenotype is generated. Undefined truth never
  # determines rankability; report missing reference associations explicitly.
  sets<-list();si<-0L
  for(ct in types)for(method in names(rankings)) {
    ix<-which(keys$cell_type==ct&common)
    ix<-ix[order(-rankings[[method]][ix],keys$gene[ix],method='radix')]
    for(f in retentions) {
      si<-si+1L;sets[[si]]<-list(cell_type=ct,method=method,retention=f,
        rankable=length(ix),ix=head(ix,floor(length(ix)*f+1e-9)))
    }
  }
  had<-exists('.Random.seed',.GlobalEnv,inherits=FALSE)
  if(had)old<-get('.Random.seed',.GlobalEnv)
  kinds<-RNGkind()
  on.exit({do.call(RNGkind,as.list(kinds));if(had)assign('.Random.seed',old,.GlobalEnv)
    else if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)},add=TRUE)
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed)
  dir.create(output_dir,recursive=TRUE);output_dir<-normalizePath(output_dir)
  dir.create(file.path(output_dir,'signal_tests'))
  csv<-function(x,nm)utils::write.csv(x,file.path(output_dir,nm),row.names=FALSE,na='')
  settings<-list(version='continuous_associations_1.0.0',source_md5=if(!is.na(.ca_script)).ca_md5(.ca_script)else NA_character_,
    inputs=files,input_md5=tools::md5sum(files),seed=seed,programs=programs,program_size=program_size,
    signal_fractions=signal_fractions,replicates=replicates,null_replicates=null_replicates,
    retentions=retentions,bh_threshold=bh_threshold,covariates=C,
    expression='log2(within-type panel CPM + 0.5)',family_size=nrow(keys),
    gam6_model_id=six$model_id,gam7_model_id=seven$model_id,
    status='exploratory_after_GAM6_validation_review',R=R.version.string)
  .ca_save(settings,file.path(output_dir,'settings.rds'))
  csv(data.frame(sample=samples,donor=inp$validation$donors),'donor_map_PRIVATE.csv')
  csv(cbind(keys,rankings,common_rankable=common),'frozen_rankings.csv')
  .ca_save(sets,file.path(output_dir,'selection_sets.rds'))
  records<-traits<-program_records<-list();counter<-0L
  for(ct in types)for(pr in seq_len(programs)) {
    # Prefer upper half of positive training variances; no validation recovery,
    # score or phenotype outcome is used to select these genes.
    ix<-which(tr$keys$cell_type==ct)
    vv<-apply(tr$x[ix,,drop=FALSE],1,stats::var)
    positive<-is.finite(vv)&vv>1e-12
    .ca_assert(sum(positive)>=program_size,'Too few variable training genes: ',ct)
    pool<-ix[positive&vv>=stats::median(vv[positive])]
    if(length(pool)<program_size)pool<-ix[positive]
    chosen<-sample(pool,program_size,replace=FALSE)
    Z<-t(tx[chosen,,drop=FALSE]);sdz<-apply(Z,2,stats::sd)
    usable<-is.finite(sdz)&sdz>1e-10
    .ca_assert(any(usable),'Selected program has no validation variation: ',ct)
    Z<-sweep(Z[,usable,drop=FALSE],2,sdz[usable],'/')
    pc<-svd(Z,nu=1,nv=1);sgn<-sign(pc$v[which.max(abs(pc$v[,1])),1])
    signal<-pc$u[,1]*sgn;signal<-signal/as.numeric(stats::sd(signal))
    pid<-paste0(ct,'_program',pr)
    program_records[[pid]]<-data.frame(program=pid,cell_type=ct,gene=keys$gene[chosen],used_in_PC1=usable)
    for(h in signal_fractions)for(rep in seq_len(replicates)) {
      counter<-counter+1L;id<-sprintf('trait_%03d',counter)
      y<-sqrt(h)*signal+sqrt(1-h)*stats::rnorm(n)
      yr<-drop(qr.resid(cq,y));a<-.ca_associate(tx,yr,df);b<-.ca_associate(ex,yr,df)
      traits[[id]]<-list(id=id,source_type=ct,program=pid,signal_fraction=h,replicate=rep,
        phenotype=stats::setNames(y,samples),signal=stats::setNames(signal,samples))
      .ca_save(cbind(keys,data.frame(reference_r=a$r,reference_p=a$p,reference_q=a$q,
        estimated_r=b$r,estimated_p=b$p,estimated_q=b$q)),file.path(output_dir,'signal_tests',paste0(id,'.rds')))
      rows<-lapply(sets,function(z){
        j<-z$ix;defined<-is.finite(a$r[j])&is.finite(b$r[j])
        ac<-is.finite(a$q[j])&a$q[j]<=bh_threshold;bc<-is.finite(b$q[j])&b$q[j]<=bh_threshold
        agree<-ac&bc&sign(a$r[j])==sign(b$r[j]);agree[is.na(agree)]<-FALSE
        source_ix<-which(keys$cell_type==z$cell_type&common)
        total_ref<-sum(is.finite(a$q[source_ix])&a$q[source_ix]<=bh_threshold)
        data.frame(trait=id,source_type=ct,program=pid,signal_fraction=h,replicate=rep,
          evaluated_type=z$cell_type,method=z$method,retention=z$retention,rankable=z$rankable,selected=length(j),
          paired_defined=sum(defined),association_r_MAE=.ca_mean(abs(a$r[j]-b$r[j])),
          association_r_concordance=.ca_cor(a$r[j],b$r[j]),
          reference_calls=sum(ac),estimated_calls=sum(bc),same_sign_overlap=sum(agree),
          reference_agreement_among_calls=.ca_div(sum(agree),sum(bc)),
          reference_recall_within_selected=.ca_div(sum(agree),sum(ac)),
          reference_recall_full_type=.ca_div(sum(agree),total_ref),
          reference_unavailable_among_calls=sum(bc&!is.finite(a$q[j])))
      })
      records[[counter]]<-do.call(rbind,rows)
      if(counter%%9L==0L)message('Completed ',counter,' signal phenotypes')
    }
  }
  results<-do.call(rbind,records);csv(results,'association_comparisons.csv')
  .ca_save(traits,file.path(output_dir,'synthetic_phenotypes.rds'))
  csv(do.call(rbind,program_records),'program_genes.csv')
  # Aggregate only within source type, signal strength and evaluated type.
  # Replicates share donors/programs; these averages are not independent-cohort CIs.
  metrics<-c('association_r_MAE','association_r_concordance','reference_agreement_among_calls',
    'reference_recall_within_selected','reference_recall_full_type','estimated_calls','same_sign_overlap')
  summary<-stats::aggregate(results[metrics],results[c('source_type','signal_fraction','evaluated_type','method','retention')],.ca_mean)
  csv(summary,'association_summary.csv')
  csv(summary[abs(summary$retention-.2)<1e-9,],'association_top20.csv')
  # All-null traits: every rejection is false. Global-null FDR equals P(any call).
  nulls<-vector('list',null_replicates);null_y<-matrix(NA_real_,n,null_replicates,dimnames=list(samples,NULL))
  for(bi in seq_len(null_replicates)) {
    y<-stats::rnorm(n);null_y[,bi]<-y;yr<-drop(qr.resid(cq,y))
    aa<-.ca_associate(tx,yr,df);bb<-.ca_associate(ex,yr,df)
    rows<-list();zcount<-0L
    for(engine in c('observed_components','CDSeq2')) {
      z<-if(engine=='CDSeq2')bb else aa
      call<-is.finite(z$q)&z$q<=bh_threshold
      zcount<-zcount+1L;rows[[zcount]]<-data.frame(replicate=bi,engine=engine,cell_type='ALL',method='unfiltered',retention=1,
        selected=length(call),false_calls=sum(call),any_false_call=any(call))
      for(ss in sets) {
        zcount<-zcount+1L;rows[[zcount]]<-data.frame(replicate=bi,engine=engine,cell_type=ss$cell_type,
          method=ss$method,retention=ss$retention,selected=length(ss$ix),false_calls=sum(call[ss$ix]),any_false_call=any(call[ss$ix]))
      }
    }
    nulls[[bi]]<-do.call(rbind,rows)
    if(bi%%50L==0L)message('Completed ',bi,' null phenotypes')
  }
  nulls<-do.call(rbind,nulls);csv(nulls,'null_replicate_results.csv')
  ns<-stats::aggregate(nulls[c('false_calls','any_false_call')],nulls[c('engine','cell_type','method','retention')],mean)
  names(ns)[names(ns)=='false_calls']<-'mean_false_calls'
  names(ns)[names(ns)=='any_false_call']<-'probability_any_false_call'
  ns$null_replicates<-null_replicates
  ns$Monte_Carlo_SE<-sqrt(ns$probability_any_false_call*(1-ns$probability_any_false_call)/null_replicates)
  csv(ns,'null_summary.csv');.ca_save(null_y,file.path(output_dir,'null_phenotypes.rds'))
  grDevices::pdf(file.path(output_dir,'association_retention.pdf'),width=10,height=7)
  tryCatch({
    colors<-c(GAM6='#0072B2',GAM7='#D55E00',contribution='#009E73',expected_count='#CC79A7')
    for(ct in types)for(h in signal_fractions) {
      z<-summary[summary$source_type==ct&summary$evaluated_type==ct&summary$signal_fraction==h,]
      graphics::par(mfrow=c(1,2),mar=c(4,4,4,1))
      for(metric in c('association_r_MAE','reference_recall_full_type')) {
        finite<-z[[metric]][is.finite(z[[metric]])];lim<-if(length(finite))range(c(0,finite))else c(0,1)
        if(diff(lim)==0)lim<-c(0,1)
        graphics::plot(NA,xlim=c(0,1),ylim=lim,xlab='Fraction retained within type',ylab=metric,
          main=paste(ct,'nominal signal',h))
        for(m in names(colors)){v<-z[z$method==m,];v<-v[order(v$retention),];graphics::lines(v$retention,v[[metric]],col=colors[m],lwd=2)}
        graphics::legend('topright',legend=names(colors),col=colors,lwd=2,cex=.75,bty='n')
      }
    }
  },finally=grDevices::dev.off())
  capture.output(utils::sessionInfo(),file=file.path(output_dir,'sessionInfo.txt'))
  writeLines(c('Completed continuous phenotype association benchmark.',
    paste(counter,'signal phenotypes;',null_replicates,'independent Gaussian null phenotypes.'),
    'GAM6/GAM7 and support rankings were frozen before phenotype generation. No fitting.',
    'Phenotypes are expression-program PC1 plus Gaussian noise; they are NOT simulated age.',
    'Default tests use log2(within-type panel CPM + 0.5) and an intercept only.',
    'BH computed over the complete gene/type family separately for each phenotype and expression source, BEFORE ranking filters.',
    'Reference associations use noisy observed single-nucleus aggregates, not error-free biological truth.',
    'Reference call agreement is NOT precision/FDR against known causal effects.',
    'Independent-null results estimate error conditional on the current expression matrices.',
    'No biological-cohort confidence intervals: programs/noise repetitions share the same donors.',
    'Exploratory: GAM6 was selected after examining these validation donors.',
    'Read association_top20.csv alongside full retention curves and null_summary.csv.'),file.path(output_dir,'SUMMARY.txt'))
  writeLines(as.character(Sys.time()),file.path(output_dir,'COMPLETED'))
  message('COMPLETED: ',output_dir)
  invisible(output_dir)
}
