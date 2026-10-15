# Three-type nominated-variant simulation v1.0.0. Source defines functions only.
.eq3_file<-local({p<-unlist(lapply(sys.frames(),function(z)z$ofile),use.names=FALSE)
 if(length(p))normalizePath(tail(p,1),mustWork=TRUE)else NA_character_})
if(!is.na(.eq3_file))sys.source(file.path(dirname(.eq3_file),'eqtl_helpers.R'),environment())
.eq3_methods<-c('GAM6','contribution','expected_count')

#' Calibrated Excit/Inhib/Oligo genetic-association simulation.
#' @param action check/prepare simulate and validate all inputs, no CDSeq2 fits;
#'   probe fits one largest-library donor in a strongest-signal scenario;
#'   run resumes all fits and analyzes; analyze reads existing fits only.
#' @param workers Fresh native PSOCK workers, one CDSeq2 thread each.
#' @param run_dir Completed direct aggregation root; the saved workflow and
#'   installed runtime must match the original experiment exactly.
#' @param gam6_file Explicit user-trained GAM6 RDS; default finds YOUR saved fit
#'   under run_dir/gam6_no_composition_loo_r2. No model is shipped or trained here.
#' @param output_dir New experiment directory; identical calls resume safely.
#' @param target_types Main injection types; defaults Excit, Inhib, Oligo.
#' @param targets_per_quintile Single-type genes per baseline GAM6 quintile/type.
#' @param shared_per_pair Genes per type-pair for EACH shared/opposite pattern.
#' @param offtype_per_type Single-type controls per remaining mixture type.
#' @param null_genes Additional wholly unperturbed nominated genes.
#' @param signal_ratios Standardized genotype-projected logCPM slope increments.
#'   These are not heritabilities or fixed log2 count effects. All affected types,
#'   including other-type controls, use the same requested signal ratio.
#' @param mafs Nominal binomial allele frequencies; default .3.
#' @param replicates Independent genotype/multinomial draws on the SAME donors.
#' @param seed Fixed seed. Use the identical configuration to resume.
#' @param expressed_fraction Minimum nonzero fraction in training and validation
#'   components for candidate targets; used only to construct the benchmark.
#' @param max_abs_beta Maximum absolute calibrated log2 count multiplier slope;
#'   infeasible signals stop before any fits; they are never silently capped.
#' @param max_source_logCPM_sd,max_source_panel_share Construction eligibility:
#'   exclude extreme baseline variability/dominance before drawing genotypes.
#'   These restrictions apply only to injected targets, not the scoring universe.
#' @param exclude_files Optional CSVs containing genes to exclude. NULL records
#'   prior pilot/alltypes/Excit target files that exist. character() excludes none.
#' @return Output path. Default: 254 unique genes, 262 injected gene/type pairs,
#'   seven mixture types, 9 scenarios x 40 donors = 360 CDSeq2 donor fits.
#'   Known components are for construction/evaluation only. Scoring never uses G.
run_eqtl_three_types<-function(action=c('check','prepare','probe','run','analyze'),workers=12L,
 run_dir='~/Downloads/CDSeq2_direct_aggregation_native_v1',gam6_file=NULL,output_dir=NULL,
 target_types=c('Excit','Inhib','Oligo'),targets_per_quintile=10L,shared_per_pair=8L,
 offtype_per_type=4L,null_genes=40L,signal_ratios=c(0,.5,.8),mafs=.3,
 replicates=3L,seed=20261015L,expressed_fraction=.8,max_abs_beta=8,
 max_source_logCPM_sd=2,max_source_panel_share=.02,exclude_files=NULL) {
 action<-match.arg(action);.eqs_assert(!is.na(.eq3_file),'Source step2_three_type_simulation.R first')
 args<-as.list(environment());args<-args[setdiff(names(formals(run_eqtl_three_types)),character())]
 args$run_dir<-normalizePath(path.expand(run_dir),mustWork=TRUE)
 if(is.null(gam6_file))args$gam6_file<-file.path(args$run_dir,'gam6_no_composition_loo_r2','GAM6_model.rds')
 if(is.null(output_dir))args$output_dir<-file.path(args$run_dir,'eqtl_three_types_v1')
 args$gam6_file<-normalizePath(path.expand(args$gam6_file),mustWork=TRUE)
 args$output_dir<-path.expand(args$output_dir)
 job<-tempfile(fileext='.rds');launcher<-tempfile(fileext='.R');on.exit(unlink(c(job,launcher)),add=TRUE)
 saveRDS(list(args=args,paths=.libPaths()),job)
 writeLines(c(paste0('source(',encodeString(.eq3_file,quote='"'),')'),
  paste0('job <- readRDS(',encodeString(job,quote='"'),')'),'.libPaths(job$paths)',
  'do.call(.eq3_run,job$args)'),launcher)
 status<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(launcher)))
 .eqs_assert(status==0L,'Three-type experiment stopped; see diagnostic above. Verified caches remain; no settings are changed automatically.')
 invisible(args$output_dir)
}

.eq3_score<-function(f,gam,e) {
 p<-e$predict_direct_variation(f,gam)
 d<-data.frame(f[c('gene','cell_type')],GAM6=p$variation_score,contribution=f$contribution_mean,
  expected_count=f$expected_log_mean_count,GAM6_outside_training_range=p$outside_training_range)
 d$common_rankable<-rowSums(!is.finite(as.matrix(d[.eq3_methods])))==0L
 for(m in .eq3_methods){d[[paste0(m,'_quintile')]]<-NA_integer_
  for(ct in unique(d$cell_type)){
   ii<-which(d$cell_type==ct&d$common_rankable);ii<-ii[order(d[[m]][ii],d$gene[ii],method='radix')]
   if(length(ii))d[[paste0(m,'_quintile')]][ii]<-pmin(5L,floor((seq_along(ii)-1)*5/length(ii))+1L)
  }
 };d
}
.eq3_equal<-function(a,b,tol=1e-8) {
 length(a)==length(b)&&identical(is.na(as.vector(a)),is.na(as.vector(b)))&&
 all(abs(as.numeric(a[!is.na(a)])-as.numeric(b[!is.na(b)]))<=tol)
}
.eq3_design<-function(inp,base,target_types,nq,np,no,nn,prevalence,max_sd,max_share,excluded,seed,root) {
 tr<-inp$training$truth;va<-inp$validation$truth;genes<-dimnames(tr)[[1]];types<-dimnames(tr)[[2]]
 eligible<-bins<-list();support<-counts<-list()
 for(k in seq_along(types)){
  ct<-types[k];a<-matrix(tr[,k,],nrow=length(genes));v<-matrix(va[,k,],nrow=length(genes))
  z<-base[base$cell_type==ct,];z<-z[match(genes,z$gene),]
  lv<-colSums(v);sdlog<-apply(log2(sweep(v,2,lv,'/')*1e6+.5),1,stats::sd)
  panelshare<-apply(sweep(v,2,lv,'/'),1,max)
  ok<-rowMeans(a)>=1&apply(a,1,stats::sd)>0&rowMeans(a>0)>=prevalence&rowMeans(v>0)>=prevalence&
   is.finite(sdlog)&sdlog>1e-10&sdlog<=max_sd&panelshare<=max_share&z$common_rankable&!genes%in%excluded
  eligible[[ct]]<-genes[ok];bins[[ct]]<-stats::setNames(z$GAM6_quintile,genes)
  support[[ct]]<-data.frame(gene=genes,cell_type=ct,training_mean=rowMeans(a),
   training_nonzero=rowMeans(a>0),validation_nonzero=rowMeans(v>0),validation_logCPM_sd=sdlog,
   validation_max_panel_share=panelshare,eligible=ok)
  for(q in 1:5)counts[[length(counts)+1L]]<-data.frame(cell_type=ct,quintile=q,
   candidates=sum(ok&z$GAM6_quintile==q,na.rm=TRUE),single_type_requested=if(ct%in%target_types)nq else 0L)
 }
 .eqs_csv(do.call(rbind,support),file.path(root,'candidate_support.csv'))
 .eqs_csv(do.call(rbind,counts),file.path(root,'candidate_counts.csv'))
 pairs<-utils::combn(target_types,2,simplify=FALSE)
 out<-.eqs_seed(seed,function(){
  rows<-effects<-list();used<-character();serial<-0L
  add<-function(chosen,cts,pattern,q=NA_integer_,opposite=FALSE){
   signs<-sample(rep(c(-1,1),length.out=length(chosen)))
   for(j in seq_along(chosen)){
    serial<<-serial+1L;id<-sprintf('THREE_SYNVAR_%04d',serial)
    rows[[length(rows)+1L]]<<-data.frame(variant=id,gene=chosen[j],pattern=pattern,
     source_set=if(length(cts))paste(cts,collapse=';')else'NONE',design_quintile=q)
    if(length(cts))for(k in seq_along(cts))effects[[length(effects)+1L]]<<-data.frame(
     variant=id,gene=chosen[j],cell_type=cts[k],effect_sign=signs[j]*(if(opposite&&k==2L)-1 else 1))
   };used<<-c(used,chosen)
  }
  # Shared pairs are allocated first to protect their more restricted intersections.
  for(pattern in c('shared_same','shared_opposite'))for(pair in pairs){
   free<-setdiff(intersect(eligible[[pair[1]]],eligible[[pair[2]]]),used)
   .eqs_assert(length(free)>=np,'Insufficient shared candidates: ',paste(pair,collapse='/'))
   add(sample(free,np),pair,pattern,opposite=pattern=='shared_opposite')
  }
  # Allocate the rarest remaining type/quintile pools first; retain lexical ties.
  slots<-do.call(rbind,lapply(target_types,function(ct)data.frame(cell_type=ct,quintile=1:5)))
  slots$n<-vapply(seq_len(nrow(slots)),function(j){ct<-slots$cell_type[j]
   sum(bins[[ct]][setdiff(eligible[[ct]],used)]==slots$quintile[j],na.rm=TRUE)},0L)
  slots<-slots[order(slots$n,slots$cell_type,slots$quintile,method='radix'),]
  for(j in seq_len(nrow(slots))){ct<-slots$cell_type[j];q<-slots$quintile[j]
   free<-setdiff(eligible[[ct]][which(bins[[ct]][eligible[[ct]]]==q)],used)
   .eqs_assert(length(free)>=nq,'Insufficient unused single-type candidates: ',ct,' Q',q,
    '. No eligibility rules relaxed; inspect candidate_counts.csv.')
   add(sample(free,nq),ct,'single_type',q)
  }
  for(ct in setdiff(types,target_types)){
   free<-setdiff(eligible[[ct]],used);.eqs_assert(length(free)>=no,'Insufficient other-type controls: ',ct)
   add(sample(free,no),ct,'other_type_control')
  }
  # Null targets are eligible by training bulk support, independently of scores/outcomes.
  free<-setdiff(genes[rowMeans(inp$training$bulk)>=1],c(used,excluded))
  .eqs_assert(length(free)>=nn,'Insufficient unused null genes');add(sample(free,nn),character(),'fully_null')
  list(design=do.call(rbind,rows),effects=do.call(rbind,effects))
 })
 rownames(out$design)<-rownames(out$effects)<-NULL
 .eqs_assert(!anyDuplicated(out$design$gene)&&!any(out$design$gene%in%excluded),'Target design overlap')
 .eqs_ids(out$effects,c('gene','cell_type'))
 out
}

.eq3_make<-function(inp,design,effects,maf,strength,replicate,seed,max_beta) {
 a<-inp$validation$truth;genes<-dimnames(a)[[1]];types<-dimnames(a)[[2]];samples<-dimnames(a)[[3]]
 lib<-apply(a,c(2,3),sum)
 .eqs_assert(all(lib>0&lib<=.Machine$integer.max&lib==floor(lib)),'Invalid component multinomial libraries')
 G<-.eqs_seed(seed,function(){
  z<-matrix(0,length(samples),nrow(design),dimnames=list(samples,design$variant))
  for(j in seq_len(ncol(z))){accepted<-FALSE
   for(attempt in seq_len(10000L)){v<-stats::rbinom(nrow(z),2,maf)
    if(sum(v>0)>=5L&&sum(v==0)>=5L){z[,j]<-v;accepted<-TRUE;break}}
   .eqs_assert(accepted,'Could not draw eligible genotypes')
  };z
 })
 expected<-a*1;nominated<-match(design$gene,genes);reserve<-setdiff(seq_along(genes),nominated)
 .eqs_assert(length(reserve)>0,'No untested reserve genes');cal<-list()
 for(j in seq_len(nrow(effects))){
  ef<-effects[j,];g<-match(ef$gene,genes);k<-match(ef$cell_type,types);x<-G[,ef$variant]
  xc<-x-mean(x);sx<-stats::sd(xc);sxx<-sum(xc^2);cpm<-a[g,k,]/lib[k,]*1e6;y0<-log2(cpm+.5)
  sy<-stats::sd(y0);sg<-ef$effect_sign
  .eqs_assert(is.finite(sy)&&sy>1e-10,'Undefined baseline variability: ',ef$gene,' / ',ef$cell_type)
  target<-strength*sy/sx;coef<-0
  objective<-function(b)sg*sum(xc*(log2(cpm*2^(sg*b*xc)+.5)-y0))/sxx-target
  if(strength>0){
   .eqs_assert(is.finite(objective(max_beta))&&objective(max_beta)>=0,
    'Signal infeasible below max_abs_beta for ',ef$gene,' / ',ef$cell_type,
    '. Preparation stops; no effects capped and no CDSeq2 fits started.')
   coef<-sg*stats::uniroot(objective,c(0,max_beta),tol=1e-10)$root
  }
  expected[g,k,]<-a[g,k,]*2^(coef*xc)
  yy<-log2(expected[g,k,]/lib[k,]*1e6+.5);delta<-sum(xc*(yy-y0))/sxx;achieved<-abs(delta)*sx/sy
  .eqs_assert(abs(achieved-strength)<1e-7,'Signal calibration failed')
  cal[[j]]<-cbind(ef,data.frame(applied_log2_count_beta=coef,baseline_logCPM_sd=sy,genotype_sd=sx,
   baseline_genotype_slope=sum(xc*y0)/sxx,target_abs_slope_increment=target,
   achieved_slope_increment=delta,achieved_signal_ratio=achieved,maximum_count_multiplier=max(2^(coef*xc))))
 }
 minreserve<-1
 for(k in seq_along(types))for(i in seq_along(samples)){
  remaining<-lib[k,i]-sum(expected[nominated,k,i]);den<-sum(a[reserve,k,i])
  .eqs_assert(is.finite(remaining)&&remaining>0&&den>0,'Injection exceeds reserve budget: ',types[k],' / ',samples[i],
   '. No effects shrunk; inspect the design/configuration. New settings require a new output_dir.')
  minreserve<-min(minreserve,remaining/lib[k,i]);expected[reserve,k,i]<-a[reserve,k,i]*remaining/den
 }
 error<-max(abs(apply(expected,c(2,3),sum)-lib)/pmax(lib,1));.eqs_assert(error<1e-12,'Expected library conservation failed')
 affected<-.eqs_key(effects)
 for(g in nominated)for(k in seq_along(types))if(!paste(genes[g],types[k],sep='\r')%in%affected)
  .eqs_assert(identical(as.numeric(expected[g,k,]),as.numeric(a[g,k,])),'Protected non-target changed')
 if(strength==0).eqs_assert(max(abs(expected-a))<1e-7,'No-injection expectation changed')
 realized<-.eqs_seed(seed+100000L,function(){
  z<-array(0,dim(a),dimnames(a))
  for(i in seq_along(samples))for(k in seq_along(types))z[,k,i]<-as.numeric(stats::rmultinom(1L,lib[k,i],expected[,k,i]));z
 })
 .eqs_assert(identical(as.numeric(apply(realized,c(2,3),sum)),as.numeric(lib)),'Realized library conservation failed')
 list(truth=realized,expected=expected,bulk=apply(realized,c(1,3),sum),genotypes=G,design=design,effects=effects,
  strength=strength,maf=maf,replicate=replicate,seed=seed,calibration=do.call(rbind,cal),
  QC=data.frame(expected_library_relative_error=error,min_reserve_fraction=minreserve,
   protected_nominated_profiles_unchanged=TRUE,realized_libraries_conserved=TRUE),
  genotype_QC=data.frame(variant=colnames(G),MAF_nominal=maf,allele_frequency=colMeans(G)/2,
   n0=colSums(G==0),n1=colSums(G==1),n2=colSums(G==2)))
}

.eq3_tests<-function(a,sim) {
 genes<-dimnames(a)[[1]];types<-dimnames(a)[[2]];lib<-apply(a,c(2,3),sum);n<-dim(a)[3];rows<-list()
 .eqs_assert(all(is.finite(a))&&all(a>=0)&&all(lib>0),'Invalid/zero component libraries')
 for(j in seq_len(nrow(sim$design)))for(k in seq_along(types)){
  g<-match(sim$design$gene[j],genes);x<-sim$genotypes[,j];y<-log2(a[g,k,]/lib[k,]*1e6+.5)
  xc<-x-mean(x);yc<-y-mean(y);sxx<-sum(xc^2);slope<-sum(xc*yc)/sxx
  err<-yc-slope*xc;se<-sqrt(sum(err^2)/(n-2)/sxx)
  p<-if(stats::sd(y)<1e-12)NA_real_ else if(se>0)2*stats::pt(abs(slope/se),n-2,lower.tail=FALSE)else 0
  rows[[length(rows)+1L]]<-data.frame(variant=sim$design$variant[j],gene=genes[g],cell_type=types[k],effect=slope,p=p)
 }
 d<-do.call(rbind,rows)
 # One BH adjustment over ALL nominated gene/type hypotheses, before selection.
 d$q<-stats::p.adjust(replace(d$p,!is.finite(d$p),1),'BH');d$q[!is.finite(d$p)]<-NA_real_;d
}
.eq3_evaluate<-function(sim,est,baseline,gam,inp,e,scenario,root) {
 sc<-.eq3_score(e$.da_observables(sim$bulk,inp$reference_mean,est$proportions),gam,e)
 tables<-list(expected=.eq3_tests(sim$expected,sim),observed=.eq3_tests(sim$truth,sim),CDSeq2=.eq3_tests(est$counts,sim))
 for(engine in names(tables)){
  d<-tables[[engine]];j<-match(d$gene,sim$design$gene);k<-match(.eqs_key(d),.eqs_key(sim$effects))
  d$injected<-sim$strength>0&!is.na(k);d$injected_sign<-0
  d$injected_sign[d$injected]<-sim$effects$effect_sign[k[d$injected]]
  d$injected_beta<-0;d$injected_beta[d$injected]<-sim$calibration$applied_log2_count_beta[k[d$injected]]
  d$pattern<-sim$design$pattern[j];d$source_set<-if(sim$strength==0)'NO_INJECTION'else sim$design$source_set[j]
  d$null_kind<-ifelse(d$injected,'injected_target',if(sim$strength==0)'no_injection_control'else
   ifelse(sim$design$source_set[j]=='NONE','fully_null_gene','other_type_of_injected_gene'))
  d$called<-is.finite(d$q)&d$q<=.05
  d$correct_detection<-d$called&d$injected&sign(d$effect)==d$injected_sign
  d$correct_detection[is.na(d$correct_detection)]<-FALSE
  d$wrong_sign_call<-d$called&d$injected&!d$correct_detection
  for(name in c('expected','observed')){
   v<-tables[[name]];d[[paste0(name,'_detectable')]]<-d$injected&is.finite(v$q)&v$q<=.05&sign(v$effect)==d$injected_sign
   d[[paste0(name,'_detectable')]][is.na(d[[paste0(name,'_detectable')]])]<-FALSE
  }
  d$expected_component_slope<-tables$expected$effect
  d$scenario<-scenario;d$engine<-engine;d$maf<-sim$maf;d$signal_ratio<-sim$strength;d$replicate<-sim$replicate
  for(scope in c('baseline','scenario')){
   s<-if(scope=='baseline')baseline else sc;ix<-match(.eqs_key(d),.eqs_key(s));.eqs_assert(!anyNA(ix),'Score keys missing')
   for(m in .eq3_methods){d[[paste(scope,m,'score',sep='_')]]<-s[[m]][ix]
    d[[paste(scope,m,'quintile',sep='_')]]<-s[[paste0(m,'_quintile')]][ix]}
   d[[paste0(scope,'_GAM6_outside_training_range')]]<-s$GAM6_outside_training_range[ix]
  };tables[[engine]]<-d
 }
 .eqs_save(sc,file.path(root,'analysis',paste0(scenario,'_scores.rds')))
 d<-do.call(rbind,tables);rownames(d)<-NULL;d
}

.eq3_attribution<-function(d,root) {
 d<-d[d$engine=='CDSeq2',];sel<-utils::read.csv(file.path(root,'selection_membership.csv'),stringsAsFactors=FALSE)
 key<-function(x)paste(x$scenario,x$gene,x$cell_type,sep='\r')
 ix<-match(key(sel),key(d));.eqs_assert(!anyNA(ix),'Attribution selection keys differ')
 cols<-c('maf','signal_ratio','replicate','pattern','source_set','injected','called','correct_detection','wrong_sign_call')
 s<-cbind(sel,d[ix,cols]);all<-d[c('scenario','gene','cell_type',cols)]
 all$score_scope<-'none';all$method<-'unfiltered';all$policy<-'all';all$group<-'all';all$retained<-TRUE
 s<-rbind(s,all[names(s)]);keys<-c('scenario','maf','signal_ratio','replicate','score_scope','method','policy','group','pattern','source_set','cell_type')
 rows<-lapply(.eqs_groups(s,keys),function(ii){z<-s[ii,];v<-z[z$retained,]
  cbind(z[1,keys,drop=FALSE],data.frame(n_tests_all=nrow(z),retained=nrow(v),
   affected_tests=sum(v$injected),correct_detected=sum(v$correct_detection),
   unaffected_tests=sum(!v$injected),false_calls=sum(v$called&!v$injected),wrong_sign_calls=sum(v$wrong_sign_call)))})
 .eqs_csv(do.call(rbind,rows),file.path(root,'source_to_recipient_by_replicate.csv'))
 pair<-s[s$pattern%in%c('shared_same','shared_opposite')&s$signal_ratio>0,]
 keys2<-c('scenario','maf','signal_ratio','replicate','score_scope','method','policy','group','pattern','source_set','gene')
 rows<-lapply(.eqs_groups(pair,keys2),function(ii){z<-pair[ii,];target<-z$injected
  cbind(z[1,keys2,drop=FALSE],data.frame(affected_types=sum(target),selected_affected_types=sum(target&z$retained),
   recovered_affected_types=sum(target&z$retained&z$correct_detection),
   both_affected_types_recovered=all(z$retained[target]&z$correct_detection[target]),
   false_recipient_calls=sum(z$retained&z$called&!target)))})
 .eqs_csv(do.call(rbind,rows),file.path(root,'shared_effect_recovery_by_gene.csv'))
}

.eq3_plot<-function(z,target_types,root) {
 grDevices::pdf(file.path(root,'three_type_detection.pdf'),width=12,height=8)
 on.exit(grDevices::dev.off(),add=TRUE)
 colors<-c(GAM6='#0072B2',contribution='#D55E00',expected_count='#009E73')
 for(maf in unique(z$maf))for(strength in sort(unique(z$signal_ratio[z$signal_ratio>0]))){
  graphics::par(mfrow=c(2,length(target_types)),mar=c(4,4,3,1),oma=c(3,0,1,0))
  for(scope in c('scenario','baseline'))for(ct in target_types){
   graphics::plot(NA,xlim=c(1,5),ylim=c(0,1),xlab='Full-panel score quintile',ylab='Correct-sign detection rate',
    main=paste(ct,scope,'MAF',maf,'signal',strength))
   for(m in .eq3_methods){v<-z[z$engine=='CDSeq2'&z$maf==maf&z$signal_ratio==strength&z$cell_type==ct&
    z$score_scope==scope&z$method==m&z$policy=='full_panel'&z$pattern=='single_type'&z$group%in%paste0('Q',1:5),]
    j<-match(paste0('Q',1:5),v$group);graphics::lines(1:5,v$detection_rate[j],col=colors[m],type='b',pch=16,lwd=2)
   }
   graphics::legend('topleft',legend=.eq3_methods,col=colors,lty=1,pch=16,cex=.65,bty='n')
  }
  graphics::mtext('Single-type targets; original all-type BH. Replicates reuse donors/genes; no independent-event inference.',side=1,outer=TRUE,cex=.8)
 }
}

.eq3_run<-function(action,workers,run_dir,gam6_file,output_dir,target_types,targets_per_quintile,
 shared_per_pair,offtype_per_type,null_genes,signal_ratios,mafs,replicates,seed,expressed_fraction,max_abs_beta,
 max_source_logCPM_sd,max_source_panel_share,exclude_files) {
 for(nm in c('workers','targets_per_quintile','shared_per_pair','offtype_per_type','null_genes','replicates','seed')){
  v<-get(nm);.eqs_assert(is.numeric(v)&&length(v)==1L&&is.finite(v)&&v>=1&&v==floor(v)&&v<1e8,'Invalid ',nm)
 }
 for(v in c(targets_per_quintile,shared_per_pair,offtype_per_type)).eqs_assert(v%%2==0,'Use even target counts for balanced directions')
 .eqs_assert(replicates<=100&&is.character(target_types)&&length(target_types)>=2L&&length(target_types)<=4L&&
  !anyNA(target_types)&&!anyDuplicated(target_types)&&all(nzchar(target_types)),'Invalid target_types/replicates')
 .eqs_assert(is.numeric(mafs)&&length(mafs)>0&&length(mafs)<=10L&&all(is.finite(mafs))&&all(mafs>=.1&mafs<=.5)&&!anyDuplicated(mafs),'Use unique MAFs in [.1,.5]')
 .eqs_assert(is.numeric(signal_ratios)&&length(signal_ratios)>1L&&length(signal_ratios)<=10L&&all(is.finite(signal_ratios))&&
  all(signal_ratios>=0&signal_ratios<=1)&&0%in%signal_ratios&&any(signal_ratios>0)&&!anyDuplicated(signal_ratios),'Invalid signal_ratios')
 .eqs_assert(is.numeric(expressed_fraction)&&length(expressed_fraction)==1L&&is.finite(expressed_fraction)&&
  expressed_fraction>0&&expressed_fraction<=1,'Invalid expressed_fraction')
 .eqs_assert(is.numeric(max_abs_beta)&&length(max_abs_beta)==1L&&is.finite(max_abs_beta)&&max_abs_beta>0&&max_abs_beta<=16,'Invalid max_abs_beta')
 .eqs_assert(is.numeric(max_source_logCPM_sd)&&length(max_source_logCPM_sd)==1L&&is.finite(max_source_logCPM_sd)&&max_source_logCPM_sd>0,'Invalid max_source_logCPM_sd')
 .eqs_assert(is.numeric(max_source_panel_share)&&length(max_source_panel_share)==1L&&is.finite(max_source_panel_share)&&max_source_panel_share>0&&max_source_panel_share<=1,'Invalid max_source_panel_share')
 .eqs_assert(requireNamespace('mgcv',quietly=TRUE),'Install mgcv in the same R library as the trained scorer')
 workflow<-file.path(run_dir,'workflow.R');e<-new.env(parent=globalenv());sys.source(workflow,e)
 old<-e$.da_load(file.path(run_dir,'experiment.rds'));runtime<-e$.da_runtime()
 .eqs_assert(identical(runtime,old$config$runtime),'Native R/package runtime differs from the direct-aggregation run. Use the same installed package/runtime; do not bypass this guard.')
 .eqs_assert(identical(e$.da_md5(workflow),old$config$script_md5),'Original workflow changed')
 sourcepaths<-c(input=file.path(run_dir,'inputs','direct_inputs.rds'),reference=file.path(run_dir,'inputs','reference.rds'),
  gam=gam6_file,features=file.path(run_dir,'analysis','validation_features_targets.rds'),
  ablation=file.path(run_dir,'gam6_no_composition_loo_r2','ablation_results.rds'))
 inp<-e$.da_load(sourcepaths['input']);gam<-e$.da_load(sourcepaths['gam']);f<-e$.da_load(sourcepaths['features'])
 e$.da_features<-setdiff(e$.da_features,'composition_loo_r2')
 .eqs_assert(identical(gam$features,e$.da_features)&&identical(gam$method,'GAM6_no_composition_loo_r2'),
  'Supply the explicit trained GAM6 ablation model used for this application')
 types<-dimnames(inp$validation$truth)[[2]];genes<-dimnames(inp$validation$truth)[[1]]
 .eqs_assert(all(target_types%in%types)&&setequal(names(gam$models),types),'Target/scorer cell types differ')
 .eqs_assert(identical(dimnames(inp$training$truth)[1:2],dimnames(inp$validation$truth)[1:2]),'Training and validation axes differ')
 e$.da_array(inp$training$truth);e$.da_array(inp$validation$truth)
 for(cohort in c('training','validation')){
  v<-inp[[cohort]]
  .eqs_assert(length(v$donors)==dim(v$truth)[3]&&!anyNA(v$donors)&&!anyDuplicated(v$donors)&&
   all(nzchar(v$donors))&&identical(dimnames(v$bulk),dimnames(v$truth)[c(1,3)])&&
   max(abs(apply(v$truth,c(1,3),sum)-v$bulk))<1e-8,'Invalid direct counts or biological donor IDs: ',cohort)
 }
 .eqs_assert(!any(inp$validation$donors%in%c(inp$training$donors,inp$reference_donors))&&
  !any(inp$training$donors%in%inp$reference_donors),'Reference/training/validation donor overlap')
 .eqs_assert(setequal(gam$donor_ids,inp$training$donors)&&setequal(gam$reference_donor_ids,inp$reference_donors),
  'Scorer biological donor provenance differs')
 replay<-e$predict_direct_variation(f,gam);saved<-e$.da_load(sourcepaths['ablation'])$validation
 .eqs_ids(replay,c('gene','cell_type'));.eqs_ids(saved,c('gene','cell_type'))
 .eqs_assert(setequal(.eqs_key(replay),.eqs_key(saved)),'Saved scorer replay keys differ')
 saved<-saved[match(.eqs_key(replay),.eqs_key(saved)),]
 .eqs_assert(.eq3_equal(replay$predicted_r,saved$predicted_r)&&.eq3_equal(replay$variation_score,saved$variation_score),
  'Saved GAM6 prediction replay failed')
 base<-.eq3_score(f,gam,e)
 if(is.null(exclude_files)){
  candidates<-file.path(run_dir,c('eqtl_pilot_v1','eqtl_alltypes_v2','eqtl_excit_v3'),'nominated_variant_gene_design.csv')
  exclude_files<-candidates[file.exists(candidates)]
 }
 .eqs_assert(is.character(exclude_files)&&!anyNA(exclude_files),'exclude_files must be character paths or NULL')
 if(length(exclude_files))exclude_files<-vapply(path.expand(exclude_files),normalizePath,'',mustWork=TRUE)
 excluded<-unique(unlist(lapply(exclude_files,function(p){z<-utils::read.csv(p,stringsAsFactors=FALSE)
  .eqs_assert('gene'%in%names(z)&&!anyNA(z$gene)&&all(nzchar(z$gene)),'Invalid prior design: ',p);as.character(z$gene)}),use.names=FALSE))
 spec<-list(version='three_type_eqtl_1.0.0',code_md5=tools::md5sum(c(.eq3_file,file.path(dirname(.eq3_file),'eqtl_helpers.R'))),
  workflow_md5=.eqs_md5(workflow),input_md5=tools::md5sum(c(sourcepaths,exclude_files)),runtime=runtime,
  target_types=target_types,targets_per_quintile=targets_per_quintile,shared_per_pair=shared_per_pair,
  offtype_per_type=offtype_per_type,null_genes=null_genes,signal_ratios=signal_ratios,mafs=mafs,
  replicates=replicates,seed=seed,expressed_fraction=expressed_fraction,max_abs_beta=max_abs_beta,
  max_source_logCPM_sd=max_source_logCPM_sd,max_source_panel_share=max_source_panel_share,excluded_genes=excluded,
  primary_policy=list(score_scope='scenario',method='GAM6',cell_type='Excit',policy='matched_tested_budget',retention=.4,BH=.05),
  oracle='expected_components_with_observed_component_comparator',
  generator='protected_nominated_genes_fixed_libraries_standardized_slope_increments_v1')
 id<-.eqs_hash(spec)
 dir.create(output_dir,recursive=TRUE,showWarnings=FALSE);root<-normalizePath(output_dir)
 lock<-file.path(root,'RUNNING.lock');.eqs_assert(dir.create(lock,showWarnings=FALSE),'An active/stale lock exists: ',lock)
 on.exit(unlink(lock,recursive=TRUE),add=TRUE)
 writeLines(c(paste('PID',Sys.getpid()),paste('HOST',Sys.info()[['nodename']])),file.path(lock,'owner.txt'))
 mp<-file.path(root,'experiment.rds')
 if(file.exists(mp)).eqs_assert(identical(.eqs_load(mp)$spec,spec),'Inputs/configuration/code changed. Use a new output_dir; caches are never mixed.')
 else .eqs_save(list(id=id,spec=spec),mp)
 for(folder in c('inputs','simulations','sample_fits','analysis','code'))dir.create(file.path(root,folder),showWarnings=FALSE)
 for(p in c(.eq3_file,file.path(dirname(.eq3_file),'eqtl_helpers.R'))){
  dest<-file.path(root,'code',basename(p))
  if(normalizePath(p)!=normalizePath(dest,mustWork=FALSE)).eqs_assert(file.copy(p,dest,overwrite=TRUE),'Cannot preserve source code')
 }
 if(!file.exists(file.path(root,'inputs','reference.rds'))).eqs_save(e$.da_load(sourcepaths['reference']),file.path(root,'inputs','reference.rds'))
 .eqs_save(base,file.path(root,'inputs','baseline_scores.rds'))
 .eqs_csv(data.frame(gene=excluded),file.path(root,'excluded_previous_genes.csv'))
 dp<-file.path(root,'inputs','design.rds')
 if(file.exists(dp))design<-.eqs_load(dp)else{
  design<-.eq3_design(inp,base,target_types,targets_per_quintile,shared_per_pair,offtype_per_type,null_genes,
   expressed_fraction,max_source_logCPM_sd,max_source_panel_share,excluded,seed,root);.eqs_save(design,dp)
 }
 .eqs_csv(design$design,file.path(root,'nominated_variant_gene_design.csv'))
 .eqs_csv(design$effects,file.path(root,'injected_gene_type_design.csv'))
 coverage<-merge(design$effects,design$design,by=c('variant','gene'),sort=FALSE)
 .eqs_csv(as.data.frame(table(coverage$pattern,coverage$cell_type)),file.path(root,'target_pattern_counts.csv'))
 jobs<-bulks<-scenarios<-qc<-cal<-preflight<-list();si<-0L
 for(mi in seq_along(mafs))for(rep in seq_len(replicates))for(strength in signal_ratios){
  si<-si+1L;name<-sprintf('scenario_%03d',si);simseed<-as.integer(seed+mi*1000L+rep)
  path<-file.path(root,'simulations',paste0(name,'.rds'))
  if(file.exists(path))sim<-.eqs_load(path)else{
   message('Preparing ',name,'; MAF ',mafs[mi],'; signal ',strength,'; replicate ',rep)
   sim<-.eq3_make(inp,design$design,design$effects,mafs[mi],strength,rep,simseed,max_abs_beta);.eqs_save(sim,path)
  }
  .eqs_assert(identical(sim$design,design$design)&&identical(sim$effects,design$effects)&&
   identical(sim$seed,simseed)&&sim$strength==strength&&sim$maf==mafs[mi]&&sim$replicate==rep,'Cached simulation context differs')
  scenarios[[si]]<-data.frame(scenario=name,maf=mafs[mi],signal_ratio=strength,replicate=rep,seed=simseed)
  bulks[[name]]<-sim$bulk;qc[[si]]<-cbind(scenario=name,sim$QC);cal[[si]]<-cbind(scenario=name,sim$calibration)
  .eqs_csv(sim$genotype_QC,file.path(root,'simulations',paste0(name,'_genotype_QC.csv')))
  for(engine in c('expected','observed')){
   tt<-.eq3_tests(if(engine=='expected')sim$expected else sim$truth,sim)
   ix<-match(.eqs_key(tt),.eqs_key(sim$effects));inj<-strength>0&!is.na(ix)
   call<-is.finite(tt$q)&tt$q<=.05;sg<-rep(0,nrow(tt));sg[inj]<-sim$effects$effect_sign[ix[inj]]
   for(ct in types){ii<-which(tt$cell_type==ct)
    preflight[[length(preflight)+1L]]<-data.frame(scenario=name,engine=engine,cell_type=ct,signal_ratio=strength,
     n_injected=sum(inj[ii]),correct_detected=sum(call[ii]&inj[ii]&sign(tt$effect[ii])==sg[ii],na.rm=TRUE),
     false_calls=sum(call[ii]&!inj[ii]))}
  }
  for(i in seq_len(ncol(sim$bulk)))jobs[[length(jobs)+1L]]<-list(cohort=name,i=i,seed=as.integer(seed+200000L+si*1000L+i))
 }
 manifest<-do.call(rbind,scenarios)
 .eqs_csv(manifest,file.path(root,'scenarios.csv'));.eqs_csv(do.call(rbind,qc),file.path(root,'simulation_QC.csv'))
 .eqs_csv(do.call(rbind,cal),file.path(root,'signal_calibration.csv'))
 .eqs_csv(do.call(rbind,preflight),file.path(root,'component_preflight_power.csv'))
 .eqs_csv(data.frame(sample=dimnames(inp$validation$truth)[[3]],donor=inp$validation$donors),file.path(root,'donors_PRIVATE.csv'))
 .eqs_csv(data.frame(check=c('training_donors','validation_donors','reference_donors','training_validation_overlap',
  'training_reference_overlap','validation_reference_overlap'),value=c(length(inp$training$donors),length(inp$validation$donors),
  length(inp$reference_donors),0,0,0)),file.path(root,'donor_checks.csv'))
 message('VERIFIED: ',nrow(design$design),' nominated genes; ',nrow(design$effects),' affected gene/type pairs; ',
  nrow(manifest),' scenarios; ',length(jobs),' donor fits. No scorer fitting.')
 writeLines(c(paste('Nominated genes:',nrow(design$design)),paste('Injected pairs in nonzero scenarios:',nrow(design$effects)),
  paste('Planned CDSeq2 donor fits:',length(jobs)),'Primary ranking: post-injection GAM6, matched-tested top 40%, Excit.',
  'Simulations verified before any fits. Input donors are previously examined; this is exploratory.',
  'Observed source components construct mixtures; estimated CDSeq2 RNA proportions enter predictors only.'),file.path(root,'PREPARED.txt'))
 if(action%in%c('check','prepare'))return(invisible(root))
 state<-list(root=root,id=id,fit=inp$fit,runtime=runtime,script_md5=.eqs_md5(workflow),types=types,bulks=bulks)
 if(action!='analyze'){
  name<-manifest$scenario[which(manifest$signal_ratio==max(manifest$signal_ratio))[1]]
  biggest<-which.max(colSums(bulks[[name]]))
  probe<-which(vapply(jobs,function(j)j$cohort==name&&j$i==biggest,TRUE))
  e$.da_run_fits(root,state,jobs[probe],1L,workflow)
  if(action=='probe')return(invisible(root))
  e$.da_run_fits(root,state,jobs,as.integer(workers),workflow)
 }
 all<-list()
 for(name in manifest$scenario){
  sim<-.eqs_load(file.path(root,'simulations',paste0(name,'.rds')))
  est<-e$.da_assemble(root,name,sim,id)
  .eqs_save(est,file.path(root,'analysis',paste0(name,'_estimates.rds')))
  all[[name]]<-.eq3_evaluate(sim,est,base,gam,inp,e,name,root)
 }
 d<-do.call(rbind,all);rownames(d)<-NULL
 .eqs_validate_details(d,.eq3_methods)
 .eqs_save(d,file.path(root,'all_test_details.rds'));.eqs_csv(d,file.path(root,'all_test_details.csv'))
 r<-.eqs_report(d,.eq3_methods,root,patterns=TRUE)
 .eq3_attribution(d,root);.eq3_plot(r$pooled,target_types,root)
 writeLines(c('Completed three-type nominated-variant eQTL simulation v1.0.0.',
  'No pretrained models supplied; this experiment explicitly reuses the user-trained application GAM6.',
  'Primary: post-injection/scenario GAM6 top40 among nominated tests, Excit. Other types/patterns and full-panel quintiles are descriptive.',
  'All seven mixture types remain in each global BH family; untested reserve genes are NOT part of the nominated testing family.',
  'Signals calibrate the increment in expected logCPM genotype slope times genotype SD / baseline expression SD.',
  'The requested ratio is not a heritability and does not force equal oracle p-values or total slopes.',
  'Shared types use the same variant/genotype vector; opposite means opposite injected slope increments.',
  'Nominated non-target profiles remain unchanged in expectation; reserve genes absorb count budgets.',
  'Target eligibility limits baseline variability and per-donor panel dominance before genotypes are drawn; inspect candidate_support.csv.',
  'No-injection components are also multinomial draws; they are not identical to direct observed aggregation.',
  'Scoring uses reference means, observed simulated bulk and estimated RNA fractions; never genotype or truth.',
  'Both expected and realized component associations are evaluated separately.',
  'New targets do not create new biological validation; all repetitions reuse donors and nominated genes.',
  'FDP is descriptive, not an FDR guarantee. Score filtering cannot increase total discoveries at fixed original q.',
  'This nominated-SNP simulation does not evaluate LD-aware genome-wide cis-eQTL discovery or fine mapping.'),file.path(root,'SUMMARY.txt'))
 capture.output(sessionInfo(),file=file.path(root,'sessionInfo.txt'))
 writeLines(as.character(Sys.time()),file.path(root,'COMPLETED'));message('COMPLETED: ',root)
 invisible(root)
}
