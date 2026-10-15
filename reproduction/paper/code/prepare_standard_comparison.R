# Usage: Rscript prepare_standard_comparison.R EXPERIMENT_ROOT NATIVE_CACHE OUTPUT_FOLDER
# Recomputes native BayesPrism tests; CDSeq2/oracle p/q and CDSeq2 selection stay unchanged.
a<-commandArgs(TRUE);p<-a[1];cache<-a[2];o<-a[3];td<-file.path(o,'supplementary_tables')
base<-read.csv(file.path(p,'eqtl_bayesprism_comparison_v1/comparison_test_details.csv'))
# Remove all archived BP support scores from delivered per-test data.
base<-base[,!grepl('^BP_',names(base))];base<-subset(base,engine!='BayesPrism');base$availability_status<-'tested';base$preprocessing<-'original 2000-gene panel'
rows<-qc<-src<-runmeta<-list()
for(i in 1:9){sc<-sprintf('scenario_%03d',i);file<-file.path(cache,paste0(sc,'.rds'));fit<-readRDS(file);sim<-readRDS(file.path(p,'eqtl_three_types_v1/simulations',paste0(sc,'.rds')))
 stopifnot(fit$sim_md5==unname(tools::md5sum(file.path(p,'eqtl_three_types_v1/simulations',paste0(sc,'.rds')))))
 runmeta[[i]]<-data.frame(scenario=sc,retained_genes=dim(fit$counts)[1],retained_nominated=sum(fit$preprocessing$retained&fit$preprocessing$nominated),conservation_error=fit$conservation_error,elapsed_seconds=fit$elapsed_seconds,cores=fit$controls$gibbs.control$n.cores,seed=fit$spec$seed,package=fit$spec$package)
 src[[i]]<-data.frame(scenario=sc,fit_md5=unname(tools::md5sum(file)),simulation_md5=fit$sim_md5,reference_md5=fit$reference_md5,seed=fit$spec$seed,package=fit$spec$package)
 z<-subset(base,scenario==sc&engine=='CDSeq2');z$engine<-'BayesPrism';z$effect<-z$p<-NA_real_
 x<-log2(sweep(fit$counts,c(2,3),apply(fit$counts,c(2,3),sum),'/')*1e6+.5)
 for(j in seq_len(nrow(z))){g<-z$gene[j];ct<-z$cell_type[j];if(!g%in%dimnames(x)[[1]])next
  yy<-x[g,ct,];gg<-sim$genotypes[,z$variant[j]];stopifnot(identical(dimnames(x)[[3]],rownames(sim$genotypes)))
  if(any(!is.finite(yy))||sd(yy)<1e-12)next
  xc<-gg-mean(gg);yc<-yy-mean(yy);ss<-sum(xc^2);b<-sum(xc*yc)/ss;se<-sqrt(sum((yc-b*xc)^2)/(length(gg)-2)/ss)
  z$effect[j]<-b;z$p[j]<-if(se>0)2*pt(abs(b/se),length(gg)-2,lower.tail=FALSE)else 0
 }
 z$q<-p.adjust(replace(z$p,!is.finite(z$p),1),'BH');z$q[!is.finite(z$p)]<-NA_real_
 z$called<-is.finite(z$q)&z$q<=.05;z$correct_detection<-z$called&z$injected&is.finite(z$effect)&sign(z$effect)==z$injected_sign;z$correct_detection[is.na(z$correct_detection)]<-FALSE
 z$wrong_sign_call<-z$called&z$injected&!z$correct_detection
 z$availability_status<-ifelse(!z$gene%in%dimnames(x)[[1]],'excluded_by_input_QC',ifelse(is.finite(z$p),'tested','constant_profile'));z$preprocessing<-'BayesPrism cleanup + default outlier QC'
 # CDSeq2-specific scoring fields are inapplicable to the BayesPrism output.
 z[,grepl('^(baseline|scenario)_(GAM6|contribution|expected_count)',names(z))]<-NA
 rows[[i]]<-z;qc[[i]]<-fit$preprocessing
}
t<-rbind(base,do.call(rbind,rows));stopifnot(!anyDuplicated(t[c('scenario','engine','gene','cell_type')]))
write.csv(t,file.path(td,'S8_eqtl_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,qc),file.path(td,'S19_BayesPrism_preprocessing.csv'),row.names=FALSE)
write.csv(do.call(rbind,src),file.path(td,'native_BayesPrism_source_manifest.csv'),row.names=FALSE)
write.csv(do.call(rbind,runmeta),file.path(td,'native_BayesPrism_run_summary.csv'),row.names=FALSE)
mem<-read.csv(file.path(p,'eqtl_bayesprism_comparison_v1/analysis/selection_membership.csv'))
mem<-subset(mem,cell_type=='Excit'&engine=='CDSeq2'&selection_policy=='own_support_or_GAM6'&ranking=='GAM6'&retention==.4)
key<-function(a)paste(a$scenario,a$gene,a$cell_type,sep=':')
et<-subset(t,cell_type=='Excit'&engine%in%c('CDSeq2','BayesPrism'))
et$selected<-TRUE;cd<-et$engine=='CDSeq2';et$selected[cd]<-mem$selected[match(key(et[cd,]),key(mem))];stopifnot(!anyNA(et$selected))
et$selection_rule<-ifelse(cd,'CDSeq2 GAM6 top 40%','No output selection');et$variant_index<-match(et$variant,sort(unique(et$variant)));et$display_index<-(et$replicate-1)*254+et$variant_index
write.csv(et,file.path(td,'S18_eqtl_before_after.csv'),row.names=FALSE)
sm<-list()
for(h in c(0,.5,.8))for(w in c('CDSeq2 before','CDSeq2 GAM6','BayesPrism standard')){
 eng<-if(w=='BayesPrism standard')'BayesPrism'else'CDSeq2';z<-subset(et,signal_ratio==h&engine==eng);nall<-sum(z$injected);if(w=='CDSeq2 GAM6')z<-subset(z,selected)
 sm[[length(sm)+1]]<-data.frame(workflow=w,engine=eng,signal_ratio=h,cell_type='Excit',n_tests=nrow(z),available_tests=sum(is.finite(z$p)),n_injected=sum(z$injected),n_injected_all=nall,available_injected=sum(z$injected&is.finite(z$p)),calls=sum(z$called),correct_detected=sum(z$correct_detection),false_calls=sum(z$called&!z$injected),wrong_sign_calls=sum(z$wrong_sign_call),pooled_FDP=if(sum(z$called))sum(z$called&!z$injected)/sum(z$called)else 0,recall_all_injected=if(nall)sum(z$correct_detection)/nall else NA,selection_policy=if(w=='CDSeq2 GAM6')'GAM6_top40'else'no_output_selection',n_replicates=3)
}
sm<-do.call(rbind,sm);write.csv(sm,file.path(td,'S6_standard_comparison.csv'),row.names=FALSE)
old<-file.path(td,'S6_matched_comparison.csv');if(file.exists(old))unlink(old)
print(sm)
