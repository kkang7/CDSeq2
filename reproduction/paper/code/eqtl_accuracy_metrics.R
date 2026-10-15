# Usage: Rscript --vanilla eqtl_accuracy_metrics.R OUTPUT_FOLDER
# Uses saved associations only. No deconvolution, scoring or BH adjustment is rerun.
# Direction-aware definitions: TP = injected target with correct inferred sign;
# FP = every other reported call (non-injected or wrong-sign); FN = all injected - TP.
args<-commandArgs(TRUE);o<-if(length(args))args[1]else'.';td<-file.path(o,'supplementary_tables')
t<-read.csv(file.path(td,'S8_eqtl_tests.csv'));t<-subset(t,engine%in%c('CDSeq2','BayesPrism'))
cts<-c('Astro','Excit','Inhib','Oligo','OPC','Immune','Vascular');wf<-c('CDSeq2 before','CDSeq2 GAM6','BayesPrism standard')
# Preserve the original matched-tested-budget eligibility and gene tie breaking.
t$GAM6_selected<-FALSE
for(sc in unique(t$scenario))for(ct in cts){ix<-which(t$scenario==sc&t$cell_type==ct&t$engine=='CDSeq2');z<-t[ix,];eligible<-rowSums(!is.finite(as.matrix(z[c('scenario_GAM6_score','scenario_contribution_score','scenario_expected_count_score')])))==0L
 j<-which(eligible);j<-j[order(-z$scenario_GAM6_score[j],z$gene[j],method='radix')];j<-head(j,floor(.4*length(j)+1e-9));t$GAM6_selected[ix[j]]<-TRUE}
et<-read.csv(file.path(td,'S18_eqtl_before_after.csv'));key<-function(d)paste(d$scenario,d$engine,d$gene,d$cell_type,sep=':')
ix<-match(key(subset(t,engine=='CDSeq2'&cell_type=='Excit')),key(et));stopifnot(identical(t$GAM6_selected[t$engine=='CDSeq2'&t$cell_type=='Excit'],et$selected[ix]))
metrics<-function(ncall,tp,fp_null,fp_sign,P,threshold){fp<-fp_null+fp_sign;stopifnot(ncall==tp+fp)
 data.frame(q_threshold=threshold,calls=ncall,TP=tp,FP=fp,FP_noninjected=fp_null,FP_wrong_sign=fp_sign,FN=P-tp,n_injected_all=P,precision=if(ncall)tp/ncall else NA_real_,recall=if(P)tp/P else NA_real_,F1=if(P)2*tp/(ncall+P)else NA_real_,directional_FDP=if(ncall)fp/ncall else 0,noninjected_FDP=if(ncall)fp_null/ncall else 0)}
curve<-function(z,P){z<-z[is.finite(z$q),];z<-z[order(z$q),];correct<-z$injected&is.finite(z$effect)&sign(z$effect)==z$injected_sign;correct[is.na(correct)]<-FALSE
 wrong<-z$injected&!correct;null<-!z$injected;end<-which(!duplicated(z$q,fromLast=TRUE));ct<-cumsum(correct);cw<-cumsum(wrong);cn<-cumsum(null)
 empty<-metrics(0,0,0,0,P,NA_real_)
 if(!length(end))return(empty)
 rows<-data.frame(q_threshold=z$q[end],calls=end,TP=ct[end],FP=cn[end]+cw[end],FP_noninjected=cn[end],FP_wrong_sign=cw[end],FN=P-ct[end],n_injected_all=P,precision=ct[end]/end,recall=ct[end]/P,F1=2*ct[end]/(end+P),directional_FDP=(cn[end]+cw[end])/end,noninjected_FDP=cn[end]/end)
 rbind(empty,rows)
}
curves<-fixed<-matched<-list();n<-0L
caps<-sort(unique(c(seq(0,.5,.01),1)))
for(h in c(.5,.8))for(rep in c(0,1:3))for(ct in c(cts,'ALL_TYPES'))for(w in wf){
 eng<-if(w=='BayesPrism standard')'BayesPrism'else'CDSeq2';z<-subset(t,signal_ratio==h&engine==eng);if(rep>0)z<-subset(z,replicate==rep);if(ct!='ALL_TYPES')z<-subset(z,cell_type==ct)
 P<-sum(z$injected);nuniverse<-nrow(z);available<-sum(is.finite(z$q));if(w=='CDSeq2 GAM6')z<-subset(z,GAM6_selected)
 tag<-data.frame(workflow=w,engine=eng,signal_ratio=h,replicate=rep,cell_type=ct,n_tests_universe=nuniverse,n_available_before_selection=available,n_retained=nrow(z),n_available_retained=sum(is.finite(z$q)))
 c<-curve(z,P);n<-n+1L;curves[[n]]<-cbind(tag,c)
 called<-is.finite(z$q)&z$q<=.05;correct<-called&z$injected&is.finite(z$effect)&sign(z$effect)==z$injected_sign;correct[is.na(correct)]<-FALSE
 fixed[[n]]<-cbind(tag,metrics(sum(called),sum(correct),sum(called&!z$injected),sum(called&z$injected&!correct),P,.05))
 for(cap in caps){ok<-c[c$directional_FDP<=cap+1e-12,];# No randomization or splitting tied q values.
  # Maximize true discoveries under the common cap; among ties choose fewer
  # incorrect calls, then the strictest threshold. The no-call row is feasible.
  ix<-order(-ok$TP,ok$FP,replace(ok$q_threshold,is.na(ok$q_threshold),-Inf))[1];r<-ok[ix,,drop=FALSE]
  matched[[length(matched)+1L]]<-cbind(tag,FDP_cap=cap,r)
 }
}
curves<-do.call(rbind,curves);fixed<-do.call(rbind,fixed);matched<-do.call(rbind,matched)
write.csv(curves,file.path(td,'S20_precision_recall_curves.csv'),row.names=FALSE)
write.csv(fixed,file.path(td,'S21_precision_recall_F1.csv'),row.names=FALSE)
write.csv(matched,file.path(td,'S22_matched_empirical_FDP.csv'),row.names=FALSE)
# Check the q<=0.05 summaries against the original CDSeq2 top-40% summaries.
old<-read.csv(file.path(td,'S5_eqtl_summary.csv'));old<-subset(old,engine=='CDSeq2'&score_scope=='scenario'&method=='GAM6'&policy=='matched_tested_budget'&group=='top40'&pattern=='all'&signal_ratio>0)
f<-subset(fixed,workflow=='CDSeq2 GAM6'&replicate==0&cell_type!='ALL_TYPES');j<-match(paste(f$cell_type,f$signal_ratio),paste(old$cell_type,old$signal_ratio));stopifnot(!anyNA(j),all(f$calls==old$calls[j]),all(f$TP==old$correct_detected[j]),all(f$FP_noninjected==old$false_calls[j]),all(f$FP_wrong_sign==old$wrong_sign_calls[j]))
print(subset(fixed,cell_type=='Excit'&replicate==0,c(workflow,signal_ratio,TP,FP,FN,precision,recall,F1)))
print(subset(matched,cell_type=='Excit'&replicate==0&FDP_cap%in%c(.05,.1,.2),c(workflow,signal_ratio,FDP_cap,q_threshold,TP,FP,directional_FDP,recall)))
