# Usage: Rscript verify_standard_comparison.R EXPERIMENT_ROOT NATIVE_CACHE OUTPUT_FOLDER
# Independent checks of native preprocessing, per-test statistics and CDSeq2 preservation.
a<-commandArgs(TRUE);p<-a[1];cache<-a[2];o<-a[3];td<-file.path(o,'supplementary_tables')
t<-read.csv(file.path(td,'S8_eqtl_tests.csv'));et<-read.csv(file.path(td,'S18_eqtl_before_after.csv'));sm<-read.csv(file.path(td,'S6_standard_comparison.csv'))
old<-read.csv(file.path(p,'eqtl_bayesprism_comparison_v1/comparison_test_details.csv'))
k<-function(d)paste(d$scenario,d$engine,d$gene,d$cell_type,sep=':')
cd<-subset(t,engine=='CDSeq2');ix<-match(k(cd),k(old));stopifnot(!anyNA(ix))
for(n in c('effect','p','q','called','correct_detection'))stopifnot(isTRUE(all.equal(cd[[n]],old[[n]][ix],tolerance=1e-12)))
stopifnot(!any(grepl('^BP_',names(t))))
qc<-read.csv(file.path(td,'S19_BayesPrism_preprocessing.csv'))
for(i in 1:9){sc<-sprintf('scenario_%03d',i);fit<-readRDS(file.path(cache,paste0(sc,'.rds')));sim<-readRDS(file.path(p,'eqtl_three_types_v1/simulations',paste0(sc,'.rds')));z<-subset(t,scenario==sc&engine=='BayesPrism');q<-p.adjust(replace(z$p,!is.finite(z$p),1),'BH');q[!is.finite(z$p)]<-NA_real_;stopifnot(nrow(z)==1778,identical(q,z$q)||isTRUE(all.equal(q,z$q,tolerance=1e-12)),sum(is.finite(z$p))<=sum(z$gene%in%dimnames(fit$counts)[[1]]),fit$conservation_error<1e-6)
 genes<-dimnames(fit$counts)[[1]];actual_error<-max(abs(apply(fit$counts,c(1,3),sum)-sim$bulk[genes,,drop=FALSE])/pmax(sim$bulk[genes,,drop=FALSE],1));stopifnot(actual_error<1e-6,all(is.finite(fit$counts)),all(fit$counts>=0),max(abs(rowSums(fit$proportions)-1))<1e-6,fit$controls$gibbs.control$seed==20261013L+i,isTRUE(fit$controls$update.gibbs))
 qcz<-subset(qc,scenario==sc);stopifnot(setequal(qcz$gene[qcz$retained],dimnames(fit$counts)[[1]]),sum(qcz$cleanup_retained)==1981,sum(qcz$retained)==length(genes),all(!is.finite(z$p[!z$gene%in%genes])))
 stopifnot(all(z$availability_status[is.finite(z$p)]=='tested'))
 y<-log2(sweep(fit$counts,c(2,3),apply(fit$counts,c(2,3),sum),'/')*1e6+.5)
 # lm supplies an independent regression calculation for every available test.
 for(j in which(is.finite(z$p))){yy<-y[z$gene[j],z$cell_type[j],];gg<-sim$genotypes[,z$variant[j]];cf<-summary(lm(yy~gg))$coefficients;stopifnot(abs(cf[2,1]-z$effect[j])<1e-10,abs(cf[2,4]-z$p[j])<1e-10)}
 s<-subset(et,scenario==sc&engine=='CDSeq2');stopifnot(sum(s$selected)==101)
 rank<-order(-s$scenario_GAM6_score,s$gene);stopifnot(setequal(s$gene[s$selected],s$gene[head(rank,101)]))
 stopifnot(all(subset(et,scenario==sc&engine=='BayesPrism')$selected))
 cat(sc,'native tests and unchanged CDSeq2 selection verified\n')
}
for(j in seq_len(nrow(sm))){r<-sm[j,];z<-subset(et,signal_ratio==r$signal_ratio&engine==r$engine);if(r$workflow=='CDSeq2 GAM6')z<-subset(z,selected);stopifnot(sum(z$called)==r$calls,sum(z$correct_detection)==r$correct_detected,sum(z$called&!z$injected)==r$false_calls,sum(z$wrong_sign_call)==r$wrong_sign_calls)}
# Outcome counts are independently recomputed above; no historical result is required.
cat('All comparison checks passed.\n')
