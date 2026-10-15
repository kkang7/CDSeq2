# Reconstruct Figure 3 after excluding archived evaluation donor index 3.
# Original simulated values/40-donor phenotype scale and frozen GAM6 ranks are retained.
args<-commandArgs(TRUE);p<-args[1];o<-args[2];td<-file.path(o,'supplementary_tables');sources<-list()
load<-function(rel){f<-file.path(p,rel);h<-unname(tools::md5sum(f));if(file.exists(paste0(f,'.md5')))stopifnot(h==strsplit(trimws(readLines(paste0(f,'.md5'),warn=FALSE)[1]),'[[:space:]]+')[[1]][1]);sources[[length(sources)+1]]<<-data.frame(file=rel,md5=h);if(grepl('csv$',f))read.csv(f)else readRDS(f)}
ph<-load('continuous_associations_v1/synthetic_phenotypes.rds');ids<-names(ph)[vapply(ph,function(x)x$source_type=='Excit'&&x$program=='Excit_program1'&&x$signal_fraction==.5&&x$replicate==1,logical(1))];stopifnot(length(ids)==1);trait<-ph[[ids]]
truth<-load('inputs/direct_inputs.rds')$validation$truth;est<-load('analysis/validation_estimates.rds')$counts;stopifnot(identical(dimnames(truth),dimnames(est)),identical(dimnames(est)[[3]],names(trait$phenotype)))
lognorm<-function(a)log2(sweep(a,c(2,3),pmax(apply(a,c(2,3),sum),1),'/')*1e6+.5)
keep<-setdiff(seq_len(dim(truth)[3]),3L)
x<-lognorm(truth)[,,keep,drop=FALSE];y<-lognorm(est)[,,keep,drop=FALSE];phen<-as.numeric(scale(trait$phenotype))[keep];phc<-phen-mean(phen)
rank<-load('continuous_associations_v1/frozen_rankings.csv');r<-subset(rank,cell_type=='Excit');r<-r[order(-r$GAM6,r$gene),];selected<-r$gene[1:400]
program<-load('continuous_associations_v1/program_genes.csv');genes<-subset(program,program==trait$program&used_in_PC1)$gene
# Choose by the frozen phenotype-independent score, not association fit.
gene<-r$gene[which(r$gene %in% genes)[1]];cat('Display trait:',ids,'gene:',gene,'score:',r$GAM6[match(gene,r$gene)],'\n')
slope<-function(m)as.vector(m%*%phc/sum(phc^2));sx<-slope(x[,'Excit',]);sy<-slope(y[,'Excit',]);z<-data.frame(trait=ids,gene=dimnames(x)[[1]],cell_type='Excit',reference_slope=sx,estimated_slope=sy,selected=dimnames(x)[[1]]%in%selected)
write.csv(z,file.path(td,'S16_continuous_slope_recovery.csv'),row.names=FALSE)
d<-data.frame(trait=ids,gene,cell_type='Excit',donor_index=keep,phenotype_SD=phen,reference_expression=x[gene,'Excit',],estimated_expression=y[gene,'Excit',]);write.csv(d,file.path(td,'S17_continuous_donor_example.csv'),row.names=FALSE)
write.csv(do.call(rbind,sources),file.path(td,'application_panels_source_manifest.csv'),row.names=FALSE)
cat('Slope concordance all/selected:',cor(sx,sy),cor(sx[z$selected],sy[z$selected]),'\n')
