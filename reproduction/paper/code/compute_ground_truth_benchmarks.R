# Usage: Rscript compute_ground_truth_benchmarks.R EXPERIMENT_ROOT OUTPUT_FOLDER NATIVE_CACHE
args<-commandArgs(TRUE);p<-args[1];o<-args[2];cache<-args[3];td<-file.path(o,'supplementary_tables')
inputs_used<-list()
load_verified<-function(path){
 hash<-unname(tools::md5sum(path)); sidecar<-paste0(path,'.md5')
 if(file.exists(sidecar)){expected<-strsplit(trimws(readLines(sidecar,warn=FALSE)[1]),'[[:space:]]+')[[1]][1];stopifnot(tolower(hash)==tolower(expected))}
 inputs_used[[length(inputs_used)+1]]<<-data.frame(source=sub(paste0(p,'/'),'',path,fixed=TRUE),md5=hash)
 readRDS(path)
}
norm<-function(a)sweep(a,c(2,3),pmax(apply(a,c(2,3),sum),1),'/')
frac<-function(a){x<-t(apply(a,c(2,3),sum));x/rowSums(x)}
prof<-list();fr<-list();summ<-list()
bench<-function(est,truth,scenario,engine,oracle){
 stopifnot(identical(dimnames(est$counts),dimnames(truth)))
 a<-norm(est$counts);b<-norm(truth);tf<-frac(truth);ef<-est$proportions[rownames(tf),colnames(tf),drop=FALSE]
 for(k in seq_len(dim(a)[2])){
 ct<-dimnames(a)[[2]][k];r<-vapply(seq_len(dim(a)[1]),function(g)suppressWarnings(cor(a[g,k,],b[g,k,])),numeric(1));rm<-sqrt(rowMeans((a[,k,]-b[,k,])^2))
 prof[[length(prof)+1]]<<-data.frame(scenario,engine,oracle,cell_type=ct,gene=dimnames(a)[[1]],Pearson=r,RMSE_CPM=rm*1e6)
 fr[[length(fr)+1]]<<-data.frame(scenario,engine,oracle,cell_type=ct,donor_index=seq_len(nrow(tf)),truth=tf[,k],estimate=ef[,k])
 summ[[length(summ)+1]]<<-data.frame(scenario,engine,oracle,cell_type=ct,fraction_r=cor(tf[,k],ef[,k]),fraction_MAE=mean(abs(tf[,k]-ef[,k])),fraction_RMSE=sqrt(mean((tf[,k]-ef[,k])^2)),expression_median_r=median(r,na.rm=TRUE),expression_defined=sum(is.finite(r)),expression_median_RMSE_CPM=median(rm*1e6))
 }
}
x<-load_verified(file.path(p,'inputs/direct_inputs.rds'));e<-load_verified(file.path(p,'analysis/validation_estimates.rds'));bench(e,x$validation$truth,'direct','CDSeq2','observed')
for(i in 1:9){s<-sprintf('scenario_%03d',i);sim<-load_verified(file.path(p,'eqtl_three_types_v1/simulations',paste0(s,'.rds')));cd<-load_verified(file.path(p,'eqtl_three_types_v1/analysis',paste0(s,'_estimates.rds')));bp<-load_verified(file.path(cache,paste0(s,'.rds')));genes<-dimnames(bp$counts)[[1]];cd_panel<-cd;mass<-t(apply(cd$counts[genes,,,drop=FALSE],c(2,3),sum)/apply(cd$counts,c(2,3),sum));cd_panel$proportions<-cd$proportions*mass[rownames(cd$proportions),colnames(cd$proportions)];cd_panel$proportions<-cd_panel$proportions/rowSums(cd_panel$proportions);cd_panel$counts<-cd$counts[genes,,,drop=FALSE];for(or in c('expected','observed')){truth<-(if(or=='expected')sim$expected else sim$truth)[genes,,,drop=FALSE];bench(cd_panel,truth,s,'CDSeq2',or);bench(bp,truth,s,'BayesPrism',or)};cat(s,'complete\n')}
write.csv(do.call(rbind,summ),file.path(td,'S12_ground_truth_accuracy.csv'),row.names=FALSE)
write.csv(do.call(rbind,fr),file.path(td,'S13_RNA_fraction_pairs.csv'),row.names=FALSE)
write.csv(do.call(rbind,prof),file.path(td,'S14_expression_accuracy.csv'),row.names=FALSE)
ss<-do.call(rbind,summ);print(aggregate(cbind(fraction_r,fraction_MAE,expression_median_r)~engine+cell_type,subset(ss,scenario!='direct'&oracle=='observed'),mean))

write.csv(do.call(rbind,inputs_used),file.path(td,'benchmark_source_manifest.csv'),row.names=FALSE)
