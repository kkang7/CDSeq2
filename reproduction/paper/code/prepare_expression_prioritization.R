# Add optional GAM6 selection to the existing common-panel expression benchmark.
# Use EXPERIMENT_ROOT to read archived scores, or '-' to reproduce from Table S23.
args <- commandArgs(TRUE)
if(length(args)!=2L)stop('Usage: Rscript prepare_expression_prioritization.R EXPERIMENT_ROOT_OR_DASH OUTPUT_FOLDER')
root<-args[1]; o<-args[2]; td<-file.path(o,'supplementary_tables')
cts<-c('Astro','Excit','Inhib','Oligo','OPC','Immune','Vascular')
p<-read.csv(file.path(td,'S14_expression_accuracy.csv'))
p<-subset(p,scenario!='direct' & oracle=='observed')
key<-function(d)paste(d$scenario,d$cell_type,d$gene,sep='|')
cd<-subset(p,engine=='CDSeq2'); bp<-subset(p,engine=='BayesPrism')
stopifnot(!anyDuplicated(key(cd)),!anyDuplicated(key(bp)),setequal(key(cd),key(bp)))
bp<-bp[match(key(cd),key(bp)),]
score_file<-file.path(td,'S23_S6_GAM6_membership.csv')
manifest<-list()
if(root=='-'){
 scores<-read.csv(score_file)
}else{
 scores<-do.call(rbind,lapply(sort(unique(cd$scenario)),function(sc){
  f<-file.path(root,'eqtl_three_types_v1','analysis',paste0(sc,'_scores.rds'))
  md5<-unname(tools::md5sum(f)); expected<-trimws(readLines(paste0(f,'.md5'),warn=FALSE)[1])
  stopifnot(identical(md5,expected))
  manifest[[sc]]<<-data.frame(source=f,md5=md5)
  d<-readRDS(f);d$scenario<-sc
  d[c('scenario','cell_type','gene','GAM6','common_rankable','GAM6_outside_training_range')]
 }))
}
stopifnot(!anyDuplicated(key(scores)),all(key(cd)%in%key(scores)))
m<-scores[match(key(cd),key(scores)),c('scenario','cell_type','gene','GAM6','common_rankable','GAM6_outside_training_range')]
m$rank_eligible<-is.finite(m$GAM6)
m$GAM6_rank<-NA_integer_;m$n_common_panel<-NA_integer_;m$n_score_eligible<-NA_integer_;m$selected<-FALSE
# Selection uses ONLY finite archived score, scenario, cell type and gene-name ties.
# Undefined observed correlations and out-of-range flags do not alter membership.
for(sc in sort(unique(m$scenario)))for(ct in cts){
 ix<-which(m$scenario==sc&m$cell_type==ct);el<-ix[m$rank_eligible[ix]]
 el<-el[order(-m$GAM6[el],m$gene[el],method='radix')]
 m$GAM6_rank[el]<-seq_along(el);m$n_common_panel[ix]<-length(ix);m$n_score_eligible[ix]<-length(el)
 m$selected[head(el,floor(.2*length(el)))]<-TRUE
}
# Independently match archived scores to existing nominated-test records.
t<-read.csv(file.path(td,'S8_eqtl_tests.csv'));t<-subset(t,engine=='CDSeq2')
ix<-match(key(t),key(m));ok<-!is.na(ix)
stopifnot(all(abs(t$scenario_GAM6_score[ok]-m$GAM6[ix[ok]])<1e-9))
write.csv(m,score_file,row.names=FALSE)
rows<-list()
for(sc in c(sort(unique(m$scenario)),'ALL'))for(ct in cts)for(w in c('BayesPrism','CDSeq2','CDSeq2+GAM6')){
 ix<-which(m$cell_type==ct & (sc=='ALL'|m$scenario==sc))
 take<-if(w=='CDSeq2+GAM6')m$selected[ix]else rep(TRUE,length(ix))
 vals<-if(w=='BayesPrism')bp$Pearson[ix]else cd$Pearson[ix]
 v<-vals[take];finite<-is.finite(v)
 rows[[length(rows)+1L]]<-data.frame(scenario=sc,cell_type=ct,workflow=w,
  n_common_panel=length(ix),n_GAM6_eligible=sum(m$rank_eligible[ix]),
  n_retained=sum(take),retained_proportion=sum(take)/length(ix),
  n_defined=sum(finite),n_undefined=sum(!finite),
  median_Pearson=median(v[finite]),q25_Pearson=unname(quantile(v[finite],.25)),q75_Pearson=unname(quantile(v[finite],.75)))
}
res<-do.call(rbind,rows)
write.csv(res,file.path(td,'S24_S6_recovery_counts.csv'),row.names=FALSE)
if(root!='-'){
 src<-file.path(td,c('S14_expression_accuracy.csv','S8_eqtl_tests.csv'))
 provenance<-rbind(do.call(rbind,manifest),data.frame(source=src,md5=unname(tools::md5sum(src))))
 write.csv(provenance,file.path(td,'expression_selection_source_manifest.csv'),row.names=FALSE)
}
print(subset(res,scenario=='ALL'),row.names=FALSE)
cat('All common-panel scores finite:',all(m$rank_eligible),'\n')
cat('All common-panel profiles common_rankable:',all(m$common_rankable),'\n')
