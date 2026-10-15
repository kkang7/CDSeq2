# Original 40-donor and all leave-one-out diagnostics. Figure 3 now shows
# the donor-index-3 omission; inference and GAM6 membership remain fixed.
args<-commandArgs(TRUE)
if(length(args)!=2L)stop('Usage: Rscript continuous_slope_sensitivity.R EXPERIMENT_ROOT OUTPUT_FOLDER')
p<-args[1];td<-file.path(args[2],'supplementary_tables');sources<-list()
load<-function(rel){
 f<-file.path(p,rel);h<-unname(tools::md5sum(f))
 if(file.exists(paste0(f,'.md5')))stopifnot(h==strsplit(trimws(readLines(paste0(f,'.md5'),warn=FALSE)[1]),'[[:space:]]+')[[1]][1])
 sources[[length(sources)+1L]]<<-data.frame(source=rel,md5=h)
 readRDS(f)
}
a<-load('inputs/direct_inputs.rds')$validation$truth
b<-load('analysis/validation_estimates.rds')$counts
trait<-load('continuous_associations_v1/synthetic_phenotypes.rds')[['trait_031']]
stopifnot(identical(dimnames(a),dimnames(b)),identical(dimnames(a)[[3]],names(trait$phenotype)),
 trait$source_type=='Excit',trait$program=='Excit_program1',trait$signal_fraction==.5,trait$replicate==1)
norm<-function(a)log2(sweep(a,c(2,3),pmax(apply(a,c(2,3),sum),1),'/')*1e6+.5)
x<-norm(a)[,'Excit',];y<-norm(b)[,'Excit',];phen<-as.numeric(scale(trait$phenotype))
s<-read.csv(file.path(td,'S16_continuous_slope_recovery.csv'));s<-s[match(rownames(x),s$gene),]
d<-read.csv(file.path(td,'S17_continuous_donor_example.csv'))
stopifnot(!anyNA(s$gene),all(s$trait=='trait_031'),sum(s$selected)==400,all(abs(phen[d$donor_index]-d$phenotype_SD)<1e-10),identical(d$donor_index,setdiff(1:40,3L)))
rows<-list()
for(exclude in 0:40){
 keep<-setdiff(1:40,exclude);z<-phen[keep]-mean(phen[keep])
 sx<-as.vector(x[,keep]%*%z/sum(z^2));sy<-as.vector(y[,keep]%*%z/sum(z^2))
 if(exclude==3){
  stopifnot(max(abs(sx-s$reference_slope))<1e-10,max(abs(sy-s$estimated_slope))<1e-10)
  # Check the matrix formula against independent lm fits for representative genes.
  for(g in c(1,200,800,1400,2000))stopifnot(abs(coef(lm(y[g,keep]~phen[keep]))[2]-sy[g])<1e-10)
 }
 rows[[length(rows)+1L]]<-data.frame(trait='trait_031',source_type='Excit',program='Excit_program1',signal_fraction=.5,noise_repetition=1,
  omitted_donor_index=exclude,n_donors=length(keep),omitted_phenotype_SD=if(exclude==0)NA else phen[exclude],
  n_genes=length(sx),n_selected=sum(s$selected),
  slope_concordance_all=cor(sx,sy),slope_concordance_selected=cor(sx[s$selected],sy[s$selected]),
  median_absolute_slope_error_all=median(abs(sy-sx)),median_absolute_slope_error_selected=median(abs(sy[s$selected]-sx[s$selected])))
}
res<-do.call(rbind,rows)
write.csv(res,file.path(td,'S25_continuous_donor_sensitivity.csv'),row.names=FALSE)
f<-file.path(td,c('S16_continuous_slope_recovery.csv','S17_continuous_donor_example.csv'))
write.csv(rbind(do.call(rbind,sources),data.frame(source=basename(f),md5=unname(tools::md5sum(f)))),file.path(td,'continuous_sensitivity_source_manifest.csv'),row.names=FALSE)
print(res[res$omitted_donor_index %in% c(0,which.min(phen)),],row.names=FALSE)
cat('Figure 3 uses the donor-3 omission; the complete 40-donor diagnostic is retained here.\n')
