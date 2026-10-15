# Usage: Rscript --vanilla run_standard_bayesprism.R EXPERIMENT_ROOT CACHE_DIR SCENARIO_INDEX [CORES]
# BayesPrism 2.2.3: documented gene cleanup + default new.prism outlier QC.
# No output reliability/support ranking is applied to BayesPrism.
a<-commandArgs(TRUE);root<-normalizePath(a[1]);out<-a[2];i<-as.integer(a[3]);cores<-if(length(a)>3)as.integer(a[4])else 1L
stopifnot(i%in%1:9,cores>=1,as.character(packageVersion('BayesPrism'))=='2.2.3')
dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out);sc<-sprintf('scenario_%03d',i)
loadv<-function(rel){f<-file.path(root,rel);h<-unname(tools::md5sum(f));if(file.exists(paste0(f,'.md5')))stopifnot(h==strsplit(trimws(readLines(paste0(f,'.md5'),warn=FALSE)[1]),'[[:space:]]+')[[1]][1]);readRDS(f)}
prep<-loadv('eqtl_bayesprism_comparison_v1/prepared_reference.rds');simrel<-paste0('eqtl_three_types_v1/simulations/',sc,'.rds');sim<-loadv(simrel)
spec<-list(version='standard_preprocessing_v1',package=as.character(packageVersion('BayesPrism')),groups=c('Rb','Mrp','other_Rb','chrM'),exp.cells=1,outlier.cut=.01,outlier.fraction=.1,pseudo.min=1e-8,chain.length=1000L,burn.in=500L,thinning=2L,alpha=1,seed=20261013L+i,key=NULL,update.gibbs=TRUE)
cache<-file.path(out,paste0(sc,'.rds'))
if(file.exists(cache)){z<-readRDS(cache);stopifnot(identical(z$spec,spec),z$sim_md5==unname(tools::md5sum(file.path(root,simrel))),z$reference_md5==unname(tools::md5sum(file.path(root,'eqtl_bayesprism_comparison_v1/prepared_reference.rds'))));cat('Verified existing cache',sc,'\n');quit(save='no')}
ref<-BayesPrism::cleanup.genes(prep$reference,input.type='GEP',species='hs',gene.group=spec$groups,exp.cells=spec$exp.cells)
p<-BayesPrism::new.prism(reference=ref,input.type='GEP',cell.type.labels=prep$types,cell.state.labels=prep$types,key=NULL,mixture=t(sim$bulk),outlier.cut=.01,outlier.fraction=.1,pseudo.min=1e-8)
genes<-colnames(p@mixture);stopifnot(identical(rownames(p@mixture),prep$samples),all(genes%in%prep$genes))
qc<-data.frame(scenario=sc,gene=prep$genes,cleanup_retained=prep$genes%in%colnames(ref),retained=prep$genes%in%genes,nominated=prep$genes%in%sim$design$gene)
write.csv(qc,file.path(out,paste0(sc,'_preprocessing.csv')),row.names=FALSE)
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
start<-proc.time()[['elapsed']]
fit<-BayesPrism::run.prism(p,n.cores=cores,update.gibbs=TRUE,gibbs.control=spec[c('chain.length','burn.in','thinning','alpha','seed')])
dir.create(tempdir(),recursive=TRUE,showWarnings=FALSE)
x<-array(0,c(length(genes),length(prep$types),length(prep$samples)),list(genes,prep$types,prep$samples))
for(ct in prep$types){v<-BayesPrism::get.exp(fit,state.or.type='type',cell.name=ct);x[,ct,]<-t(v[prep$samples,genes,drop=FALSE])}
rho<-BayesPrism::get.fraction(fit,which.theta='final',state.or.type='type')[prep$samples,prep$types,drop=FALSE]
err<-max(abs(apply(x,c(1,3),sum)-sim$bulk[genes,,drop=FALSE])/pmax(sim$bulk[genes,,drop=FALSE],1))
stopifnot(all(is.finite(x)),all(x>=0),err<1e-6,max(abs(rowSums(rho)-1))<1e-6)
z<-list(scenario=sc,spec=spec,counts=x,proportions=rho,controls=fit@control_param,preprocessing=qc,conservation_error=err,elapsed_seconds=proc.time()[['elapsed']]-start,sim_md5=unname(tools::md5sum(file.path(root,simrel))),reference_md5=unname(tools::md5sum(file.path(root,'eqtl_bayesprism_comparison_v1/prepared_reference.rds'))),session=capture.output(sessionInfo()))
saveRDS(z,cache);writeLines(paste(unname(tools::md5sum(cache)),basename(cache)),paste0(cache,'.md5'));cat('Completed',sc,'retained genes',length(genes),'seconds',z$elapsed_seconds,'\n')
