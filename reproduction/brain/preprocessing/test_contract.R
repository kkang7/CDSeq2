self<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);here<-dirname(normalizePath(self))
source(file.path(here,'contract.R'))
g<-c('g1','g2');k<-c('A','B');donors<-c('d1','d2')
x<-matrix(1:4,2,dimnames=list(g,k));a<-array(1:8,c(2,2,2),list(g,k,donors))
b<-list(gene_panel=g,major_celltypes=k,pseudobulk_counts=apply(a,c(1,3),sum),
 truth_expression_counts=setNames(list(x,x),donors),truth_cell_counts=x,
 target_metadata=data.frame(sample=donors,donor=donors),reference_subjects='r1')
s<-SingleCellExperiment::SingleCellExperiment(assays=list(counts=methods::as(x,'dgCMatrix')),
 colData=data.frame(subject=c('r1','r1'),cellType=k,row.names=k))
v<-list(counts=a,metadata=b$target_metadata,n_cells=x,pooled_reference=x,reference_donors='r1')
r<-list(reference_mean=x/10,reference_donors='r1',fit=list(mcmc_iterations=300L))
f<-tempfile()
z<-rosmap_signatures(b,s,v,r);write.csv(data.frame(input=names(z),sha256=unname(z)),f,row.names=FALSE)
rosmap_verify_objects(b,s,v,r,f)
for(change in c('count','donor','fit','reference')) {
 bb<-b;ss<-s;vv<-v;rr<-r
 if(change=='count')vv$counts[1]<-vv$counts[1]+1
 if(change=='donor')vv$metadata$donor[1]<-'different'
 if(change=='fit')rr$fit$mcmc_iterations<-301L
 if(change=='reference')SummarizedExperiment::assay(ss,'counts')[1,1]<-9
 stopifnot(inherits(try(rosmap_verify_objects(bb,ss,vv,rr,f),silent=TRUE),'try-error'))
}
unlink(f);cat('PASS: valid objects accepted; changed count, donor, fit and reference rejected.\n')
