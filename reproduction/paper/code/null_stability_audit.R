# Diagnostic sensitivities only. No primary donor/normalization/selection changes.
a<-commandArgs(TRUE); v<-readRDS(a[1]);out<-a[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
types<-dimnames(v$counts$CDSeq2)[[2]];genes<-dimnames(v$counts$CDSeq2)[[1]];keys<-data.frame(gene=rep(genes,length(types)),cell_type=rep(types,each=length(genes)))
norm<-function(z,method){do.call(rbind,lapply(seq_along(types),function(k){x<-z[,k,];if(method=='logCounts')return(log2(x+.5));x<-sweep(x,2,colSums(x),'/')*1e6;if(method=='logCPM')log2(x+.5)else x}))}
center<-function(x)x-rowMeans(x)
test<-function(x,y){xx<-center(x);yy<-y-mean(y);ss<-rowSums(xx^2);r<-drop(xx%*%yy)/sqrt(ss*sum(yy^2));r[!is.finite(r)|ss<1e-20]<-NA;r<-pmax(-1,pmin(1,r));p<-2*pt(-abs(r)*sqrt((ncol(x)-2)/pmax(1-r^2,1e-15)),ncol(x)-2);q<-p.adjust(replace(p,!is.finite(p),1),'BH');list(r=r,p=p,q=q,call=is.finite(p)&q<=.05)}
selected<-rep(FALSE,nrow(keys));for(s in v$selection_sets)if(s$method=='GAM6'&&s$retention==.2)selected[s$ix]<-TRUE
rows<-profiles<-list();appendrow<-function(...)rows[[length(rows)+1]]<<-data.frame(...,stringsAsFactors=FALSE)
for(engine in names(v$counts))for(method in c('logCPM','CPM','logCounts'))for(cohort in c(39L,40L)){
 keep<-if(cohort==39)setdiff(1:40,3)else 1:40;x<-norm(v$counts[[engine]],method)[,keep];xx<-center(x);ss<-rowSums(xx^2);ok<-is.finite(ss)&ss>=1e-20;zx<-xx[ok,,drop=FALSE]/sqrt(ss[ok]);pcs<-svd(zx,nu=0,nv=3);ev<-pcs$d^2/sum(pcs$d^2);pc1<-pcs$v[,1]
 appendrow(record='expression_structure',engine=engine,normalization=method,n=cohort,trait=NA,omitted_donor=if(cohort==39)3 else 0,calls=NA,selected_calls=NA,pc1_fraction=ev[1],effective_rank=exp(-sum(ev[ev>0]*log(ev[ev>0]))),eligible_profiles=sum(ok),phenotype_pc1_r=NA)
 for(j in 1:200){
  y<-v$null_phenotypes[keep,j];t<-test(x,y)
  appendrow(record='null_trait',engine=engine,normalization=method,n=cohort,trait=j,omitted_donor=if(cohort==39)3 else 0,calls=sum(t$call),selected_calls=sum(t$call&selected),pc1_fraction=NA,effective_rank=NA,eligible_profiles=sum(ok),phenotype_pc1_r=cor(y,pc1))
  if(engine=='CDSeq2'&&method=='logCPM'&&cohort==39&&sum(t$call)>1000){
   # Profile alignment, variation and donor concentration for each burst.
   pc_r<-drop(xx%*%pc1)/sqrt(ss*sum(pc1^2));max_share<-apply(xx^2,1,max)/ss
   profiles[[length(profiles)+1]]<-data.frame(trait=j,keys,called=t$call,selected=selected,expression_sd=sqrt(ss/(cohort-1)),r=t$r,q=t$q,pc1_alignment=pc_r,max_centered_donor_share=max_share,dominant_donor=keep[max.col(xx^2,ties.method='first')])
   for(donor in keep){kk<-which(keep!=donor);tt<-test(x[,kk],y[kk]);appendrow(record='leave_one_more_donor',engine=engine,normalization=method,n=38,trait=j,omitted_donor=donor,calls=sum(tt$call),selected_calls=sum(tt$call&selected),pc1_fraction=NA,effective_rank=NA,eligible_profiles=NA,phenotype_pc1_r=NA)}
  }
 }
 message(engine,' ',method,' ',cohort)
}
write.csv(do.call(rbind,rows),file.path(out,'S39_null_stability.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,profiles),file.path(out,'S40_null_profile_diagnostics.csv'),row.names=FALSE,na='')
support<-list();for(engine in names(v$counts)){z<-v$counts[[engine]];for(j in 1:40)for(k in seq_along(types))support[[length(support)+1]]<-data.frame(engine=engine,donor_index=j,cell_type=types[k],modeled_library=sum(z[,k,j]),zero_genes=sum(z[,k,j]==0),gene_count=length(genes))}
write.csv(do.call(rbind,support),file.path(out,'S42_donor_expression_support.csv'),row.names=FALSE)
cat('Null structure, normalization and donor-omission diagnostics complete.\n')
