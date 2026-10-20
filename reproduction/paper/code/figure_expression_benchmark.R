# Reproduce Figure S4; A-B data/geometry are unchanged, workflow labels updated.
# Rebuild new supplementary benchmarks and revised main Figure 3 from bundled CSVs.
args<-commandArgs(TRUE);o<-if(length(args))args[1]else'.';td<-file.path(o,'supplementary_tables')
source(file.path(o,'cell_type_labels.R'))
read<-function(s)read.csv(file.path(td,s));s<-read('S12_ground_truth_accuracy.csv');cts<-c('Astro','Excit','Inhib','Oligo','OPC','Immune','Vascular')
fig<-function(name,w,h,fun){for(ext in c('pdf','png')){file<-file.path(o,'figures',paste0(name,'.',ext));if(ext=='pdf')pdf(file,width=w,height=h,useDingbats=FALSE)else png(file,width=w,height=h,units='in',res=300,type='cairo');par(mar=c(4.4,4.6,3,1),mgp=c(2.7,.8,0),las=1);fun();dev.off()}}
lab<-function(x)mtext(x,3,line=1.1,adj=if(par("mfrow")[2]==1)-.028 else -.12,font=2,cex=1.2)
# Gene/scenario distributions accompany nine scenario-level fraction metrics.
profiles<-read('S14_expression_accuracy.csv')
profiles<-subset(profiles,scenario!='direct'&oracle=='observed')
violin<-function(v,x,col,bounds){v<-v[is.finite(v)];if(length(unique(v))>1){d<-density(v,n=512,from=max(bounds[1],min(v)),to=min(bounds[2],max(v)),cut=0);w<-.16*d$y/max(d$y);polygon(c(x-w,rev(x+w)),c(d$x,rev(d$x)),col=adjustcolor(col,.65),border=col,lwd=.8)}else segments(x-.12,v[1],x+.12,v[1],col=col,lwd=2);q<-quantile(v,c(.25,.5,.75));segments(x,q[1],x,q[3],lwd=2.3,col='#333333');points(x,q[2],pch=21,bg='white',cex=.65);if(length(v)<=10)points(x+seq(-.08,.08,length.out=length(v)),sort(v),pch=16,cex=.45,col=col)}
membership<-read('S23_S6_GAM6_membership.csv')
key<-function(d)paste(d$scenario,d$cell_type,d$gene,sep='|')
cd<-subset(profiles,engine=='CDSeq2')
ix<-match(key(cd),key(membership));stopifnot(!anyNA(ix))
selected<-cd[membership$selected[ix],]
fig('Figure_S4_ground_truth_comparison',10,8.4,function(){
 par(mfrow=c(3,1),mar=c(4.2,4.7,3.0,1),cex=.92)
 z<-subset(s,scenario!='direct'&oracle=='observed')
 for(j in 1:2){
  metric<-c('fraction_r','fraction_MAE','Pearson')[j]
  bounds<-if(j==2)c(0,ceiling(max(z$fraction_MAE,na.rm=TRUE)*100*1.08))else if(j==3)c(-1,1)else c(min(0,floor(min(z$fraction_r,na.rm=TRUE)*10)/10),1)
  plot(NA,xlim=c(.5,7.5),ylim=bounds,xaxt='n',xlab='',ylab=c('RNA-fraction correlation','RNA-fraction MAE (pp)','Across-donor expression r')[j],main=c('RNA-fraction correlation: 9 scenario estimates','RNA-fraction absolute error: 9 scenario estimates','Individual-expression recovery: all gene/scenario profiles')[j],cex.main=.92)
  axis(1,1:7,ct_name(cts,TRUE),cex.axis=.73);lab(LETTERS[j])
  for(e in c('CDSeq2','BayesPrism'))for(i in 1:7){
   col<-if(e=='CDSeq2')'#0072B2'else'#D55E00';off<-if(e=='CDSeq2')-.19 else .19
   v<-z[z$engine==e&z$cell_type==cts[i],metric]*if(j==2)100 else 1
   violin(v,i+off,col,bounds)
  }
  if(j==1)legend('bottomleft',c('CDSeq2 (without gene-type prioritization)','BayesPrism'),fill=c('#0072B2','#D55E00'),bty='n',horiz=TRUE,cex=.75)
 }
 # Membership subsets existing, common-panel correlations.
 par(mar=c(5.1,4.7,3.0,1))
 plot(NA,xlim=c(.5,7.5),ylim=c(-1,1),xaxt='n',xlab='',ylab='Across-donor expression r',main='Individual-expression recovery: inference and optional prioritization',cex.main=.92)
 axis(1,1:7,ct_name(cts,TRUE),cex.axis=.73);lab('C');abline(h=0,lty=3,col='#BBBBBB')
 palette<-c('#D55E00','#999999','#0072B2')
 for(w in 1:3)for(i in 1:7){
  d<-if(w==3)selected else subset(profiles,engine==c('BayesPrism','CDSeq2')[w])
  v<-d$Pearson[d$cell_type==cts[i]];v<-v[is.finite(v)];x<-i+c(-.27,0,.27)[w]
  if(length(unique(v))>1){
   den<-density(v,n=512,from=max(-1,min(v)),to=min(1,max(v)),cut=0)
   width<-.115*den$y/max(den$y)
   polygon(c(x-width,rev(x+width)),c(den$x,rev(den$x)),col=adjustcolor(palette[w],.65),border=palette[w],lwd=.8)
  }
  q<-quantile(v,c(.25,.5,.75));segments(x,q[1],x,q[3],lwd=2.3,col='#333333');points(x,q[2],pch=21,bg='white',cex=.65)
 }
 lx<-c(.62,2.7,5);points(lx,rep(-.85,3),pch=22,bg=palette,cex=.8);text(lx+.12,rep(-.85,3),c('BayesPrism\n100%','CDSeq2 (without gene-type\nprioritization); 100%','CDSeq2 (with gene-type\nprioritization); 20.0%'),adj=c(0,.5),cex=.72)
 counts<-read('S24_S6_recovery_counts.csv'); cc<-subset(counts,scenario=='ALL'&cell_type=='Excit'); nn<-cc$n_retained[match(c('BayesPrism','CDSeq2','CDSeq2+GAM6'),cc$workflow)]
 mtext(paste('Retained per type across 9 scenarios (same workflow order):',paste(format(nn,big.mark=',',trim=TRUE),collapse='; ')),1,line=3.7,cex=.73)
})
