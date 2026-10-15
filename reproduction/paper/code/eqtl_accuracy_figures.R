# Usage: Rscript --vanilla eqtl_accuracy_figures.R OUTPUT_FOLDER
args<-commandArgs(TRUE);o<-if(length(args))args[1]else'.';td<-file.path(o,'supplementary_tables')
source(file.path(o,'cell_type_labels.R'))
read<-function(n)read.csv(file.path(td,n));cur<-read('S20_precision_recall_curves.csv');fixed<-read('S21_precision_recall_F1.csv');matched<-read('S22_matched_empirical_FDP.csv')
wf<-c('BayesPrism standard','CDSeq2 before','CDSeq2 GAM6');labels<-c('BayesPrism','CDSeq2 (without gene-type prioritization)','CDSeq2 (with gene-type prioritization; top 40%)');cols<-c('#D55E00','#777777','#0072B2');ltys<-c(3,2,1);pchs<-c(17,15,16)
fig<-function(name,w,h,fun){for(ext in c('pdf','png')){f<-file.path(o,'figures',paste0(name,'.',ext));if(ext=='pdf')pdf(f,width=w,height=h,useDingbats=FALSE)else png(f,width=w,height=h,units='in',res=300,type='cairo');fun();dev.off()}}
lab<-function(x)mtext(x,side=3,line=1.0,adj=-.04,font=2,cex=1.25)
sharedlegend<-function(){par(oma=c(0,0,0,0),fig=c(0,1,0,.09),new=TRUE,mar=c(0,0,0,0));plot.new();legend('center',labels,col=cols,lty=ltys,lwd=2,pch=pchs,ncol=1,bty='n',cex=.72)}
pr<-function(h,r=0){plot(NA,xlim=c(0,100),ylim=c(0,103),xaxs='i',yaxs='i',xlab='Recall of all injected targets (%)',ylab='Precision (%)',main=paste('Precision-recall | ratio',h),cex.main=.96);abline(h=c(20,40,60,80,100),v=c(20,40,60,80,100),col='#EEEEEE',lwd=.7)
 for(j in seq_along(wf)){z<-subset(cur,cell_type=='Excit'&replicate==r&signal_ratio==h&workflow==wf[j]&calls>0);z<-z[order(z$q_threshold),];lines(z$recall*100,z$precision*100,col=cols[j],lty=ltys[j],lwd=1.6);q<-subset(fixed,cell_type=='Excit'&replicate==r&signal_ratio==h&workflow==wf[j]);points(q$recall*100,q$precision*100,pch=21,bg=cols[j],col='white',cex=1.2,lwd=1)}
 if(r==0)mtext('Circles: original q <= 0.05 threshold',3,line=.05,cex=.72,col='#555555')else mtext(paste('Repetition',r,'| circles: q <= 0.05'),3,line=.05,cex=.7,col='#555555')
}
fig('Figure_5_precision_recall_F1_matched_FDP',11,10.6,function(){par(mfrow=c(3,2),mar=c(4.2,4.4,3,1.0),oma=c(5.7,0,0,0),mgp=c(2.65,.8,0),las=1,cex=1.0)
 for(h in c(.5,.8)){pr(h);lab(if(h==.5)'A'else'B')}
 for(h in c(.5,.8)){z<-subset(fixed,cell_type=='Excit'&replicate==0&signal_ratio==h);z<-z[match(wf,z$workflow),];xx<-barplot(z$F1,col=cols,border=NA,ylim=c(0,.8),names.arg=rep('',3),ylab='F1 score',main=paste('F1 at q <= 0.05 | ratio',h),cex.main=.96,cex.names=.63);axis(1,at=xx,labels=c('BayesPrism','CDSeq2\n(without gene-type\nprioritization)','CDSeq2\n(with gene-type\nprioritization)'),tick=FALSE,line=1.1,cex.axis=.63);text(xx,z$F1+.035,sprintf('%.3f',z$F1),cex=.86);lab(if(h==.5)'C'else'D')}
 for(h in c(.5,.8)){m<-sapply(c(.05,.1,.2),function(cap){z<-subset(matched,cell_type=='Excit'&replicate==0&signal_ratio==h&abs(FDP_cap-cap)<1e-10);z$recall[match(wf,z$workflow)]*100});xx<-barplot(m,beside=TRUE,col=cols,border=NA,names.arg=c('5%','10%','20%'),ylim=c(0,32),ylab='Full-target recall (%)',xlab='Common empirical directional FDP cap',main=paste('Recall under FDP caps | ratio',h),cex.main=.96);text(xx,m+1,sprintf('%.1f',m),cex=.69);lab(if(h==.5)'E'else'F')}
 sharedlegend()
})
fig('Figure_S9_all_type_precision_recall_F1',11,11.4,function(){par(mfrow=c(3,2),mar=c(8.0,4.4,2.8,1),oma=c(5.7,0,0,0),mgp=c(2.65,.8,0),las=1,cex=.95);cts<-c('Astro','Excit','Inhib','Oligo','OPC','Immune','Vascular');i<-0
 for(metric in c('precision','recall','F1'))for(h in c(.5,.8)){i<-i+1;plot(NA,xlim=c(.5,7.5),ylim=c(0,1.05),xaxt='n',xlab='',ylab=switch(metric,precision='Precision',recall='Full-target recall',F1='F1 score'),main=paste(switch(metric,precision='Precision',recall='Recall',F1='F1'),'| ratio',h),cex.main=.95);axis(1,1:7,ct_name(cts,TRUE),las=2,cex.axis=.70);abline(h=seq(0,1,.2),col='#EEEEEE');for(j in seq_along(wf)){z<-subset(fixed,replicate==0&signal_ratio==h&workflow==wf[j]);z<-z[match(cts,z$cell_type),];points(1:7+(j-2)*.17,z[[metric]],pch=pchs[j],col=cols[j],cex=1.0)};lab(LETTERS[i])}
 sharedlegend()
})
fig('Figure_S10_replicate_precision_recall',11,7.6,function(){par(mfrow=c(2,3),mar=c(4.2,4.2,3,1),oma=c(5.7,0,0,0),mgp=c(2.6,.8,0),las=1,cex=.88);i<-0;for(h in c(.5,.8))for(r in 1:3){i<-i+1;pr(h,r);lab(LETTERS[i])};sharedlegend()})
