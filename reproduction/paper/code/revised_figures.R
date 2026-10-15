# Rebuild revised Figure 2 directly from the accompanying tables.
args<-commandArgs(trailingOnly=TRUE)
o<-if(length(args))args[1]else'.'
source(file.path(o,'cell_type_labels.R'))
v<-read.csv(file.path(o,'supplementary_tables/S1_profile_predictions.csv'))
cur<-read.csv(file.path(o,'supplementary_tables/S2_recovery_curves.csv'))
cts<-c('Astro','Excit','Inhib','Oligo','OPC','Immune','Vascular')
cols<-c(GAM6='#0072B2',contribution='#D55E00',expected_count='#009E73')
lt<-c(GAM6=1,contribution=2,expected_count=3)
fig<-function(name,w,h,fun){for(ext in c('pdf','png')){f<-file.path(o,'figures',paste0(name,'.',ext));if(ext=='pdf')pdf(f,width=w,height=h,family='Helvetica',useDingbats=FALSE)else png(f,width=w,height=h,units='in',res=300,type='cairo');par(mar=c(4.2,4.5,3.1,1),mgp=c(2.6,.7,0),las=1,cex=.93);fun();dev.off()}}
label<-function(x)mtext(x,side=3,line=1.25,adj=-.11,font=2,cex=1.2)
fig('Figure_2_expression_recovery',10,7.7,function(){
 layout(matrix(c(1,1,1,2,3,4),2,3,byrow=TRUE),heights=c(1.12,1))
 par(mar=c(4.2,4.5,3.3,1));plot(NA,xlim=c(.4,7.6),ylim=c(-1,1),xaxt='n',xlab='',ylab='Across-donor recovery r',main='Prioritization enriches profiles with stronger recovery',cex.main=1.02)
 abline(h=0,lty=3,col='#B0B0B0');axis(1,seq_along(cts),ct_name(cts,TRUE),cex.axis=.73);mtext('A',side=3,line=1.25,adj=-.028,font=2,cex=1.2)
 violin<-function(z,x,color){z<-z[is.finite(z)];d<-density(z,from=-1,to=1,n=512,cut=0);w<-.18*d$y/max(d$y);polygon(c(x-w,rev(x+w)),c(d$x,rev(d$x)),col=adjustcolor(color,.72),border=color,lwd=.7);q<-quantile(z,c(.25,.5,.75));segments(x,q[1],x,q[3],lwd=2.5,col='#333333');points(x,q[2],pch=21,bg='white',col='#333333',cex=.75)}
 for(i in seq_along(cts)){d<-v[v$cell_type==cts[i],];d<-d[order(-d$variation_score,d$gene),];violin(d$observed_r,i-.20,'#A0A0A0');violin(d$observed_r[1:400],i+.20,cols['GAM6'])}
 legend('bottomright',c('CDSeq2 (without gene-type prioritization)','CDSeq2 (with gene-type prioritization); top 20%'),fill=c('#BDBDBD',cols['GAM6']),bty='n',cex=.85,horiz=FALSE)
 lims<-list(Excit=c(.85,1),Inhib=c(.40,.66),Oligo=c(.40,.96))
 for(i in seq_along(lims)){ct<-names(lims)[i];par(mar=c(4.6,4.6,3.5,.8));plot(NA,xlim=c(0,100),ylim=lims[[ct]],xlab='Profiles retained (%)',ylab='Median recovery r',main=ct_name(ct),cex.main=1.0);label(LETTERS[i+1]);abline(v=seq(0,100,20),col='#EEEEEE');for(m in names(cols)){d<-subset(cur,stage=='validation'&cell_type==ct&method==m&retention>0);d<-d[order(d$retention),];lines(d$retention*100,d$median_Pearson,col=cols[m],lty=lt[m],lwd=if(m=='GAM6')2.6 else 1.8)};if(ct=='Excit')legend('bottomleft',c('CDSeq2 (with gene-type\nprioritization)','Contribution','Expected count'),col=cols,lty=lt,lwd=2,bty='n',cex=.75);mtext('Left: fewer retained | Higher: better recovery',side=3,line=.25,cex=.59,col='#555555')}
})
