#!/usr/bin/env
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(calibrate)
library(ggExtra)

dat=read.table(args[1],sep='\t',header=T)
cols <- rep(alpha('grey', 0.2), nrow(dat))

polydat=read.table(args[2],sep='\t',header=T)
mm=subset(polydat,ID=='SNP')
extras=subset(polydat,ID=='indel')
dmm=subset(polydat,ID=='dSNP')

extrasP=extras$Position
mmP=mm$Position
dmmP=dmm$Position
Up_mmp=mmP-1
Down_mmp=mmP+1
Nei_mmp= as.vector(rbind(Up_mmp,Down_mmp)) 

cols[mmP] <- 'red' # Viral clones, Sanger calls, SNPs only (47)
cols[dmmP] <- 'blue' # Viral clones, Sanger calls, dense SNPs
cols[extrasP] <- 'black' # Viral clones, Sanger calls, diffused signal and indels
cols[Nei_mmp] <- alpha('green',0.5)

pdf(args[3],width=10, height=6)
plot(log2(dat$minDS),log2(dat$minDS/dat$minDS0),col=cols,xlab='log2(minor Diff-seq coverage frequencies)',ylab=args[4])
legend("topleft",c("SNV","dense SNV", "SNV neighbors", 'indels',"Non-variant positions") ,pch=c(20,20,20,20) ,col=c("red","blue",alpha('green',0.5),"black",alpha('grey',0.2)),cex=1)

p = ggplot() + geom_point(data=dat,aes(x=log2(minDS),y=log2(minDS/minDS0)),col=cols) + labs(x='log2(Minor allele Diff-seq coverage frequencies)',y=args[4])+ theme(panel.background = element_rect(fill ="transparent", colour = "black"),panel.grid.minor = element_blank(), panel.grid.major = element_blank())+theme(axis.text=element_text(size=16),axis.title=element_text(size=18))

ggExtra::ggMarginal(
  p, data=dat,
  margins = 'both',
  size = 5,
  col = 'black',type='histogram',
  fill = 'green'
)

p = ggplot() + geom_point(data=dat,aes(x=log2(DS),y=log2(DS/DS0)),col=cols) + labs(x='log2(Diff-seq coverage frequencies)',y=args[4])+ theme(panel.background = element_rect(fill ="transparent", colour = "black"),panel.grid.minor = element_blank(), panel.grid.major = element_blank())+theme(axis.text=element_text(size=16),axis.title=element_text(size=18))

ggExtra::ggMarginal(
  p, data=dat,
  margins = 'both',
  size = 5,
  col = 'black',type='histogram',
  fill = 'green'
)
