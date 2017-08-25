#!/usr/bin/env

args = commandArgs(trailingOnly=TRUE)
library(ggplot2)

IDf=read.table(args[1],sep="\t",header=T)
mm=subset(IDf,ID2=='SNP')
extras=subset(IDf,ID2=='indel')
extrasP=extras$Position
mmP=mm$Position

dat=read.table(args[2],sep="\t",header=T)
cols <- rep(alpha('grey', 0.2), nrow(dat))

for (i in c(0:args[4])){
  cols[mmP+2668*i] <- 'red' # Viral clones, Sanger calls, SNPs only
  cols[extrasP+2668*i] <- 'black' # Viral clones, indels
}

# Makes stacked graphs of the coverage frequencies (figure 3A)
pdf(args[3],width=10, height=6)
par(mar=c(4,4,1,1))
ggplot()+geom_bar(data=dat,aes(x=Position,y=Frequency),stat="identity",color=cols)+facet_grid(Library~.,scales="free",space="free")+labs(x="Reference position",y = "Diff-seq Coverage Frequency") +theme_bw(base_size = 16) + theme(panel.background = element_rect(fill ="transparent", colour = "black"),panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(axis.text=element_text(size=8),axis.title=element_text(size=16),legend.text=element_text(size=14))
dev.off()