#!/usr/bin/env
args = commandArgs(trailingOnly=TRUE)

dat=read.table(args[1],sep='\t',header=T)

pdf(args[2],width=10,height=6)
par(mar=c(4,4,1,1))
L=length(dat$DilutionFactor)
pt <- rep(19, L)
cols=rep("red", L)
blu = c() 
count = 0
for (i in c(1,6,8,9)){
  count = count+1
  blu[count] = i
}
cols[blu] <- 'blue'

plot(log2(dat$DilutionFactor),log2(dat$DS/dat$NX),xlab="log2(dilution factor)", ylab = "log2(diff-seq/nextera non-reference allele frequencies (average for SNP positions))", pch=pt,col=cols, cex=1.5,
cex.lab=1.5, cex.axis=1.5)
legend("bottomright",legend=c('Dilution series A',"Dilution series B"), col=c("red","blue"), pch=c(19,19), cex=1.5)

plot(dat$DilutionFactor,dat$DS/dat$NX,xlab="dilution factor", ylab = "diff-seq/nextera non-reference allele frequencies (average for SNP positions))", pch=pt,col=cols, cex=1.5,
cex.lab=1.5, cex.axis=1.5)
legend("topright",legend=c('Dilution series A',"Dilution series B"), col=c("red","blue"), pch=c(19,19), cex=1.5)

dev.off()