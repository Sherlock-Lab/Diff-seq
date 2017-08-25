#!/usr/bin/env
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(calibrate)

# Diffseq files
dataDS1=read.table(args[1],sep='\t',header=T)
dataDS2=read.table(args[2],sep='\t',header=T)

# nextera files
dataNX1=read.table(args[3],sep='\t',header=T)
dataNX2=read.table(args[4],sep='\t',header=T)

# mark mismatches
cols <- rep(alpha('grey',0.2), nrow(dataDS1))
polydat=read.table(args[5],sep='\t',header=T)
mm=subset(polydat,ID=='SNP')
extras=subset(polydat,ID=='indel')
dmm=subset(polydat,ID=='dSNP')
extrasP=extras$Position
mmP=mm$Position
dmmP=dmm$Position
cols[mmP] <- 'red' # Viral clones, Sanger calls, SNPs only (47)
cols[dmmP] <- 'blue' # Viral clones, Sanger calls, dense SNPs
cols[extrasP] <- 'black' # Viral clones, Sanger calls, diffused signal and indels

# outfile
pdf(args[6],width=6, height=6)
par(mar=c(4,4,1,1))

# Average Diffseq replicates
DSnrT=(dataDS1$NonRef+dataDS2$NonRef)/2

# Average Nextera replicates
NXnrT=(dataNX1$NonRef+dataNX2$NonRef)/2

# log2
plot(log2(NXnrT),log2(DSnrT),col=cols,main = args[7], xlab="Nextera", ylab='Diff-seq', xlim = c(-25,0), ylim = c(-25,0),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex=1.2)
#, xlim = c(-25,-8), ylim = c(-25,-8)
abline(0,1)
legend("bottomright",c("SNPs","dense SNPs", "matched") ,pch=c(1,1,1) ,col=c("red", "blue",alpha('grey',0.2)), cex=1.2)

#legend("bottomright",c("SNPs","dense SNPs", "indels","matched") ,pch=c(1,1,1,1) ,col=c("red", "blue","black",alpha('grey',0.2)), cex=1.2)

# plot(NXnrT,DSnrT,col=cols,main = args[7], xlab="Nextera", ylab='Diff-seq', xlim = c(-30,-5), ylim = c(-30,-5),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex=1.2)
# abline(0,1)
# legend("bottomright",c("SNPs","dense SNPs", "indels","matched") ,pch=c(1,1,1,1) ,col=c("red", "blue","black",alpha('grey',0.2)), cex=1.2)

dev.off()