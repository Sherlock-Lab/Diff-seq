#!/usr/bin/env
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(calibrate)

data1=read.table(args[1],sep='\t',header=T)
pdf(args[2],width=6, height=6)
par(mar=c(4,4,1,1))

boxplot(Af~REFALT,data=data1, xlab='Ref/Alt', ylab="Fraction of N at position", ylim=c(0,1), cex.lab=1.2, cex.axis=0.9, col="green")
boxplot(Cf~REFALT,data=data1, xlab='Ref/Alt', ylab="Fraction of N at position", ylim=c(0,1), cex.lab=1.2, cex.axis=0.9, col="blue")
boxplot(Gf~REFALT,data=data1, xlab='Ref/Alt', ylab="Fraction of N at position", ylim=c(0,1), cex.lab=1.2, cex.axis=0.9, col="black")
boxplot(Tf~REFALT,data=data1, xlab='Ref/Alt', ylab="Fraction of N at position", ylim=c(0,1), cex.lab=1.2, cex.axis=0.9, col="red")
boxplot(Ar~REFALT,data=data1, xlab='Ref/Alt', ylab="Fraction of N at position", ylim=c(0,1), cex.lab=1.2, cex.axis=0.9, col="green")
boxplot(Cr~REFALT,data=data1, xlab='Ref/Alt', ylab="Fraction of N at position", ylim=c(0,1), cex.lab=1.2, cex.axis=0.9, col="blue")
boxplot(Gr~REFALT,data=data1, xlab='Ref/Alt', ylab="Fraction of N at position", ylim=c(0,1), cex.lab=1.2, cex.axis=0.9, col="black")
boxplot(Tr~REFALT,data=data1, xlab='Ref/Alt', ylab="Fraction of N at position", ylim=c(0,1), cex.lab=1.2, cex.axis=0.9, col="red")

dev.off()
