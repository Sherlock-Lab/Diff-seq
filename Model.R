#!/usr/bin/env
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(AUC)
library(ROCR)

pdf(args[1],width=10, height=6)
# pdf("DataFolder/Model/V4240population.pdf",width=10, height=6)

dat1=read.table(args[2],sep='\t',header=T)
# dat1=read.table("DataFolder/Model/V4240_S3CandidateParameters.txt",sep='\t',header=T)
roc1=roc(dat1$TotMinfZ+dat1$TotMinrZ,dat1$Status)
auc(roc1)
plot(roc1,cex.axis="1.4",cex.lab="1.4",col='blue')

cc = 1
K = args[3]

col8 = c('blue','red', 'yellow', 'black', 'purple', 'green', 'light blue', 'pink', 'dark green')
# for (i in args[3])
for (i in 6:K){
  cc = cc + 1
  dat1=read.table(args[i],sep='\t',header=T)
  roc1=roc(dat1$TotMinfZ+dat1$TotMinrZ,dat1$Status)
  auc(roc1)
  print (auc(roc1))
  plot(roc1,add=TRUE,col=col8[cc])
}
leg = c() 
pal = c()
ll = c()
count = 0
for (i in args[4]:args[5]){
  count = count+1
  leg[count] = args[i]
  pal[count] = col8[count]
  ll[count] = 1
}
legend("bottomright", legend = leg, col = pal, lty = ll)
dev.off()

# Rscript Model.R DataFolder/Model/V4240population.pdf DataFolder/Model/V4240_S3CandidateParameters.txt 8 9 12 DataFolder/Model/V4240_S10CandidateParameters.txt DataFolder/Model/v5248V4240_S2CandidateParameters.txt DataFolder/Model/v5248_S1CandidateParameters.txt replicate1 replicate2 diluted control
# Rscript Model.R DataFolder/Model/M4150population.pdf DataFolder/Model/M4150_S4CandidateParameters.txt 6 7 8 DataFolder/Model/M4150_S11CandidateParameters.txt replicate1 replicate2
