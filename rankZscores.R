#!/usr/bin/env
args = commandArgs(trailingOnly=TRUE)
library(calibrate)
library(ggplot2)

dat=read.table(args[1],sep='\t',header=T)
pdf(args[2],width=10, height=6)

# mark mismatches
cols <- rep(alpha('grey',0.2), nrow(dat))
polydat=read.table(args[3],sep='\t',header=T)
mm=subset(polydat,ID=='SNP')
extras=subset(polydat,ID=='indel')
dmm=subset(polydat,ID=='dSNP')
extrasP=extras$Position-1
mmP=mm$Position-1
dmmP=dmm$Position-1
cols[mmP] <- 'red' # Viral clones, Sanger calls, SNPs only (47)
cols[dmmP] <- 'blue' # Viral clones, Sanger calls, dense SNPs
cols[extrasP] <- 'black' # Viral clones, Sanger calls, diffused signal and indels

Total=dat$TotfZ + dat$TotrZ
MinorTot=dat$TotMinfZ + dat$TotMinrZ
Major=dat$MajfZ + dat$MajrZ
Minor=dat$MinfZ + dat$MinrZ
TriN1=dat$TriN1f + dat$TriN1r
TriN2=dat$TriN2f + dat$TriN2r
Rnk_T=rank(Total)
Rnk_MinT=rank(MinorTot)
Rnk_Maj=rank(Major)
Rnk_Min=rank(Minor)
Rnk_1=rank(TriN1)
Rnk_2=rank(TriN2)

plot(Rnk_T,Total,col=cols,xlab = "Rank", ylab="Total", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
legend("topleft",legend=c("SNPs","dense SNPs","match"),col=c("red","blue",alpha('grey',0.3)), pch=c(1,1,1), cex = 1.3)

plot(Rnk_MinT,MinorTot,col=cols,xlab = "Rank", ylab="Minor allele", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
legend("topleft",legend=c("SNPs","dense SNPs","match"),col=c("red","blue",alpha('grey',0.3)), pch=c(1,1,1), cex = 1.3)

plot(Rnk_Maj,Major,col=cols,xlab = "Rank", ylab="Major allele (allele dependent)", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
legend("topleft",legend=c("SNPs","dense SNPs","match"),col=c("red","blue",alpha('grey',0.3)), pch=c(1,1,1), cex = 1.3)

plot(Rnk_Min,Minor,col=cols,xlab = "Rank", ylab="Minor allele (allele dependent)", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
legend("topleft",legend=c("SNPs","dense SNPs","match"),col=c("red","blue",alpha('grey',0.3)), pch=c(1,1,1), cex = 1.3)

plot(Rnk_1,TriN1,col=cols,xlab = "Rank", ylab="Major allele (trinucleotidic context dependent)", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
legend("topleft",legend=c("SNPs","dense SNPs","match"),col=c("red","blue",alpha('grey',0.3)), pch=c(1,1,1), cex = 1.3)

plot(Rnk_2,TriN2,col=cols,xlab = "Rank", ylab="Minor allele (trinucleotidic context dependent)", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
legend("topleft",legend=c("SNPs","dense SNPs","match"),col=c("red","blue",alpha('grey',0.3)), pch=c(1,1,1), cex = 1.3)

# plot(Rnk_T,Total,col=cols,xlab = "Rank", ylab="Total", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
# legend("topleft",legend=c("SNPs","dense SNPs","indels","match"),col=c("red","blue", "black",alpha('grey',0.3)), pch=c(1,1,1,1), cex = 1.3)

# plot(Rnk_MinT,MinorTot,col=cols,xlab = "Rank", ylab="Minor allele", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
# legend("topleft",legend=c("SNPs","dense SNPs","indels","match"),col=c("red","blue", "black",alpha('grey',0.3)), pch=c(1,1,1,1), cex = 1.3)

# plot(Rnk_Maj,Major,col=cols,xlab = "Rank", ylab="Major allele (allele dependent)", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
# legend("topleft",legend=c("SNPs","dense SNPs","indels","match"),col=c("red","blue", "black",alpha('grey',0.3)), pch=c(1,1,1,1), cex = 1.3)

# plot(Rnk_Min,Minor,col=cols,xlab = "Rank", ylab="Minor allele (allele dependent)", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
# legend("topleft",legend=c("SNPs","dense SNPs","indels","match"),col=c("red","blue", "black",alpha('grey',0.3)), pch=c(1,1,1,1), cex = 1.3)

# plot(Rnk_1,TriN1,col=cols,xlab = "Rank", ylab="Major allele (trinucleotidic context dependent)", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
# legend("topleft",legend=c("SNPs","dense SNPs","indels","match"),col=c("red","blue", "black",alpha('grey',0.3)), pch=c(1,1,1,1), cex = 1.3)

# plot(Rnk_2,TriN2,col=cols,xlab = "Rank", ylab="Minor allele (trinucleotidic context dependent)", cex.axis="1.7",cex.lab="1.7", cex=1.3, xlim = c(0,2700))
# legend("topleft",legend=c("SNPs","dense SNPs","indels","match"),col=c("red","blue", "black",alpha('grey',0.3)), pch=c(1,1,1,1), cex = 1.3)


dev.off()