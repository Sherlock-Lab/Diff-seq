#!/usr/bin/env
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(heatmap3)

dat=read.table(args[1],sep='\t',header=T)
datMB = with(dat, dat[order(RefAlt, triN, Position) , ])

pdf(args[2],width=8, height=12)

my_pal <- colorRampPalette(c("yellow", "orange", "brown"))(n = 299)

heatmap3(datMB[c(4:9)], Rowv = NA, Colv = NA, scale="none", cexRow = 0.01 + 1/log10(nrow(datMB)), labRow=c(datMB$Position), main="Biases")

heatmap3(datMB[c(4:9)], Rowv = NA, Colv = NA, scale="none", cexRow = 0.01 + 1/log10(nrow(datMB)), labRow=datMB$RefAlt, main="Biases")

heatmap3(datMB[c(4:9)], Rowv = NA, Colv = NA, scale="none", cexRow = 0.01 + 1/log10(nrow(datMB)), labRow=datMB$triN, main="Biases")

my_palette <- colorRampPalette(c("yellow", "red"))(n = 299)
heatmap3(datMB[c(10:15)], Rowv = NA, Colv = NA, scale="none", cexRow = 0.01 + 1/log10(nrow(datMB)), labRow=c(datMB$Position), col=my_palette, main="Log2 Coverage Frequencies")

dev.off()