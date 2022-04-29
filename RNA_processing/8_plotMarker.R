library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(viridis)


##============================================================================
## marker list
## hESC
# POU5F1 (OCT4)
# SOX2
# NANOG
# UTF1

## ME
# MIXL1
# TBXT (T, BRACHYURY)
# EOMES

## EC 
# CDH5 (CD144, VE-cadherin)
# VWF (vWF)
# KDR (VEGFR)
# CD34
# ENG (CD105)
# CLDN5

## SMC
# SMC1A (alpha-SMC)
# PDGFRB (CD140b)
# MYH10 (Myosin IIB)
# TAGLN (SM22a)
##============================================================================
sampleList <- c('hESC', 'ME', 'EC01', 'EC02', 'EC04', 'EC06', 'EC08', 'EC12', 'EC24', 'EC48', 'nonEC01', 'nonEC02', 'nonEC04', 'nonEC06', 'nonEC08', 'nonEC12', 'nonEC24', 'nonEC48')
geneList <- c('POU5F1','SOX2','NANOG','UTF1','MIXL1','TBXT','EOMES','CDH5','VWF','KDR','CD34','ENG','CLDN5','SMC1A','PDGFRB','MYH10','TAGLN')


# Heatmap
data1 <- read.table('CountMatrix/RNA.TPM.qq.stage.txt', header=T)

valMat <- NULL
for (geneID in geneList){
    tmpVec <- data1[which(data1$geneID==geneID),]
    valMat <- rbind(valMat, tmpVec)
}

rownames(valMat) <- valMat$geneID
valMat <- valMat[, 3:ncol(valMat)]

#
get_zscore <- function(x){(x - mean(x)) /sd(x)}
zmat <- apply(valMat, 1, get_zscore)   

maxVal <- 3
minVal <- -1
zmat[which(zmat>maxVal)] <- maxVal
zmat[which(zmat<minVal)] <- minVal

#
valColor <- rev(rocket(15))
pdf('Plots/markers.TPM.zscore.heatmap.pdf', width=5, height=8, pointsize=3)
pheatmap(zmat, color= valColor,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1)
dev.off()



