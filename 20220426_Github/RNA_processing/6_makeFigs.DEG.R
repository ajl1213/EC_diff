library(ggplot2)
library(pheatmap)
library(RColorBrewer)


#===========================================================================
## load data
data1 <- read.table('DEGs/RNA.DegLabel.zscore.txt', header=T, row.names=1)
valMat <- data1[, 3:ncol(data1)]
metadata <- read.table('DataInfo.txt', header=T)


#===========================================================================
valMat <- as.matrix(valMat)

maxVal <- 3
minVal <- -1

valMat[which(valMat>maxVal)] <- maxVal
valMat[which(valMat<minVal)] <- minVal

summary(as.vector(valMat))

valColor <- brewer.pal(9, 'YlGnBu')

pdf('Plots/RNA.DEGs.zscore.pdf')

pheatmap(valMat, color= valColor,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=F,
cutree_rows=1, cutree_cols=1)
dev.off()











