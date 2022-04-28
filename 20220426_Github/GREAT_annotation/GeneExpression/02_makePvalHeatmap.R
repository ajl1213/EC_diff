

library(pheatmap)
library(RColorBrewer)



dir.create('Plots')

#=====================================================================

## EC Diff
inFile <- paste0('CountMatrix/EC_cRE.overlapPval.txt')
data1 <- read.table(inFile, header=T, row.names=1)
data1[data1>40] <- 40
data1 <- as.matrix(data1)

valColor <- brewer.pal(15, 'YlOrBr')

outFile <- paste0('Plots/EC.LogPval.heatmap.pdf')
pdf(outFile)
pheatmap(data1, color= valColor,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1,
display_numbers=T, number_format='%.1f')

dev.off()


## nonEC 
inFile <- paste0('CountMatrix/nonEC_cRE.overlapPval.txt')
data1 <- read.table(inFile, header=T, row.names=1)
data1[data1>40] <- 40
data1 <- as.matrix(data1)

valColor <- brewer.pal(15, 'YlOrBr')

outFile <- paste0('Plots/nonEC.LogPval.heatmap.pdf')
pdf(outFile)
pheatmap(data1, color= valColor,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1,
display_numbers=T, number_format='%.1f')

dev.off()







