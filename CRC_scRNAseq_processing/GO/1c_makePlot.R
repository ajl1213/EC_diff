library(pheatmap)
library(RColorBrewer)
library(viridis)


dir.create('Plots')

maxPval <- 20

#
col_range <- brewer.pal(9, 'YlOrBr')

data1 <- read.table('StageGOBP.merged.selected.txt', header=T, row.names=1)
print(data1)

data1[data1 > maxPval] <- maxPval

pdf('Plots/StageDEG.GOBP.pdf')
pheatmap(data1,
color= col_range,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1,
display_numbers = T, number_format='%.1f'
)
dev.off()




