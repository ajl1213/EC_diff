library(pheatmap)
library(RColorBrewer)
library(viridis)



dir.create('Plots')

#=====================================================================

inFile<-paste0('GREAT_GOBP/EC_cRE.GOTABLE.txt')
data1<-read.table(inFile, header=T, row.names=1)
data1[data1>250] <- 250
data1<-as.matrix(data1)

valColor<-brewer.pal(15, 'GnBu')

outFile<-paste0('Plots/EC_cRE.GoCombined.heatmap.pdf')
pdf(outFile)
pheatmap(data1, color= valColor,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1,
display_numbers=T, number_format='%.1f'
)
dev.off()



