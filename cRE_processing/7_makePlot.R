library(ggplot2)
library(RColorBrewer)
library(umap)


dir.create('Plots')

##=============================================================================================================================
## H3K27ac ChIP-seq read stat
mapStat <- read.table('/home/ajl1213/Projects/Endothelial/data/ChIP/Mapping/FRIP/H3K27ac.FRIP.txt', header=T, row.names=1)
sampleOrder <- rownames(mapStat)[order(mapStat$peakReads, decreasing=T)]

#
nonPeakReads <- mapStat$totalReads-mapStat$peakReads
FRIP <- mapStat$peakReads/mapStat$totalReads*100

#
a <- data.frame(sampleID=rownames(mapStat), readCount=nonPeakReads, category='nonPeakReads')
b <- data.frame(sampleID=rownames(mapStat), readCount=mapStat$peakReads, category='peakReads')
mergedDF <- rbind(a,b)

#
mergedDF$sampleID <- factor(mergedDF$sampleID, levels=rev(sampleOrder))
mergedDF$category <- factor(mergedDF$category, levels=c('nonPeakReads','peakReads'))

#
color1 <- brewer.pal(9, 'YlOrBr')[7]
color2 <- brewer.pal(9, 'Greys')[6]

#
g1 <- ggplot() +
    geom_bar(data=mergedDF, aes(x=sampleID, y=readCount, fill=category), stat='identity') +
    scale_fill_manual(values=c(color2, color1)) +
    coord_flip() +
    labs(x='', y='Number of total mapped ChIP-seq reads') +
    theme_classic()

pdf('Plots/MapStat.H3K27ac.pdf')
plot(g1)
dev.off()

#
mapStat[sampleOrder,]


##=============================================================================================================================
## UMAP plot
#
data1 <- read.table('CountMatrix/H3K27ac.RPM.qq.txt', header=T, row.names=1)
valMat <- data1[, 1:ncol(data1)]
metadata <- read.table('DataInfo.txt', header=T, row.names=1)

# umap
umap.defaults
umapObject <- umap(t(valMat), dims=1:20, n_neighbors=5, n_epochs=200)

umapDF <- umapObject$layout
colnames(umapDF) <- c('umap1','umap2')
umapDF <- cbind(umapDF, metadata)

# plot
umapDF$dataType <- factor(umapDF$dataType, level=c('hESC','Mesoderm','EarlyEC','MidEC','LateEC','FullEC','EarlyNonEC','MidNonEC','LateNonEC','FullNonEC'))
umapDF$cellType <- factor(umapDF$cellType, level=c('hESC','Mesoderm','EC','nonEC'))

g1 <- ggplot(umapDF) +
    geom_point(aes(x=umapDF$umap2 ,y=umapDF$umap1, color=umapDF$dataType, shape=umapDF$cellType)) +
#    geom_text(aes(x=umapDF$umap2 ,y=umapDF$umap1, label=rownames(umapDF)), size=1.5) +
    scale_color_manual(values=c(
    brewer.pal(8, 'Dark2')[4],
    brewer.pal(9, 'YlOrRd')[5],
    brewer.pal(9, 'PuBu')[5],
    brewer.pal(9, 'PuBu')[6],
    brewer.pal(9, 'PuBu')[7],
    brewer.pal(9, 'PuBu')[8],
    brewer.pal(9, 'Greys')[5],
    brewer.pal(9, 'Greys')[6],
    brewer.pal(9, 'Greys')[7],
    brewer.pal(9, 'Greys')[8]
    )) +
    scale_shape_manual(values=c(1,2,5,0)) +
    labs(x='UMAP2',y='UMAP1') +
    theme_classic()

pdf("Plots/H3K27ac.allPeak.UMAP.pdf", height=4, width=8, pointsize=3)
plot(g1)
dev.off()





