library(ggplot2)
library(RColorBrewer)
library(umap)


dir.create('Plots')


##=============================================================================================================================
## RNA-seq read stat
mapStat <- read.table('MapStat.txt', header=T, row.names=1)

#
countMat <- read.table('CountMatrix/RNA.rawCount.txt', header=T, row.names=1)
countMat <- countMat[,2:ncol(countMat)]
geneCounts <- colSums(countMat)
sampleOrder <- names(geneCounts)[order(geneCounts, decreasing=T)]

#
UnmappedToGenes <- mapStat$UniqMapRead-geneCounts
geneBodyFrac <- geneCounts/mapStat$UniqMapRead*100

#
a <- data.frame(sampleID=rownames(mapStat), readCount=UnmappedToGenes, category='UnmappedToGenes')
b <- data.frame(sampleID=names(geneCounts), readCount=geneCounts, category='MappedToGeneBody')
mergedDF <- rbind(a,b)

#
mergedDF$sampleID <- factor(mergedDF$sampleID, levels=rev(sampleOrder))
mergedDF$category <- factor(mergedDF$category, levels=c('UnmappedToGenes','MappedToGeneBody'))

#
color1 <- brewer.pal(9, 'YlGnBu')[7]
color2 <- brewer.pal(9, 'Greys')[6]

#
g1 <- ggplot() +
    geom_bar(data=mergedDF, aes(x=sampleID, y=readCount, fill=category), stat='identity') +
    scale_fill_manual(values=c(color2, color1)) +
    coord_flip() +
    labs(x='', y='Number of total mapped RNA-seq reads') +
    theme_classic()

pdf('Plots/MapStat.RNA.pdf')
plot(g1)
dev.off()

#
DF <- data.frame(mapStat, geneBodyFrac, MappedToGenebody=geneCounts)
orderedDF <- DF[names(geneCounts)[order(geneCounts, decreasing=T)],]
orderedDF



##=============================================================================================================================
## UMAP plot
#
data1 <- read.table('CountMatrix/RNA.TPM.qq.txt', header=T, row.names=1)
valMat <- data1[, 2:ncol(data1)]
metadata <- read.table('DataInfo.txt', header=T, row.names=1)

# umap
umap.defaults
umapObject <- umap(t(valMat), dims=1:10, n_neighbors=10, n_epochs=300)

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

pdf("Plots/RNA.allGene.UMAP.pdf", height=4, width=8, pointsize=3)
plot(g1)
dev.off()








