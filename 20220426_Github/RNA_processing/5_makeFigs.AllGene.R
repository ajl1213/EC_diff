library(ggplot2)
library(pheatmap)
library(RColorBrewer)



dir.create('Plots')


#===========================================================================
## load data
data1 <- read.table('CountMatrix/RNA.TPM.qq.zscore.txt', header=T, row.names=1)
valMat <- data1[, 2:ncol(data1)]
metadata <- read.table('DataInfo.txt', header=T)

#===========================================================================
data1_adj <- cor(as.matrix(valMat))
print(summary(as.vector(data1_adj)))


## Set annotation
annoCol <- data.frame(dataType=metadata$dataType)
rownames(annoCol) <- metadata$sampleID
annoCol

setColColor <- list(
dataType=c(
hESC=brewer.pal(8, 'Dark2')[4],
Mesoderm=brewer.pal(9, 'YlOrRd')[5],
EarlyEC=brewer.pal(9, 'PuBu')[5],
MidEC=brewer.pal(9, 'PuBu')[6],
LateEC=brewer.pal(9, 'PuBu')[7],
FullEC=brewer.pal(9, 'PuBu')[8],
EarlyNonEC=brewer.pal(9, 'Greys')[5],
MidNonEC=brewer.pal(9, 'Greys')[6],
LateNonEC=brewer.pal(9, 'Greys')[7],
FullNonEC=brewer.pal(9, 'Greys')[8]
)
)

valColor <- brewer.pal(9, 'BuPu')

## Heatmap
pdf('Plots/RNA.allGene.corHeatmap.pdf')
pheatmap(data1_adj, annotation_col=annoCol,
color= valColor, annotation_colors=setColColor,
cluster_rows=T, cluster_cols=T, 
show_colnames=T, show_rownames=F, 
cutree_rows=1, cutree_cols=1)
dev.off()


#===========================================================================
## PCA
data1 <- read.table('CountMatrix/RNA.TPM.qq.txt', header=T, row.names=1)
valMat <- data1[,2:ncol(data1)]
data.pca <- as.data.frame(predict(prcomp(t(valMat)),scale=FALSE))

# Variance in each PC
a <- (predict(prcomp(t(valMat))))^2
b <- colSums(a)
totalVar <- sum(b)
var_PC <- b/totalVar*100

# Combine Data
df <- data.frame(data.pca, DataType=metadata$dataType, CellType=metadata$cellType)

# Make Plot
TargetFactor1 <- factor(df$DataType, level=c('hESC','Mesoderm','EarlyEC','MidEC','LateEC','FullEC','EarlyNonEC','MidNonEC','LateNonEC','FullNonEC'))
TargetFactor2 <- factor(df$CellType, level=c('hESC','Mesoderm','EC','nonEC'))

TargetFactor1
TargetFactor2

g1 <- ggplot(df) +
    geom_point(aes(x=df$PC2 ,y=df$PC1, color=TargetFactor1, shape=TargetFactor2)) +
    #geom_text(aes(x=df$PC2 ,y=df$PC1, label=rownames(df)), size=1.5) +
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
    scale_shape_manual(values=c(1, 2, 5, 0)) +
    labs(x=paste0('PC2 (', var_PC['PC2'],'%)'), y=paste0('PC1 (', var_PC['PC1'],'%)')) +
    theme_classic()

pdf("Plots/RNA.allGene.PCA.pdf", height=3, width=6, pointsize=3)
plot(g1)
dev.off()






