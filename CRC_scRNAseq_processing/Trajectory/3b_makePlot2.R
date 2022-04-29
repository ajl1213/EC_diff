
library(Seurat)
library(harmony)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(viridis)


##===================================================================================================================================================================================================
my_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  legend.position="none",
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

disease_colors <- c(brewer.pal(8, 'Dark2')[1], brewer.pal(8, 'Dark2')[4])
subEC_colors <- c(brewer.pal(8, 'Accent')[5], brewer.pal(8, 'Accent')[7], brewer.pal(8, 'Accent')[8])
heatmap_color <- colorRampPalette(c('lightblue','white','red'))(12)

##===================================================================================================================================================================================================
# Density histogram by trajectory
metadata <- read.table('pseudotime.metadata.txt', header=T)
metadata$disease <- factor(metadata$disease, levels=c('Normal','Tumor'))
metadata$sub_ident <- factor(metadata$sub_ident, levels=c('TipLikeECs', 'StalkLikeECs', 'ProliferatingECs'))

print(table(metadata$disease))
print(table(metadata$sub_ident))

#
p <- ggplot(metadata) +
    geom_density(aes(x=pseudotime, fill=disease), alpha=0.3) +
    geom_rug(aes(x=pseudotime, y=0, col=disease), alpha=0.5, position=position_jitter(height=0)) +
    scale_color_manual(values=disease_colors) +
    scale_fill_manual(values=disease_colors) +
    theme_classic() +
    labs(title=NULL, x='Pseudotime', y='Density') +
    theme(legend.position='bottom')

pdf('Plots/pseudotime.densityHist.disease.pdf', height=3, width=9, pointsize=3)
plot(p)
dev.off()

#
idx <- 1
for (i in c('TipLikeECs', 'StalkLikeECs', 'ProliferatingECs')){
    cur_data <- metadata[which(metadata$sub_ident==i),]

    p <- ggplot() +
	geom_density(data=metadata, aes(x=pseudotime, fill=disease), alpha=0.3) +
	geom_rug(data=cur_data, aes(x=pseudotime, y=0, col=sub_ident), alpha=0.5, position=position_jitter(height=0)) +
	scale_color_manual(values=subEC_colors[idx]) +
	scale_fill_manual(values=disease_colors) +
	theme_classic() +
	labs(title=i, x='Pseudotime', y='Density') +
	theme(legend.position='bottom')

    pdf(paste0('Plots/pseudotime.densityHist.subEC.', i,'.pdf'), height=3, width=9, pointsize=3)
    plot(p)
    dev.off()

    idx <- idx+1
}

##===================================================================================================================================================================================================
minVal <- -1
maxVal <- 1.5

# 
DF <- read.table(paste0('ValMat/targetGene.log2fc.order.smooth.txt'), header=T, row.names=1)
tmpLabel <- c()
tmpLabel[DF$log2 < 0] <- 'DOWN'
tmpLabel[DF$log2 > 0] <- 'UP'
print(table(tmpLabel))

valMat <- DF[,5:ncol(DF)]

# make heatmap
valMat[valMat < minVal] <- minVal
valMat[valMat > maxVal] <- maxVal

png(paste0('Plots/pseudotime.zscore.heatmap.png'), width=8, height=8, res=1000, units='in')
pheatmap(valMat, color= heatmap_color,
cluster_rows=F, cluster_cols=F,
show_colnames=F, show_rownames=F,
treeheight_row=0, treeheight_col=0,
legend=F
)
dev.off()

## label
minLog2fc <- -0.4
maxLog2fc <- 0.4

DF <- read.table(paste0('ValMat/targetGene.log2fc.order.smooth.txt'), header=T, row.names=1)
valMat <- DF[,4:ncol(DF)]
DF$log2fc[DF$log2fc > maxLog2fc] <- maxLog2fc
DF$log2fc[DF$log2fc < minLog2fc] <- minLog2fc

# make heatmap
annoRow <- data.frame(peakLabel=DF$peakLabel, DEG=DF$DEG, trajGene=DF$trajGene, log2fc=DF$log2fc)
rownames(annoRow) <- rownames(DF)

setAnnoColor <- list(
peakLabel = c(EarlyEC=brewer.pal(8, 'Dark2')[2], MidEC=brewer.pal(8, 'Dark2')[6], LateEC=brewer.pal(8, 'Dark2')[5], FullEC=brewer.pal(8, 'Dark2')[1]),
DEG = c(DEG_intersect=brewer.pal(8, 'Dark2')[8], none='white'),
trajGene = c(Traj_intersect=brewer.pal(8, 'Dark2')[8], none= 'white'),
log2fc =  colorRampPalette(c('blue','black','red'))(30)
)

valMat[valMat < minVal] <- minVal
valMat[valMat > maxVal] <- maxVal

pdf(paste0('Plots/pseudotime.zscore.heatmap.label.pdf'), width=8, height=8)
pheatmap(valMat, 
color= heatmap_color, annotation_row = annoRow, annotation_colors = setAnnoColor,
cluster_rows=F, cluster_cols=F,
show_colnames=F, show_rownames=T,
treeheight_row=0, treeheight_col=0,
legend=T
)
dev.off()



##===================================================================================================================================================================================================
seuset <- readRDS('/home/ajl1213/Projects/Endothelial/data/snRNA/ColonCancer_GSE132465/2_Analysis.harmony/SeuratObjects/CRC_snRNA.postAlign.Anno.rds')

#
cur_celltype <- 'EndothelialCells'
EC_c <- subset(seuset, idents = cur_celltype)
EC_c <- FindVariableFeatures(EC_c, selection='vst', nfeatures=1000)
EC_c <- ScaleData(EC_c, verbose = FALSE)
EC_c <- RunPCA(EC_c, npcs = 50, verbose = FALSE)
EC_c <- RunHarmony(object = EC_c, group.by.vars = 'orig.ident', project.dim = F)
EC_c <- RunUMAP(EC_c, reduction = "harmony", dims = 1:5)


dir.create('UMAPPlot')
DefaultAssay(EC_c) <- 'RNA'
#for (geneID in c('APLN','APLNR','COL4A1','COL4A2','LAMC1')){
#for (geneID in c('EMP1', 'PLAT', 'SH3BGRL3', 'NRP2', 'SOX17', 'APLN', 'APLNR', 'TM4SF1')){
for (geneID in c('TM4SF18','TP53I11','UBE2J1')){
    print(geneID)
    png(paste0('UMAPPlot/', geneID, '.snRNA.png'), width=4.5, height=8, res=1000, units='in')
    #p <- FeaturePlot(EC_c, features=geneID, cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q50') + my_theme + ggtitle('')
    #p <- FeaturePlot(EC_c, features=geneID, cols=viridis(256), order=T, min.cutoff='q10', max.cutoff='q90') + my_theme + ggtitle('')
    #p <- FeaturePlot(EC_c, features=geneID, cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q95') + my_theme + ggtitle('')
    p <- FeaturePlot(EC_c, features=geneID, cols=viridis(256), order=T, min.cutoff='q20', max.cutoff='q75') + my_theme + ggtitle('')
    plot(p)
    dev.off()
}




