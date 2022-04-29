#!/home/ajl1213/anaconda2/bin/R



library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)

seuset <- readRDS('SeuratObjects/CRC_snRNA.postAlign.rds')

##=======================================================================================================
# theme:
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
# color:
nCtrl <- 10
nDisease <- 23
sample_colors <- c(rev(brewer.pal(9, 'YlGnBu'))[1:nCtrl], colorRampPalette(c('brown','red','purple','orange'))(nDisease))
celltype_colors <- c(brewer.pal(9, 'Set1')[1:5], brewer.pal(9, 'Set1')[7:8], brewer.pal(8, 'Set2')[1:3])
subcelltype_colors <- c(brewer.pal(8, 'Dark2'))

##=======================================================================================================
#
# EpithelialCells - 0, 2, 3, 16, 17, 21, 27, 31, 33, 36
# EntericGlialCells - 28
# EndothelialCells - 18, 34, 37
# StromalCells - 13, 15
# SmoothMuscleCells - 23
# Tcells - 1, 5, 7, 8, 10, 11, 14, 22, 32
# Myelods - 9, 12, 19, 35
# MastCells - 30
# IgA+IgG+PlasmaBcells - 6, 20, 25, 26
# CD19+CD20+Bcells - 4 
# undetermined - 24, 29

##=======================================================================================================
# celltype assignment
seusetAnno <- SetIdent(seuset, WhichCells(seuset, idents=c(0, 2, 3, 16, 17, 21, 27, 31, 33, 36)), value='EpithelialCells')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(28)), value='EntericGlialCells')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(18, 34, 37)), value='EndothelialCells')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(13, 15)), value='StromalCells')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(23)), value='SmoothMuscleCells')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(1, 5, 7, 8, 10, 11, 14, 22, 32)), value='Tcells')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(9, 12, 19, 35)), value='Myelods')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(30)), value='MastCells')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(6, 20, 25, 26)), value='IgA+IgG+PlasmaBcells')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(4)), value='CD19+CD20+Bcells')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(24, 29)), value='undetermined')

seusetAnno@active.ident <- factor(seusetAnno@active.ident, level=c('EpithelialCells','EntericGlialCells','EndothelialCells','StromalCells','SmoothMuscleCells','Tcells','Myelods','MastCells','IgA+IgG+PlasmaBcells','CD19+CD20+Bcells', 'undetermined'))
seusetAnno$ident <- seusetAnno@active.ident

##=======================================================================================================
# Sub-EC population 
library(harmony)

EC_c <- subset(seusetAnno, idents='EndothelialCells')

#EC_c <- NormalizeData(EC_c)
EC_c <- FindVariableFeatures(EC_c, selection='vst', nfeatures=1000)
EC_c <- ScaleData(EC_c, verbose = FALSE)
EC_c <- RunPCA(EC_c, npcs = 50, verbose = FALSE)
EC_c <- RunHarmony(object = EC_c, group.by.vars = 'orig.ident', project.dim = F)
EC_c <- RunUMAP(EC_c, reduction = "harmony", dims = 1:5)
EC_c <- FindNeighbors(EC_c, reduction = "harmony", dims = 1:5)
EC_c <- FindClusters(EC_c, resolution= 1.0, graph.name='RNA_snn')

#
dir.create('subECs')

# Tip-like ECs: RGCC, RAMP3
# Stalk-like ECs: ACKR1, SELP
# Proliferating ECs: BIRC5, CKS1B, MKI67, PCNA
# Lymphatic ECs: LYVE1, PROX1
# Pericytes: RGS5

png('subECs/cluster.png', width=5, height=5, res=1000, units='in')
DimPlot(EC_c, reduction='umap', label=T) + my_theme + ggtitle('')
dev.off()

png('subECs/sample.png', width=5, height=5, res=1000, units='in')
DimPlot(EC_c, reduction='umap', group.by='orig.ident', label=F, cols=sample_colors) + my_theme + ggtitle('')
dev.off()

DefaultAssay(EC_c) <- 'RNA'
for(cur_gene in c('RGCC', 'RAMP3', 'ACKR1', 'SELP', 'BIRC5', 'CKS1B', 'MKI67', 'PCNA', 'LYVE1', 'PROX1','RGS5','PDGFRB')){
  print(cur_gene)
  png(paste0('subECs/', cur_gene, '.snRNA.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(EC_c, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q10', max.cutoff='q90') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

# cell-type assignment
# TipLikeECs - 0, 2, 3, 5, 7, 12
# StalkLikeECs - 1, 8, 9
# ProliferatingECs - 10
# LymphaticECs - 6
# Pericytes - 4, 11

#
EC_c.anno <- SetIdent(EC_c, WhichCells(EC_c, idents=c(0, 2, 3, 5, 7, 12)), value='TipLikeECs')
EC_c.anno <- SetIdent(EC_c.anno, WhichCells(EC_c, idents=c(1, 8, 9)), value='StalkLikeECs')
EC_c.anno <- SetIdent(EC_c.anno, WhichCells(EC_c, idents=c(10)), value='ProliferatingECs')
EC_c.anno <- SetIdent(EC_c.anno, WhichCells(EC_c, idents=c(6)), value='LymphaticECs')
EC_c.anno <- SetIdent(EC_c.anno, WhichCells(EC_c, idents=c(4, 11)), value='Pericytes')

EC_c.anno$sub_ident <- EC_c.anno@active.ident
EC_c.anno.filter <- subset(EC_c.anno, idents=c('TipLikeECs', 'StalkLikeECs', 'ProliferatingECs', 'LymphaticECs', 'Pericytes'))
EC_c.anno.filter@active.ident <- factor(EC_c.anno.filter@active.ident, levels=c('TipLikeECs', 'StalkLikeECs', 'ProliferatingECs', 'LymphaticECs', 'Pericytes'))
EC_c.anno.filter$sub_ident <- EC_c.anno.filter@active.ident

##
png('subECs/subEC.png', width=5, height=5, res=1000, units='in')
DimPlot(EC_c.anno.filter, reduction='umap', group.by='sub_ident', cols=subcelltype_colors, label=F) + my_theme + ggtitle('')
dev.off()

png('subECs/subEC.label.png', width=5, height=5, res=1000, units='in')
DimPlot(EC_c.anno.filter, reduction='umap', group.by='sub_ident', cols=subcelltype_colors, label=T) + my_theme + ggtitle('')
dev.off()

DefaultAssay(EC_c.anno.filter) <- 'RNA'
for(cur_gene in c('RGCC','RAMP3','ACKR1','SELP','BIRC5','CKS1B','MKI67','PCNA','LYVE1','PROX1','RGS5','PDGFRB')){
  print(cur_gene)
  png(paste0('subECs/', cur_gene, '.snRNA.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(EC_c.anno.filter, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q0', max.cutoff='q90') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

##=======================================================================================================
# Add subEC info
idx <- match(rownames(EC_c.anno@meta.data), rownames(seusetAnno@meta.data))
seusetAnno$sub_ident <- as.character(seusetAnno$ident)
seusetAnno$sub_ident[idx] <- as.character(EC_c.anno$sub_ident)

# Filter
cellList <- c('EpithelialCells','EntericGlialCells','EndothelialCells','StromalCells','SmoothMuscleCells','Tcells','Myelods','MastCells','IgA+IgG+PlasmaBcells','CD19+CD20+Bcells')
seusetAnno.filter <- subset(seusetAnno, subset= ident %in% cellList)
seusetAnno.filter$ident <- factor(seusetAnno.filter$ident, levels=cellList)
seusetAnno.filter@active.ident <- seusetAnno.filter$ident
seusetAnno.filter$disease <- factor(seusetAnno.filter$disease, levels=c('Normal','Tumor'))

##=======================================================================================================
# Umap
png('Plots/snRNA_umap_celltype.FINAL.png', width=6, height=6, res=1000, units='in')
DimPlot(seusetAnno.filter, reduction='umap', group.by='ident', cols=celltype_colors, label=F) + my_theme + ggtitle('')
dev.off()

png('Plots/snRNA_umap_celltype.FINAL.label.png', width=6, height=6, res=1000, units='in')
DimPlot(seusetAnno.filter, reduction='umap', group.by='ident', cols=celltype_colors, label=T) + my_theme + ggtitle('')
dev.off()

# Cell number report
seusetAnno.filter$orig.ident <- factor(seusetAnno.filter$orig.ident, levels=c(
'SMC01.N', 'SMC02.N', 'SMC03.N', 'SMC04.N', 'SMC05.N',
'SMC06.N', 'SMC07.N', 'SMC08.N', 'SMC09.N', 'SMC10.N',
'SMC01.T', 'SMC02.T', 'SMC03.T', 'SMC04.T', 'SMC05.T',
'SMC06.T', 'SMC07.T', 'SMC08.T', 'SMC09.T', 'SMC10.T',
'SMC11.T', 'SMC14.T', 'SMC15.T', 'SMC16.T', 'SMC17.T',
'SMC18.T', 'SMC19.T', 'SMC20.T', 'SMC21.T', 'SMC22.T',
'SMC23.T', 'SMC24.T', 'SMC25.T'))

nCellTable <- table(seusetAnno.filter$ident, seusetAnno.filter$orig.ident)
sumVec <- apply(nCellTable, 2, sum)
fracTable<-NULL
for (i in 1:nrow(nCellTable)){
    tmpVec=nCellTable[i,]
    fracVec=tmpVec/sumVec*100
    fracTable=rbind(fracTable, fracVec)
}
sumVec
rownames(fracTable) <- rownames(nCellTable)
pdf('Plots/CellFrac.rna.bar.pdf', height=5, width=10, pointsize=2)
barplot(fracTable, col=celltype_colors, legend=rownames(fracTable), ylab='Proportion (%)', las=2)
dev.off()

#
table(seusetAnno.filter$disease)
table(seusetAnno.filter$ident)
table(seusetAnno.filter$ident)/sum(table(seusetAnno.filter$ident))*100


##=======================================================================================================
# Save the work
saveRDS(seusetAnno.filter, 'SeuratObjects/CRC_snRNA.postAlign.Anno.rds')





