#!/home/ajl1213/anaconda2/bin/R


library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)


seuset <- readRDS('SeuratObjects/CRC_snRNA.postAlign.rds')

table(seuset@meta.data$orig.ident)
table(seuset@meta.data$disease)

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

##=======================================================================================================
# UMAP
dir.create('Plots')
# by cluster
png('Plots/snRNA_umap_cluster.png', width=10, height=10, res=1000, units='in')
DimPlot(seuset, reduction='umap', label=T) + my_theme + ggtitle('')
dev.off()
# by pathology
png('Plots/snRNA_umap_disease.png', width=10, height=10, res=1000, units='in')
DimPlot(seuset, reduction='umap', group.by='disease', label=F) + scale_color_manual(values=c(brewer.pal(9, 'Set1')[2], brewer.pal(9, 'Set1')[1])) + my_theme + ggtitle('')
dev.off()
# by sample
png('Plots/snRNA_umap_sample.png', width=10, height=10, res=1000, units='in')
DimPlot(seuset, reduction='umap', group.by='orig.ident', label=F, cols=sample_colors) + my_theme + ggtitle('')
dev.off()
# cell type
png('Plots/snRNA_umap_celltype.png', width=10, height=10, res=1000, units='in')
DimPlot(seuset, reduction='umap', group.by='celltype', label=T, cols=celltype_colors) + my_theme + ggtitle('')
dev.off()
# sub cell type
png('Plots/snRNA_umap_subcelltype.png', width=10, height=10, res=1000, units='in')
DimPlot(seuset, reduction='umap', group.by='subtype', label=T, cols=rainbow(length(unique(seuset$subtype)))) + my_theme + ggtitle('')
dev.off()
# curated cell type
png('Plots/snRNA_umap_curatedcelltype.png', width=10, height=10, res=1000, units='in')
DimPlot(seuset, reduction='umap', group.by='curatedcelltype', label=T, cols=rainbow(length(unique(seuset$curatedcelltype)))) + my_theme + ggtitle('')
dev.off()


##=======================================================================================================
# CellMarker Distribution
seuset@active.assay <- 'RNA'
dir.create('CellMarkerDistrib')

# Epithelial
png('CellMarkerDistrib/Epithelial.EPCAM.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='EPCAM', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
# EntericGlialCells
png('CellMarkerDistrib/EntericGlial.PLP1.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='PLP1', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
# Stromal
png('CellMarkerDistrib/Stromal.DCN.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='DCN', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
# Endo
png('CellMarkerDistrib/Endo.ENG.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='ENG', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Endo.PECAM1.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='PECAM1', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Endo.CLDN5.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='CLDN5', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
# Pericyte
png('CellMarkerDistrib/Peri.RGS5.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='RGS5', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
# SMC
png('CellMarkerDistrib/SMC.PDGFRB.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='PDGFRB', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
# Myeloid
png('CellMarkerDistrib/Myeloid.CD68.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='CD68', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
# Tcell
png('CellMarkerDistrib/Tcell.CD3D.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='CD3D', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
# Bcell
png('CellMarkerDistrib/Bcell.CD79A.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='CD79A', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()
# Mast
png('CellMarkerDistrib/Mast.KIT.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='KIT', cols=viridis(256), order=T, min.cutoff='q5', max.cutoff='q98') + my_theme + ggtitle('')
dev.off()









