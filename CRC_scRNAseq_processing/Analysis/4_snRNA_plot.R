
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(dplyr)


seuset <- readRDS('SeuratObjects/CRC_snRNA.postAlign.Anno.rds')
seuset$orig.ident <- factor(seuset$orig.ident, levels=c(
'SMC01.N', 'SMC02.N', 'SMC03.N', 'SMC04.N', 'SMC05.N',
'SMC06.N', 'SMC07.N', 'SMC08.N', 'SMC09.N', 'SMC10.N',
'SMC01.T', 'SMC02.T', 'SMC03.T', 'SMC04.T', 'SMC05.T',
'SMC06.T', 'SMC07.T', 'SMC08.T', 'SMC09.T', 'SMC10.T',
'SMC11.T', 'SMC14.T', 'SMC15.T', 'SMC16.T', 'SMC17.T',
'SMC18.T', 'SMC19.T', 'SMC20.T', 'SMC21.T', 'SMC22.T',
'SMC23.T', 'SMC24.T', 'SMC25.T'
))

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
# feature plot
DefaultAssay(seuset) <- 'RNA'
dir.create('CellMarkerDistribFinal')
for(cur_gene in c('EPCAM', 'PLP1', 'DCN', 'CLDN5', 'PDGFRB', 'CD68', 'CD3D', 'CD79A', 'KIT','IGHG1','IGHA1','MS4A1','CD19')){
  print(cur_gene)
  png(paste0('CellMarkerDistribFinal/', cur_gene, '.snRNA.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(seuset, features=cur_gene, cols=viridis(256), order=T, min.cutoff='q20', max.cutoff='q98') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

##=======================================================================================================
# umap by sample
png('Plots/snRNA_umap_sample.FINAL.png', width=8, height=8, res=1000, units='in')
DimPlot(seuset, reduction='umap', group.by='orig.ident', label=F, cols=sample_colors) + my_theme + ggtitle('')
dev.off()

# sample barplot
seuset_meta <- seuset@meta.data
variable <- 'orig.ident'
clusters <- unique(seuset_meta$ident)
df <- data.frame()
for(i in 1:length(clusters)){
  idx=which(seuset_meta$ident==clusters[i])
  cur_df <- table(seuset_meta$orig.ident[idx])
  cur_df <- as.data.frame(cur_df / table(seuset_meta[[variable]])[names(cur_df)])
  cur_df$Freq <- cur_df$Freq * 1/(sum(cur_df$Freq))
  cur_df$cluster <- clusters[i]
  df <- rbind(df, cur_df)
}
df$Var1 <- factor(df$Var1, levels=c(
'SMC01.N', 'SMC02.N', 'SMC03.N', 'SMC04.N', 'SMC05.N',
'SMC06.N', 'SMC07.N', 'SMC08.N', 'SMC09.N', 'SMC10.N',
'SMC01.T', 'SMC02.T', 'SMC03.T', 'SMC04.T', 'SMC05.T',
'SMC06.T', 'SMC07.T', 'SMC08.T', 'SMC09.T', 'SMC10.T',
'SMC11.T', 'SMC14.T', 'SMC15.T', 'SMC16.T', 'SMC17.T',
'SMC18.T', 'SMC19.T', 'SMC20.T', 'SMC21.T', 'SMC22.T',
'SMC23.T', 'SMC24.T', 'SMC25.T'
))

pdf('Plots/snRNA_barplot_sample.pdf', height=4, width=9)
p <- ggplot(df, aes(y=Freq, x=cluster, fill=Var1)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=sample_colors) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
print(p)
dev.off()

##=======================================================================================================
# snRNA-seq umap
png('Plots/snRNA_umap_disease.FINAL.png', width=8, height=8, res=1000, units='in')
DimPlot(seuset, reduction='umap', group.by='disease', label=F) + scale_color_manual(values=c(brewer.pal(9, 'Set1')[2], brewer.pal(9, 'Set1')[1])) + my_theme + ggtitle('')
dev.off()

##=======================================================================================================
# QC violin plots

#
pdf(paste0('Plots/snRNA_nCount.cellType.violin.pdf'), width=8, height=4)
VlnPlot(seuset, features='nCount_RNA', cols=celltype_colors, pt.size=0, group.by='ident', ncol=1) +
  geom_boxplot(fill='white', outlier.shape=NA) + RotatedAxis() + NoLegend()
dev.off()

pdf(paste0('Plots/snRNA_nFeature.cellType.violin.pdf'), width=8, height=4)
VlnPlot(seuset, features='nFeature_RNA', cols=celltype_colors, pt.size=0, group.by='ident', ncol=1) +
  geom_boxplot(fill='white', outlier.shape=NA) + RotatedAxis() + NoLegend()
dev.off()

pdf(paste0('Plots/snRNA_nCount.sample.violin.pdf'), width=8, height=4)
VlnPlot(seuset, features='nCount_RNA', pt.size=0, group.by='orig.ident', cols=sample_colors, ncol=1) +
  geom_boxplot(fill='white', outlier.shape=NA) + RotatedAxis() + NoLegend()
dev.off()

pdf(paste0('Plots/snRNA_nFeature.sample.violin.pdf'), width=8, height=4)
VlnPlot(seuset, features='nFeature_RNA', pt.size=0, group.by='orig.ident', ncol=1, cols=sample_colors) +
  geom_boxplot(fill='white', outlier.shape=NA) + RotatedAxis() + NoLegend()
dev.off()


##=======================================================================================================
# Cell type proportions
cellTypeList <- c('EpithelialCells','EntericGlialCells','EndothelialCells','StromalCells','Pericytes','Tcells','Myelods','MastCells','IgA+IgG+PlasmaBcells','CD19+CD20+Bcells')

seuset <- subset(seuset, ident %in% cellTypeList)
seuset$ident <- factor(seuset$ident, levels=cellTypeList)
meta_list <- seuset@meta.data %>% dplyr::group_split(orig.ident)

temp <- lapply(meta_list, function(meta){
  print(table(meta$disease))
  df <- as.data.frame(meta$ident %>% table / nrow(meta))
  colnames(df) <- c('cluster', 'proportion')
  df$sampleID <- unique(meta$orig.ident)
  df$disease <- unique(meta$disease)
  df  
})
proportion_df <- Reduce(rbind, temp)
proportion_df$cluster_num <- as.numeric(proportion_df$cluster)

# relevel
proportion_df$cluster<- factor(
  as.character(proportion_df$cluster), levels=cellTypeList)

proportion_df$cluster_num <- as.numeric(proportion_df$cluster)

proportion_df$disease<- factor(proportion_df$disease, levels=rev(c('Normal','Tumor')))

# box plot
p <- ggplot(proportion_df, aes(y=proportion, x=reorder(cluster, -cluster_num), fill=disease)) +
  geom_boxplot(outlier.shape=NA, color='black') +
  coord_flip() +
  stat_compare_means(method='wilcox.test', label='p.signif', label.y=0.9) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position="bottom",
  ) + xlab('') + theme_bw()

pdf('Plots/Celltype_composition_boxplot.pdf', width=5, height=4)
plot(p)
dev.off()



