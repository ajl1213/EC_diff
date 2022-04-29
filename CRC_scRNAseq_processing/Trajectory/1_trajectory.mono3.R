
library(Seurat)
library(monocle3)
library(harmony)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)


dir.create('Plots')

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

nCtrl <- 10
nDisease <- 23
sample_colors <- c(rev(brewer.pal(9, 'YlGnBu'))[1:nCtrl], colorRampPalette(c('brown','red','purple','orange'))(nDisease))
celltype_colors <- c(brewer.pal(9, 'Set1')[1:5], brewer.pal(9, 'Set1')[7:8], brewer.pal(8, 'Set2')[1:3])
disease_colors <- c(brewer.pal(9, 'Set1')[2], brewer.pal(9, 'Set1')[1])
subcelltype_colors <- c(brewer.pal(8, 'Dark2'))

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

cds <- as.cell_data_set(EC_c)
cds <- cluster_cells(cds, cluster_method='louvain')
cds <- learn_graph(cds)

#
png('Plots/partition.png', width=4.5, height=8, res=1000, units='in')
p <- plot_cells(cds, label_branch_points=F, label_principal_points=F, show_trajectory_graph=F, label_leaves=F, label_roots=F, label_groups_by_cluster=F, label_cell_groups=T, color_cells_by = "partition", cell_size = 1.2) + my_theme
plot(p)
dev.off()

png('Plots/cluster.png', width=4.5, height=8, res=1000, units='in')
p <- plot_cells(cds, label_branch_points=F, label_principal_points=F, show_trajectory_graph=F, label_leaves=F, label_roots=F, label_groups_by_cluster=F, label_cell_groups=T, color_cells_by = "cluster", cell_size = 1.2) + my_theme
plot(p)
dev.off()

#
png('Plots/sample.png', width=4.5, height=8, res=1000, units='in')
p <- plot_cells(cds, label_branch_points=F, label_principal_points=F, show_trajectory_graph=F, label_leaves=F, label_roots=F, label_groups_by_cluster=F, label_cell_groups=F, color_cells_by = "orig.ident", cell_size = 1.2) + scale_color_manual(values=sample_colors) + my_theme
plot(p)
dev.off()

png('Plots/sample.png', width=4.5, height=8, res=1000, units='in')
p <- plot_cells(cds, label_branch_points=F, label_principal_points=F, show_trajectory_graph=F, label_leaves=F, label_roots=F, label_groups_by_cluster=F, label_cell_groups=F, color_cells_by = "orig.ident", cell_size = 1.2) + scale_color_manual(values=sample_colors) + my_theme
plot(p)
dev.off()

cds$disease <- factor(cds$disease, levels=c('Normal','Tumor'))
png('Plots/disease.png', width=4.5, height=8, res=1000, units='in')
p <- plot_cells(cds, label_branch_points=F, label_principal_points=F, show_trajectory_graph=F, label_leaves=F, label_roots=F, label_groups_by_cluster=F, label_cell_groups=F, color_cells_by = "disease", cell_size = 1.2) +
    scale_color_manual(values=disease_colors) + my_theme
plot(p)
dev.off()

cds$sub_ident <- factor(cds$sub_ident, levels=c('TipLikeECs', 'StalkLikeECs', 'ProliferatingECs', 'LymphaticECs', 'Pericytes'))
png('Plots/subEC_celltype.png', width=4.5, height=8, res=1000, units='in')
p <- plot_cells(cds, label_branch_points=F, label_principal_points=F, show_trajectory_graph=F, label_leaves=F, label_roots=F, label_groups_by_cluster=F, label_cell_groups=F, color_cells_by = "sub_ident", cell_size = 1.2) +
    scale_color_manual(values=subcelltype_colors) + my_theme
plot(p)
dev.off()

# select root cluster
root_group <- colnames(cds)[which(cds@clusters$UMAP$clusters==c('13','16'))];
cds <- order_cells(cds, root_cells = root_group);
png('Plots/pseudotime.png', width=4.5, height=8, res=1000, units='in')
p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=F, label_leaves=F, label_branch_points=F, label_principal_points = F, label_roots=F, graph_label_size = 0.3, cell_size = 1.2) + my_theme
#p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=F, label_leaves=F, label_branch_points=F, label_principal_points = F, label_roots=F, graph_label_size = 0.3, cell_size = 1.2) 
plot(p)
dev.off()

##===================================================================================================================================================================================================
# Obtain metadata
table(partitions(cds))
idx <- which(partitions(cds)==1)
metadata <- data.frame(barcodeID=colnames(cds)[idx], cds@colData[idx,], pseudotime=cds@principal_graph_aux$UMAP$pseudotime[idx])
metadata <- metadata[order(metadata$pseudotime),]
rownames(metadata) <- NULL

write.table(metadata, 'pseudotime.metadata.txt', row.names=F, col.names=T, sep='\t', quote=F)

##===================================================================================================================================================================================================
# Obtain trajectory genes
pseudotime_res <- graph_test(cds, neighbor_graph='principal_graph', cores=4)
pseudotime_res <- data.frame(geneID=rownames(pseudotime_res), pseudotime_res) ; rownames(pseudotime_res) <- NULL; pseudotime_res$V1 <- NULL
write.table(pseudotime_res, 'pseudotime_res.all.txt', row.names=F, col.names=T, sep='\t', quote=F)



