#!/home/ajl1213/anaconda2/bin/R


library(Seurat)
library(dplyr)
library(data.table)
library(mclust)
library(harmony)


data1 <- data.frame(fread('/home/ajl1213/Projects/Endothelial/data/snRNA/ColonCancer_GSE132465/1_downLoad/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz', header=T, sep='\t', stringsAsFactors=F))
rownames(data1) <- data1$Index
valMat <- data1[,2:ncol(data1)]

seuset <- CreateSeuratObject(counts = valMat, project = "CRC.snRNA")

# Data alignment
seuset@meta.data$orig.ident <- factor(seuset@meta.data$orig.ident, levels=c(
'SMC01.N', 'SMC02.N', 'SMC03.N', 'SMC04.N', 'SMC05.N',
'SMC06.N', 'SMC07.N', 'SMC08.N', 'SMC09.N', 'SMC10.N',
'SMC01.T', 'SMC02.T', 'SMC03.T', 'SMC04.T', 'SMC05.T',
'SMC06.T', 'SMC07.T', 'SMC08.T', 'SMC09.T', 'SMC10.T',
'SMC11.T', 'SMC14.T', 'SMC15.T', 'SMC16.T', 'SMC17.T',
'SMC18.T', 'SMC19.T', 'SMC20.T', 'SMC21.T', 'SMC22.T',
'SMC23.T', 'SMC24.T', 'SMC25.T'))

table(seuset@meta.data$orig.ident)

disease <- c()
disease[grep('T', seuset@meta.data$orig.ident)] <- 'Tumor'
disease[grep('N', seuset@meta.data$orig.ident)] <- 'Normal'

disease <- factor(disease, level=c('Normal','Tumor'))
seuset@meta.data$disease <- disease

##========================================================================================================================================================================================
# Load predetermined cell label
metaData <- read.table('/home/ajl1213/Projects/Endothelial/data/snRNA/ColonCancer_GSE132465/1_downLoad/GSE132465_GEO_processed_CRC_10X_cell_annotation.mod.txt', header=T)
metaData <- data.frame(row.names=metaData$barcodeID, celltype=metaData$celltype, subtype=metaData$subcelltype, curatedcelltype=metaData$curatedcelltype)
metaData <- metaData[rownames(metaData) %in% rownames(seuset@meta.data),]
seuset <- AddMetaData(seuset, metadata=metaData)

##========================================================================================================================================================================================
# standard seurat pipeline
seuset <- NormalizeData(seuset)
seuset <- FindVariableFeatures(seuset, selection='vst', nfeatures=3000)
seuset <- ScaleData(seuset, verbose = FALSE)
seuset <- RunPCA(seuset, npcs = 50, verbose = FALSE)

# Harmony
seuset <- RunHarmony(object = seuset, group.by.vars = 'orig.ident')

#
seuset <- RunUMAP(seuset, reduction = "harmony", dims = 1:50)
seuset <- FindNeighbors(seuset, reduction = "harmony", dims = 1:50)
seuset <- FindClusters(seuset, resolution= 1.5)

# Save the work
dir.create('SeuratObjects')
saveRDS(seuset, 'SeuratObjects/CRC_snRNA.postAlign.rds')
#


