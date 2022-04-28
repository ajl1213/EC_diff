
library(Vennerable)
library(GeneOverlap)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(viridis)


minPct <- 0.2
minQval <- 0.05
minLog2fc <- 0

##===================================================================================================================================================================================================
## Gene overlap

# Get gene list
# GTF
gtf <-read.table('/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt', header=T) 

# snRNA-seq DEGs 
snDegTotal <- read.table('/home/ajl1213/Projects/Endothelial/data/snRNA/ColonCancer_GSE132465/2_Analysis.harmony/CellTypeDEGs/snRNA.cellTypeDEGs.txt', header=T)
EC_DF <- snDegTotal[which(snDegTotal$celltype=='EndothelialCells'),]

#
totalGeneList <- intersect(EC_DF$gene, gtf$geneID) 
#EC_DF <- EC_DF[which(EC_DF$pct.1 > minPct & EC_DF$pct.2 > minPct),]
EC_DF <- EC_DF[which(EC_DF$pct.1 > minPct | EC_DF$pct.2 > minPct),]
EC_DF <- EC_DF[which(abs(EC_DF$avg_log2FC)> minLog2fc),]
validGenes <- EC_DF$gene

#
snDegGenes <- snDegTotal[which(snDegTotal$celltype=='EndothelialCells' & snDegTotal$p_val_adj < minQval),]$gene
snDegGenes <- snDegGenes[snDegGenes %in% totalGeneList]
snDegGenes <- snDegGenes[snDegGenes %in% validGenes]
print(length(snDegGenes))

# trajectory-associated genes
trajDF <- read.table('pseudotime_res.all.txt', header=T)
trajGenes <- trajDF$geneID[which(trajDF$q_value < minQval)]
trajGenes <- trajGenes[trajGenes %in% totalGeneList]
trajGenes <- trajGenes[trajGenes %in% validGenes]
print(length(trajGenes))

# GWAS-SNP target genes
gwasDF <- read.table('/home/ajl1213/Projects/Endothelial/Analysis2/GwasTarget/02_GetTargetGenes/TargetGeneList/Colorectal_cancer.MergedTarget.abcScore.txt', header=T)
gwasGenes <- unique(gwasDF$geneID)
gwasGenes <- gwasGenes[gwasGenes %in% totalGeneList]
print(length(gwasGenes))

# VENN
data.list <- list(snDegGenes, trajGenes, gwasGenes)
names(data.list) <- c('snDEGs', 'trajGenes','gwasGenes')
venn.ds <- Venn(data.list)

pdf('Plots/snRNA.GwasTarget.venn.pdf')
#plot(vennData, doWeights = TRUE, type = "ChowRuskey",)
plot(venn.ds, doWeights = TRUE, type = "circles")
dev.off()

#
overlap_obj <- newGeneOverlap(snDegGenes, gwasGenes, genome.size = length(totalGeneList))
overlap_obj <- testGeneOverlap(overlap_obj)

print(overlap_obj@pval)
print(overlap_obj@odds.ratio)
a <- overlap_obj@intersection
print(length(a))

#
overlap_obj <- newGeneOverlap(trajGenes, gwasGenes, genome.size = length(totalGeneList))
overlap_obj <- testGeneOverlap(overlap_obj)

print(overlap_obj@pval)
print(overlap_obj@odds.ratio)
b <- overlap_obj@intersection
print(length(b))

# make output
dir.create('GenesOverlap')

df <- data.frame(geneID=a, label='DEG_intersect')
write.table(df, file='GenesOverlap/DEG_intersect.txt', row.names=F, col.names=T, sep="\t", quote=F)

df <- data.frame(geneID=b, label='Traj_intersect')
write.table(df, file='GenesOverlap/Traj_intersect.txt', row.names=F, col.names=T, sep="\t", quote=F)


##===================================================================================================================================================================================================
#
dir.create('ValMat')

targetgenelist <- union(a, b)
seuset <- readRDS('/home/ajl1213/Projects/Endothelial/data/snRNA/ColonCancer_GSE132465/2_Analysis.harmony/SeuratObjects/CRC_snRNA.postAlign.Anno.rds')
seuset@active.ident <- factor(seuset$sub_ident)

EC_seuset <- subset(seuset, idents =c('TipLikeECs', 'StalkLikeECs', 'ProliferatingECs'))
valMat <- EC_seuset@assays$RNA@data

zScore <- function(x){(x - mean(x)) /sd(x)}
zmat <- apply(valMat, 1, zScore)
zmat <- t(zmat)
zmat <- zmat[targetgenelist, ]

df <- data.frame(geneID=rownames(zmat), zmat)
write.table(df, file='ValMat/zmat.txt', row.names=F, col.names=T, sep="\t", quote=F)


##===================================================================================================================================================================================================
## Gene overlap with stage-speicific DEGs
stageDEGs <- read.table('/home/ajl1213/Projects/Endothelial/data/RNA/Processing/DEGs/RNA.DegLabel.txt', header=T)

for (cur_stage in c('hESC','Mesoderm','EarlyEC','MidEC','LateEC','FullEC')){
    print(cur_stage)

    cur_DF <- stageDEGs[which(stageDEGs$maxDataType==cur_stage),]
    cur_stageDEGs <- cur_DF$geneID

    data.list <- list(snDegGenes, trajGenes, cur_stageDEGs)
    names(data.list) <- c('snDEGs', 'trajGenes', paste0(cur_stage, '_DEGs'))
    venn.ds <- Venn(data.list)

    pdf(paste0('Plots/', cur_stage, '_DEGs.venn.pdf'))
    #plot(vennData, doWeights = TRUE, type = "ChowRuskey",)
    plot(venn.ds, doWeights = TRUE, type = "circles")
    dev.off()

    #
    overlap_obj <- newGeneOverlap(snDegGenes, cur_stageDEGs, genome.size = length(totalGeneList))
    overlap_obj <- testGeneOverlap(overlap_obj)

    print(paste0('snDEG_pval: ', overlap_obj@pval))
    print(paste0('snDEG_OR: ', overlap_obj@odds.ratio))

    #
    overlap_obj <- newGeneOverlap(trajGenes, cur_stageDEGs, genome.size = length(totalGeneList))
    overlap_obj <- testGeneOverlap(overlap_obj)

    print(paste0('trajGene_pval: ', overlap_obj@pval))
    print(paste0('trajGene_OR: ', overlap_obj@odds.ratio))

}

## Heatmap
library(pheatmap)

DF <- NULL
for (cur_stage in c('hESC','Mesoderm','EarlyEC','MidEC','LateEC','FullEC','EarlyNonEC','MidNonEC','LateNonEC','FullNonEC')){
    print(cur_stage)
    cur_DF <- stageDEGs[which(stageDEGs$maxDataType==cur_stage),]
    cur_stageDEGs <- cur_DF$geneID

    #
    overlap_obj <- newGeneOverlap(snDegGenes, cur_stageDEGs, genome.size = length(totalGeneList))
    overlap_obj <- testGeneOverlap(overlap_obj)

    snDEG_pval <- overlap_obj@pval

    #
    overlap_obj <- newGeneOverlap(trajGenes, cur_stageDEGs, genome.size = length(totalGeneList))
    overlap_obj <- testGeneOverlap(overlap_obj)

    traj_pval <- overlap_obj@pval

    #
    vals <- c(snDEG_pval, traj_pval)
    DF <- rbind(DF, vals)

}

colnames(DF) <- c('snDEG_pval','traj_pval')
rownames(DF) <- c('hESC','Mesoderm','EarlyEC','MidEC','LateEC','FullEC','EarlyNonEC','MidNonEC','LateNonEC','FullNonEC')

logDF <- -log10(DF)
print(logDF)

logDF[logDF > 40] <- 40
logDF <- as.matrix(logDF)

valColor <- brewer.pal(15, 'YlGnBu')

outFile <- paste0('Plots/StageDEGs.LogPval.heatmap.pdf')
pdf(outFile)
pheatmap(logDF, color= valColor,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1,
display_numbers=T, number_format='%.1f')
dev.off()









