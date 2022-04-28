library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(enrichR)



dir.create('GOBP')

##===============================================================================================================================================================================
minPct <- 0.2 
minQval <- 0.05
minLog2fc <- 0

## Gene list
# GTF
gtf <-read.table('/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt', header=T) 

# snRNA-seq DEGs 
snDegTotal <- read.table('/home/ajl1213/Projects/Endothelial/data/snRNA/ColonCancer_GSE132465/2_Analysis.harmony/CellTypeDEGs/snRNA.cellTypeDEGs.txt', header=T)
EC_DF <- snDegTotal[which(snDegTotal$celltype=='EndothelialCells'),]

#
totalGeneList <- intersect(EC_DF$gene, gtf$geneID) 
EC_DF <- EC_DF[which(EC_DF$pct.1 > minPct | EC_DF$pct.2 > minPct),]
EC_DF <- EC_DF[which(abs(EC_DF$avg_log2FC)> minLog2fc),]
validGenes <- EC_DF$gene

#
snDegGenes <- snDegTotal[which(snDegTotal$celltype=='EndothelialCells' & snDegTotal$p_val_adj < minQval),]$gene
snDegGenes <- snDegGenes[snDegGenes %in% totalGeneList]
snDegGenes <- snDegGenes[snDegGenes %in% validGenes]

# Trajectory-associated genes
trajDF <- read.table('/home/ajl1213/Projects/Endothelial/data/snRNA/ColonCancer_GSE132465/3_Trajectory/pseudotime_res.all.txt', header=T)
trajGenes <- trajDF$geneID[which(trajDF$q_value < minQval)]
trajGenes <- trajGenes[trajGenes %in% totalGeneList]
trajGenes <- trajGenes[trajGenes %in% validGenes]

##===============================================================================================================================================================================
#a=listEnrichrDbs() 

#dbs <- c('GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021')
dbs <- c('GO_Biological_Process_2021')
#dbs <- c('GO_Biological_Process_2018')

# snRNA-seq DEGs
enrichTable <- enrichr(snDegGenes, dbs)
df <- NULL
for (db in dbs){
    enrichTable[[db]]$dataType <- 'snDEGs'
    enrichTable[[db]]$db <- db

    df <- rbind(df, enrichTable[[db]])
}
df.filtered <- df[which(df$P.value<0.05),]
write.table(df.filtered, 'GOBP/snDEGs.txt', row.names=F, col.names=T, sep='\t', quote=F)

# Trajectory-associated genes
enrichTable <- enrichr(trajGenes, dbs)
df <- NULL
for (db in dbs){
    enrichTable[[db]]$dataType <- 'trajGenes'
    enrichTable[[db]]$db <- db

    df <- rbind(df, enrichTable[[db]])
}
df.filtered <- df[which(df$P.value<0.05),]
write.table(df.filtered, 'GOBP/trajGenes.txt', row.names=F, col.names=T, sep='\t', quote=F)

# Stage-specific DEGs
stageDEGs <- read.table('/home/ajl1213/Projects/Endothelial/data/RNA/Processing/DEGs/RNA.DegLabel.txt', header=T)
for (cur_stage in c('EarlyEC','MidEC','LateEC','FullEC')){
    print(cur_stage)
    cur_DF <- stageDEGs[which(stageDEGs$maxDataType==cur_stage),]
    cur_stageDEGs <- cur_DF$geneID
    
    #
    a <- intersect(cur_stageDEGs, snDegGenes)
    b <- intersect(a, trajGenes)
    print(b)

    #
    enrichTable <- enrichr(cur_stageDEGs, dbs)
    df <- NULL
    for (db in dbs){
	enrichTable[[db]]$dataType <- 'gwasTargets'
	enrichTable[[db]]$db <- db

	df <- rbind(df, enrichTable[[db]])
    }
    df.filtered <- df[which(df$P.value<0.05),]
    write.table(df.filtered, paste0('GOBP/', cur_stage, '.txt'), row.names=F, col.names=T, sep='\t', quote=F)
}



