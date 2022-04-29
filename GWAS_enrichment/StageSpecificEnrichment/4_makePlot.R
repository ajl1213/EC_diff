library(pheatmap)
library(RColorBrewer)
library(viridis)


dir.create('Plots')


cellTypeList <- c('hESC','Mesoderm','EarlyEC','MidEC','LateEC','FullEC')
gwasList <- c(
'Crohns_disease',
'Asthma',
'Inflammatory_bowel_disease',
'Multiple_sclerosis',
'Eczema',
'Atrial_fibrillation',
'Atopic_asthma',
'Coronary_artery_disease',
'Varicose_veins',
'Alzheimers_disease_late_onset',
'Colorectal_cancer',
'Ulcerative_colitis',
'Rheumatoid_arthritis'
)

minPval <- 1e-6

#
data1 <- read.table('RESULT/MergedDataList.txt', header=T)
data1$empiricalPval[data1$empiricalPval < minPval] <- minPval
data1$adjPval <- p.adjust(data1$empiricalPval, method='BH')
write.table(data1, file='RESULT/MergedDataList.adjPval.txt', row.names=F, col.names=T, sep="\t", quote=F)

#
data1$logQval <- -log10(data1$adjPval)

df <- NULL
for (gwasID in gwasList){
    valList <- c()
    for (cellType in cellTypeList){
	logQval <- data1$logQval[which(data1$gwasID==gwasID & data1$cellType==cellType)]
	valList <- c(valList, logQval)
    }
    df <- rbind(df, valList)
}

colnames(df) <- cellTypeList
rownames(df) <- gwasList

#
valColor <- c('white', brewer.pal(9, 'OrRd'))

pdf('Plots/gwasEnrich.byStage.logQval.pdf')
pheatmap(df, 
color=valColor,
cluster_rows=T, cluster_cols=F,
show_colnames=T, show_rownames=T, 
breaks=seq(0, max(df), length.out=9),
display_numbers = T, number_format='%.1f'
)
dev.off()



#
a <- data1[which(data1$adjPval < 0.05 & data1$adjPval > 0.01),]
print(a)
b <- data1[which(data1$adjPval < 0.01 & data1$adjPval > 0.001),]
print(b)
c <- data1[which(data1$adjPval < 0.001),]
print(c)




