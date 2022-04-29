

library(ggplot2)
library(RColorBrewer)


dir.create('Plots')

adjPval_thresh <- 1e-3
fcVal_thresh <- 1.5

labelList <- c(
'Coronary_artery_disease',
'Multiple_sclerosis',
'Asthma',
'Inflammatory_bowel_disease',
'Colorectal_cancer',
'Atrial_fibrillation',
'Crohns_disease',
'Eczema',
'Rheumatoid_arthritis',
'Atopic_asthma',
'Alzheimers_disease_late_onset',
'Varicose_veins',
'Systolic_blood_pressure',
'Pulse_pressure',
'PR_interval',
'Mean_arterial_pressure'
)

##=========================================================================================================================== 
## MergedEC
# overall data
data <- read.table('/home/ajl1213/Projects/Endothelial/Analysis2/GwasEnrich/1_getGwasCatalog/RESULT/MergedEC.DataList.txt', header=T)
data$adjPval <- p.adjust(data$empiricalPval, method='BH')
insigDF <- data[-which(data$adjPval < adjPval_thresh & data$fcVal > fcVal_thresh), ]
dim(insigDF)

# Significant data (labeled)
sigDF <- read.table('DataLabel/MergedEC.SigGwasList.Label.txt', header=T)
dim(sigDF)
sigDF$categoryID <- factor(sigDF$categoryID, level=c('Disease','HematologicalTrait','OtherTrait'))

labelDF <- sigDF[sigDF$gwasID %in% labelList,]

#
g1 <- ggplot() +
    geom_point(aes(insigDF$fcVal, -log10(insigDF$adjPval), size=insigDF$nPeak), shape=1, col='darkgrey') +
    geom_point(aes(sigDF$fcVal, -log10(sigDF$adjPval), size=sigDF$nPeak, col=sigDF$categoryID)) +
    geom_text(aes(labelDF$fcVal, -log10(labelDF$adjPval), label=labelDF$gwasID, col=labelDF$categoryID), check_overlap=F, vjust=0, size=1) +
    scale_color_manual(values=c(brewer.pal(8, 'Accent')[6], brewer.pal(8, 'Accent')[5], brewer.pal(8, 'Accent')[1])) +
    scale_fill_manual(values=c(brewer.pal(8, 'Accent')[6], brewer.pal(8, 'Accent')[5], brewer.pal(8, 'Accent')[1])) +
    scale_size(breaks=seq(0, 160, by=40)) +
    geom_hline(yintercept=-log10(adjPval_thresh), linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=fcVal_thresh, linetype="dashed", color = "darkgrey") +
    ylim(0, 5.5) +
    labs(x='EC-specific cRE enrichment (foldchange over expectaion)', y='-log10(adjusted P)') +
    theme_classic()

#
pdf('Plots/MergedEC.GwasEnrich.final.pdf', width=6, height=4.5, pointsize=4)
plot(g1)
dev.off()
##=========================================================================================================================== 





