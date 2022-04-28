
library(ggplot2)



dir.create('Plots')

adjPval_thresh <- 1e-3
fcVal_thresh <- 1.5
maxFcVal <- 3

##=========================================================================================================================== 
## MergedEC
##
data <- read.table('RESULT/MergedEC.DataList.txt', header=T)
data$adjPval <- p.adjust(data$empiricalPval, method='BH')
data$fcVal[data$fcVal > maxFcVal] <- maxFcVal

dataLabel <- data[which(data$adjPval < adjPval_thresh & data$fcVal > fcVal_thresh), ]
dim(dataLabel)
dataLabel

##
g1 <- ggplot() +
    geom_point(aes(data$fcVal, -log10(data$adjPval), size=data$nPeak), shape=1) +
    geom_text(aes(dataLabel$fcVal, -log10(dataLabel$adjPval), label=dataLabel$gwasID), check_overlap=F, vjust=0, size=1, col='red') +
    scale_size(breaks=seq(0, 180, by=60)) +
    geom_hline(yintercept=-log10(adjPval_thresh), linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=fcVal_thresh, linetype="dashed", color = "darkgrey") +
    ylim(0, 5.5) +
    labs(x='EC-specific cRE enrichment (Foldchange over expectaion)', y='-log10(adjusted P)') +
    theme_classic() 

##
pdf('Plots/MergedEC.GwasEnrich.pdf', width=8, height=5, pointsize=4)
plot(g1)
dev.off()


##
write.table(dataLabel, 'RESULT/MergedEC.SigGwasList.txt', col.names=T, row.names=F, sep='\t', quote=F)



