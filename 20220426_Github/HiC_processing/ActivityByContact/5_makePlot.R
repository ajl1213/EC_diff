

library(ggplot2)
library(RColorBrewer)
library(viridis)



dir.create('Plots')

cellTypeList <- c('Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC')
celltype_colors <- c(brewer.pal(9, 'YlOrRd')[5], brewer.pal(8, 'Dark2')[2], brewer.pal(8, 'Dark2')[6], brewer.pal(8, 'Dark2')[5], brewer.pal(8, 'Dark2')[1])

#===============================================================================================================================================================================
nGroup <- 10
idx <- 1

# Increasing order of contact score
for (cellType in cellTypeList){
    print(cellType)
    data1 <- read.table(paste0('RnaEffect/', cellType, '.ContactValuePerGene.txt'), header=T)
    data1Adj <- data1[order(data1$totalContactVal, decreasing=F),]

    intervalLabel <- c(rep(1:nGroup, each=as.integer(nrow(data1Adj)/nGroup)), rep(nGroup, times=nrow(data1Adj)-length(rep(1:nGroup, each=as.integer(nrow(data1Adj)/nGroup)))))
    df <- data.frame(data1Adj, intervalLabel=as.factor(intervalLabel))

    g <- ggplot() +
        geom_boxplot(aes(df$intervalLabel, log2(df$rnaCount+1), fill=df$intervalLabel), outlier.shape=NA, notch=TRUE) +
        scale_fill_manual(values=colorRampPalette(c('white', celltype_colors[idx]))(nGroup+2)[3:(nGroup+2)]) +
        labs(title=cellType, x='Genes (increasing order of contact score)', y='Log2RnaExpression') +
        theme_classic() +
        theme(legend.position='none') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

    pdf(paste0('Plots/', cellType, '.ContactValRNA.pdf'), width=6, height=3, pointsize=3)
    plot(g)
    dev.off()

    idx <- idx +1
}





