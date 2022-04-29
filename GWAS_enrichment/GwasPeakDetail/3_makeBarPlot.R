
library(ggplot2)
library(RColorBrewer)




data1=read.table('GwasPeakData/GwasPeakData.txt', header=T)

data1$gwasID <- factor(data1$gwasID, level=rev(unique(data1$gwasID)))
data1$cellType <- factor(data1$cellType, level=rev(c('EarlyEC','MidEC','LateEC','FullEC')))


g1 <- ggplot(data=data1, aes(x=data1$gwasID, y=data1$nPeak, fill=data1$cellType)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c(brewer.pal(8, 'Dark2')[1], brewer.pal(8, 'Dark2')[5], brewer.pal(8, 'Dark2')[6], brewer.pal(8, 'Dark2')[2])) +
    coord_flip() +
    labs(x=NULL, y='Peak count') +
    theme_classic() 

pdf('Plots/GwasPeakData.barplot.pdf', height=6, width=6, pointsize=3)
plot(g1)
dev.off()



