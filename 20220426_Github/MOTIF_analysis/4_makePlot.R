
library(ggplot2)



dir.create('Figs')
colList=c('lightyellow','yellow','orange','red','brown')


data1=read.table('StageSpecificTF/StageSpecificTfList.ValMat.final.txt', header=T)

data1$TfID=factor(data1$TfID, level=unique(data1$TfID))
data1$CellType=factor(data1$CellType, level=rev(unique(data1$CellType)))
data1$PvalLabel=factor(data1$PvalLabel, level=c('i0-10','i10-30', 'i30-60', 'i60-100', 'above100'))
data1$RnaLabel=factor(data1$RnaLabel, level=c('i0-1', 'i1-3', 'i3-5', 'i5-7', 'above7'))

unique(data1$PvalLabel)

head(data1)

g1=ggplot(data1)
g2=geom_point(aes(x=TfID, y=CellType, size=PvalLabel, col=RnaLabel))
g3=scale_color_manual(values=colList)
g4=scale_size_manual(values=c(0, 3, 5, 7, 9))
g5=labs(x='', y='')
g6=theme_classic()
g7=theme(axis.text.x = element_text(angle = 90, hjust = 1))


pdf('Figs/StageSpecificTF.pdf', width=12, height=4, pointsize=3)
plot(g1+g2+g3+g4+g5+g6+g7)
dev.off()




