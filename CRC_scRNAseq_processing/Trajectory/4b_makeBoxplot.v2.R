
library(ggplot2)
library(RColorBrewer)



data1 <- read.table('BySubtype/zval.bysubtype.txt', header=T)
data1$cms <- factor(data1$cms, levels=c('Normal','CMS1','CMS2','CMS3','CMS4'))

##======================================================================================================== 
#

DOWN <- data1[which(data1$geneLabel=='down'),]
UP <- data1[which(data1$geneLabel=='up'),]

# DOWN
Normal <- DOWN$zval[which(DOWN$cms=='Normal')]
CMS1 <- DOWN$zval[which(DOWN$cms=='CMS1')]
CMS2 <- DOWN$zval[which(DOWN$cms=='CMS2')]
CMS3 <- DOWN$zval[which(DOWN$cms=='CMS3')]
CMS4 <- DOWN$zval[which(DOWN$cms=='CMS4')]

pdf('BySubtype/DOWN.CMS_subtype.box.pdf')
boxplot(Normal, CMS1, CMS2, CMS3, CMS4,
names=c('Normal','CMS1','CMS2','CMS3','CMS4'),
col=c(brewer.pal(8, 'Dark2')[1], brewer.pal(8, 'Dark2')[2], brewer.pal(8, 'Dark2')[3], brewer.pal(8, 'Dark2')[4], brewer.pal(8, 'Dark2')[5]),
outline=F, las=2, notch=T, ylab='Zscore expression', main='DOWN')
abline(h=0, col='darkgrey', lty=2, lwd=1.5
)
dev.off()

# UP
Normal <- UP$zval[which(UP$cms=='Normal')]
CMS1 <- UP$zval[which(UP$cms=='CMS1')]
CMS2 <- UP$zval[which(UP$cms=='CMS2')]
CMS3 <- UP$zval[which(UP$cms=='CMS3')]
CMS4 <- UP$zval[which(UP$cms=='CMS4')]

pdf('BySubtype/UP.CMS_subtype.box.pdf')
boxplot(Normal, CMS1, CMS2, CMS3, CMS4,
names=c('Normal','CMS1','CMS2','CMS3','CMS4'),
col=c(brewer.pal(8, 'Dark2')[1], brewer.pal(8, 'Dark2')[2], brewer.pal(8, 'Dark2')[3], brewer.pal(8, 'Dark2')[4], brewer.pal(8, 'Dark2')[5]),
outline=F, las=2, notch=T, ylab='Zscore expression', main='UP')
abline(h=0, col='darkgrey', lty=2, lwd=1.5
)
dev.off()


a <- t.test(Normal, CMS1, paired=FALSE, var.equal=FALSE)
print(a$p.value)
a <- t.test(Normal, CMS2, paired=FALSE, var.equal=FALSE)
print(a$p.value)
a <- t.test(Normal, CMS3, paired=FALSE, var.equal=FALSE)
print(a$p.value)
a <- t.test(Normal, CMS4, paired=FALSE, var.equal=FALSE)
print(a$p.value)


##======================================================================================================== 
#
#p <- ggplot() +
#    geom_boxplot(aes(data1$cms, data1$zval, fill=data1$cms), outlier.shape=NA, notch=TRUE) + 
#    scale_fill_brewer(palette='Dark2') + 
#    ylim(-3, 3) + 
#    labs(x='', y='Zscore expression')+
#    theme_classic() +
#    theme(legend.position='none') +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#
#pdf('BySubtype/CMS_subtype.zval.box.v2.pdf')
#plot(p)
#dev.off()




