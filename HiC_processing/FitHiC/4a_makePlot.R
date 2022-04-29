library(RColorBrewer)


dir.create('Plots')
res <- '10kb'
qval <- '1e-1'


#=================================================================================================================================================
## SigInter Distance Histogram
sampleList <- c('ME','EC12','EC48')

for (i in sampleList){
    data1 <- read.table(paste0('SigInter/', i, '.', res, '.', qval, '.FitHiC.all.txt'), header=T)
    pdf(paste0('Plots/', i, '.SigInterDist.hist.pdf'), width=6, height=3, pointsize=3)
    hist(data1$dist, breaks=100, col='darkgrey', xlab='Genomic distance (bp)', ylab='Frequency', main=i, freq=F)
    lines(density(data1$dist), col='brown')

    abline(v=mean(data1$dist), lty=2, lwd=1.5)
    abline(v=median(data1$dist), lty=2, lwd=1.5)
    print(length(data1$dist))
    print(summary(data1$dist))
}

##
data1 <- read.table(paste0('SigInter/MergedEC.', res, '.', qval, '.FitHiC.ProCre.txt'), header=T)
pdf(paste0('Plots/MergedEC.SigInterDist.hist.pdf'), width=6, height=3, pointsize=3)
hist(data1$dist, breaks=100, col='darkgrey', xlab='Genomic distance (bp)', ylab='Frequency', main='MergedEC', ylim=c(0, 60000))
abline(v=mean(data1$dist), lty=2, lwd=1.5)
abline(v=median(data1$dist), lty=2, lwd=1.5)
print(length(data1$dist))
print(summary(data1$dist))


#=================================================================================================================================================
## SigInter Categorization PieChart
data1 <- read.table(paste0('SigInter/MergedEC.', res, '.', qval, '.FitHiC.ProCre.txt'), header=T)

df <- data.frame(table(data1$interLabel))
df$FracVal <- round(df$Freq/nrow(data1)*100, digits=2)
df

ccFrac <- df$FracVal[which(df$Var1=='CC')]
coFrac <- df$FracVal[which(df$Var1=='CO')]
ooFrac <- df$FracVal[which(df$Var1=='OO')]
poFrac <- df$FracVal[which(df$Var1=='PO')]
ppFrac <- df$FracVal[which(df$Var1=='PP')]
pcFrac <- df$FracVal[which(df$Var1=='PC')]

vals <- c(ccFrac, coFrac, ooFrac, poFrac, ppFrac, pcFrac)
label <- c('CC','CO','OO','PO','PP','PC')
label <- paste(label, '(', vals,'%)', sep='')

pdf('Plots/SigInter.pie.pdf', height=4, width=4, pointsize=3)
pie(vals, label, border='white')
dev.off()

## 3D pie
library(plotrix)

col_cc <- brewer.pal(9, 'Oranges')[6]
col_co <- brewer.pal(9, 'Oranges')[4]
col_oo <- brewer.pal(8, 'Set2')[8]
col_po <- brewer.pal(9, 'Blues')[4]
col_pp <- brewer.pal(9, 'Blues')[6]
col_pc <- brewer.pal(9, 'PuBuGn')[9]

pdf('Plots/SigInter.3dPie.pdf', height=4, width=4, pointsize=3)
pie3D(vals, labels=label, radius=1,  explode=0.1, col=c(col_cc, col_co, col_oo, col_po, col_pp, col_pc), main='Long-range interaction type')
dev.off()





