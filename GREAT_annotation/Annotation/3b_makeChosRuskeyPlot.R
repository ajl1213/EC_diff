
library(RColorBrewer)
library(Vennerable)


dir.create('ChowRuskeyPlot')


#======================================== Chow-Ruskey Plot
for (i in c('EC_cRE','nonEC_cRE')){

    inputName <- paste0('ChowInput/',i, '.VennInput.txt')
    data1 <- read.table(inputName, header=T)

    ##
    vennData <- Venn(SetNames = c('Early', 'Mid', 'Late', 'Full'), Weight = c(data1$nVal))

    outName <- paste0('ChowRuskeyPlot/', i, '.ChowRuskey.pdf')

    pdf(outName)
    plot(vennData, doWeights = TRUE, type = "ChowRuskey",)
    dev.off()

}





