
library(Vennerable)


#==============================================================================================================================
inputName <- paste0('VennInput.txt')
data1 <- read.table(inputName, header=T)

vennData <- Venn(SetNames = c('ME', 'EC12', 'EC48'), Weight = c(data1$nVal))

outName <- paste0('Plots/SigInter.venn.pdf')

pdf(outName)
#plot(vennData, doWeights = TRUE, type = "ChowRuskey",)
plot(vennData, doWeights = TRUE, type = "circles")
dev.off()


