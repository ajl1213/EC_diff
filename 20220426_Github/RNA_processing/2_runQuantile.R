library(preprocessCore)


#
data <- read.table("CountMatrix/RNA.rawCount.txt",head=T,row.names=1)
valMat <- as.matrix(data[,2:(ncol(data))])

qq_values <- normalize.quantiles.robust(valMat)

#
colnames(qq_values) <- colnames(data)[2:ncol(data)]
qq_values <- data.frame(ensembleID=rownames(data), geneID=data$geneID, qq_values)


write.table(qq_values, file='CountMatrix/RNA.rawCount.qq.txt', row.name=F, col.names=T, sep="\t", quote=F)








