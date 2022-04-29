library(preprocessCore)



data <- read.table("CountMatrix/H3K27ac.readCount.txt",head=T,row.names=1)

values <- as.matrix(data[,1:(ncol(data))])
qq_values <- normalize.quantiles.robust(values)
colnames(qq_values) <- colnames(data)[1:ncol(data)]
qq_values <- data.frame(peakID=row.names(data), qq_values)

write.table(qq_values, file='CountMatrix/H3K27ac.readCount.qq.txt', row.name=F, col.names=T, sep="\t", quote=F)






