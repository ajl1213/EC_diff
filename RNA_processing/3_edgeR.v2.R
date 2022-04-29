
library(edgeR)



# Load Data
data1 <- read.table("CountMatrix/RNA.rawCount.qq.txt",header=T, row.names=1)
valMat <- data1[, 2:ncol(data1)]

metadata <- read.table('DataInfo.txt', header=T)
dataType <- factor(metadata$dataType)

# Create a DGEList object
y <- DGEList(counts=valMat, group=dataType)
y

# TPM
rpm_df <- data.frame(ensembleID=rownames(cpm(y)), geneID= data1$geneID, cpm(y))
write.table(rpm_df, file='CountMatrix/RNA.TPM.qq.txt', row.names=F, col.names=T, sep="\t", quote=F)


#============================================================================================================================
dir.create('DEGs')

# Differentiated cells
for (i in c('hESC','Mesoderm','FullEC','FullNonEC')){
    print(i)

    testVec <- c()
    testVec[dataType==i] <- 1
    testVec[dataType!=i] <- 0

    # Design matrix
    design <- model.matrix(~testVec)
    rownames(design) <- colnames(y)
    design

    # Estimate dispersion
    y_nd <- estimateDisp(y, design, robust=TRUE)
    y_nd$common.dispersion

    # GLM fitting
    fit <- glmQLFit(y_nd, design, robust=TRUE)
    head(fit$coefficients)

    # Identify differential peaks
    qlf<- glmQLFTest(fit)
    qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
    valTable <- qlf$table[which(qlf$table$FDR < 0.05),]
    valTable <- data.frame(ensembleID=rownames(valTable), valTable)
    write.table(valTable, file=paste0('DEGs/', i,'.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)
}

#============================================================================================================================
# "EarlyEC" "EarlyNonEC" "FullEC" "FullNonEC" "hESC" "LateEC" "LateNonEC" "Mesoderm" "MidEC" "MidNonEC"  

# Differentiating cells
design <- model.matrix(~0+dataType, data=y$samples)
colnames(design) <- levels(dataType)

# EarlyEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(1,0,0,0,0,0,0,-1,0,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DEGs/EarlyEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

# EarlyNonEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,1,0,0,0,0,0,-1,0,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DEGs/EarlyNonEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

# MidEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,0,0,0,0,0,0,-1,1,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DEGs/MidEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

# MidNonEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,0,0,0,0,0,0,-1,0,1))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DEGs/MidNonEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

# LateEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,0,0,0,0,1,0,-1,0,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DEGs/LateEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

# LateNonEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,0,0,0,0,0,1,-1,0,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DEGs/LateNonEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)





