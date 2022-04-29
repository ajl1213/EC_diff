
library(edgeR)



# Load Data
data1 <- read.table("CountMatrix/H3K27ac.readCount.qq.txt",header=T, row.names=1)
valMat <- data1

metadata <- read.table('DataInfo.txt', header=T)
dataType <- factor(metadata$dataType)

# Create a DGEList object
y <- DGEList(counts=valMat, group=dataType)
y

##
# RPM
rpm_df <- data.frame(peakID=rownames(cpm(y)), cpm(y))
write.table(rpm_df, file='CountMatrix/H3K27ac.RPM.qq.txt', row.names=F, col.names=T, sep="\t", quote=F)


#============================================================================================================================
dir.create('DiffPeaks/')

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
    valTable <- data.frame(peakID=rownames(valTable), valTable)
    write.table(valTable, file=paste0('DiffPeaks/', i,'.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)
}

#============================================================================================================================
# "EarlyEC" "EarlyNonEC" "FullEC" "FullNonEC" "hESC" "LateEC" "LateNonEC" "Mesoderm" "MidEC" "MidNonEC"  

# Differentiating cells
design <- model.matrix(~0+dataType, data=y$samples)
colnames(design) <- levels(dataType)

## EarlyEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(1,0,0,0,0,0,0,-1,0,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DiffPeaks/EarlyEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

## EarlyNonEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,1,0,0,0,0,0,-1,0,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DiffPeaks/EarlyNonEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

## MidEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,0,0,0,0,0,0,-1,1,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DiffPeaks/MidEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

## MidNonEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,0,0,0,0,0,0,-1,0,1))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DiffPeaks/MidNonEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

## LateEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,0,0,0,0,1,0,-1,0,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DiffPeaks/LateEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)

## LateNonEC
y_nd <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y_nd, design, robust=TRUE)
qlf<- glmQLFTest(fit, contrast=c(0,0,0,0,0,0,1,-1,0,0))
qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")
valTable=qlf$table[which(qlf$table$FDR < 0.05),]
valTable=data.frame(peakID=rownames(valTable), valTable)
write.table(valTable, file=paste0('DiffPeaks/LateNonEC.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)






