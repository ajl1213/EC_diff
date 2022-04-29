library(pheatmap)
library(RColorBrewer)



##
gwasList <- c(
'Crohns_disease',
'Asthma',
'Inflammatory_bowel_disease',
'Multiple_sclerosis',
'Eczema',
'Atrial_fibrillation',
'Atopic_asthma',
'Coronary_artery_disease',
'Varicose_veins',
'Alzheimers_disease_late_onset',
'Colorectal_cancer',
'Ulcerative_colitis',
'Rheumatoid_arthritis'
)

col_range <- colorRampPalette(c('white','yellow','darkorange','red','darkred'))(20)
celltype_colors <- c(brewer.pal(9, 'YlOrRd')[5], brewer.pal(9, 'PuBu')[5], brewer.pal(9, 'PuBu')[6], brewer.pal(9, 'PuBu')[7], brewer.pal(9, 'PuBu')[8]) 
maxVal <- 60

dir.create('Plots')


for (gwasID in gwasList){
    print(gwasID)

    #
    data1 <- read.table(paste0('ValMat/',gwasID,'.txt'), header=T, row.names=1)

    # row annotation
    annoRow <- data.frame(maxCell=data1$maxCell)
    rownames(annoRow) <- rownames(data1)

    setAnnoColor <- list(
    maxCell=c(
    Mesoderm=celltype_colors[1],
    EarlyEC=celltype_colors[2],
    MidEC=celltype_colors[3],
    LateEC=celltype_colors[4],
    FullEC=celltype_colors[5]
    )
    )

    #
    valMat <- data1[, 2:ncol(data1)]
    summary(as.vector(as.matrix(valMat)))


    pdf(paste0('Plots/', gwasID, '.heatmap.pdf'))
    pheatmap(valMat,
    color= col_range, 
    annotation_row=annoRow, annotation_colors=setAnnoColor,
    cluster_rows=F, cluster_cols=F,
    show_colnames=T, show_rownames=T,
    breaks=seq(0, maxVal, length.out=20)
    )
    dev.off()

}


##














