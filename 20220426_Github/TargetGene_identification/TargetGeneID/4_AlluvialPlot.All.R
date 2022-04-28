
library(ggalluvial)
library(RColorBrewer)


dir.create('AlluvialPlots')

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


cellTypeList <- c('EarlyEC', 'MidEC', 'LateEC', 'FullEC')
colorList <- c(brewer.pal(8, 'Dark2')[2], brewer.pal(8, 'Dark2')[1], brewer.pal(8, 'Dark2')[5], brewer.pal(8, 'Dark2')[6])



# SNP to Peak
for (gwasID in gwasList){
    print(gwasID)

    #
    data1 <- read.table(paste0('RESULT/', gwasID, '.SnpToPeak.final.txt'), header=T)
    newDF <- NULL
    for (i in cellTypeList){
	tmp <- data1[which(data1$peakLabel==i),]
	tmp <- tmp[order(tmp$tagSNP),]

	newDF <- rbind(newDF, tmp)
    }

    newDF$tagSNP <- factor(newDF$tagSNP, levels=unique(newDF$tagSNP))
    newDF$peakID <- factor(newDF$peakID, levels=unique(newDF$peakID))
    newDF$peakLabel <- factor(newDF$peakLabel, levels=cellTypeList)

    #
    p1 <- ggplot(data = newDF, aes(axis1 = tagSNP, axis2 = peakID, y = frac)) +
	scale_x_discrete(limits = c('SNP','Peak'), expand = c(.2, .05)) +
	geom_alluvium(aes(fill=peakLabel)) +
	geom_stratum(aes(fill=peakLabel)) +
	scale_fill_manual(values=colorList) +
	#geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=2) +
	theme_minimal() 

    pdf(paste0('AlluvialPlots/', gwasID,'.SnpToPeak.pdf'))
    plot(p1)
    dev.off()

}


# Peak to Gene
for (gwasID in gwasList){
    print(gwasID)

    #
    data1 <- read.table(paste0('RESULT/', gwasID, '.PeakToGene.final.txt'), header=T)
    newDF <- NULL
    for (i in cellTypeList){
	tmp <- data1[which(data1$peakLabel==i),]
	tmp <- tmp[order(tmp$peakID),]

	newDF <- rbind(newDF, tmp)
    }

    newDF$peakID <- factor(newDF$peakID, levels=unique(newDF$peakID))
    newDF$geneID <- factor(newDF$geneID, levels=unique(newDF$geneID))
    newDF$peakLabel <- factor(newDF$peakLabel, levels=cellTypeList)

    #
    p1 <- ggplot(data = newDF, aes(axis1 = peakID, axis2 = geneID, y = frac)) +
	scale_x_discrete(limits = c('Peak','Gene'), expand = c(.2, .05)) +
	geom_alluvium(aes(fill=peakLabel)) +
	geom_stratum(aes(fill=peakLabel)) +
	scale_fill_manual(values=colorList) +
	#geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=2) +
	theme_minimal() 

    pdf(paste0('AlluvialPlots/', gwasID,'.PeakToGene.pdf'))
    plot(p1)
    dev.off()

}




