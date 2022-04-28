
library(rGREAT)


inputDIR <- '/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/PeakBED'


dir.create('GREAT_GOBP')
dir.create('GREAT_GOMF')
dir.create('GREAT_MousePheno')

dir.create('GREAT_GeneAnno')

dataTypeList <- c(
'hESC',
'Mesoderm',
'EarlyEC',
'MidEC',
'LateEC',
'FullEC',
'EarlyNonEC',
'MidNonEC',
'LateNonEC',
'FullNonEC'
)



for (dataType in dataTypeList){
    print(dataType)

    peakFile <- read.table(paste0(inputDIR, '/', dataType, '.bed'))
    peakList <- peakFile[, 1:3]
    greatJob <- submitGreatJob(peakList, species='hg19') 

    #availableOntologies(greatJob)
    valTable <- getEnrichmentTables(greatJob, download_by='tsv', ontology=c("GO Molecular Function", "GO Biological Process", "GO Cellular Component", "Mouse Phenotype" ))

    # GOBP
    allList <- valTable[['GO Biological Process']]
    selectedList <- allList[, c('Ontology','ID','Desc', 'BinomP', 'RegionFoldEnrich','GeneFoldEnrich','ObsGenes','TotalGenes','Regions','Genes')]
    selectedList$adjPval <- p.adjust(selectedList$BinomP, method='BH')
    write.table(selectedList, paste0('GREAT_GOBP/', dataType, '.GOBP.txt'), row.names=F, col.names=T, sep='\t', quote=F) 

    # GOMF
    allList <- valTable[['GO Molecular Function']]
    selectedList <- allList[, c('Ontology','ID','Desc', 'BinomP', 'RegionFoldEnrich','GeneFoldEnrich','ObsGenes','TotalGenes','Regions','Genes')]
    write.table(selectedList, paste0('GREAT_GOMF/', dataType, '.GOMF.txt'), row.names=F, col.names=T, sep='\t', quote=F) 

    # MousePheno
    allList <- valTable[['Mouse Phenotype']]
    selectedList <- allList[, c('Ontology','ID','Desc', 'BinomP', 'RegionFoldEnrich','GeneFoldEnrich','ObsGenes','TotalGenes','Regions','Genes')]
    write.table(selectedList, paste0('GREAT_MousePheno/', dataType, '.MousePheno.txt'), row.names=F, col.names=T, sep='\t', quote=F) 

    # Gene Annotation
    pdf(paste0('GREAT_GeneAnno/', dataType,'.GeneAnno.pdf'))
    genomeRanges <- plotRegionGeneAssociationGraphs(greatJob)
    dev.off()
    df <- data.frame(chrID=seqnames(genomeRanges), pt1=start(genomeRanges), pt2=end(genomeRanges), geneID=genomeRanges$gene, distToTSS=genomeRanges$distTSS)
    write.table(df, file=paste0('GREAT_GeneAnno/', dataType, '.GeneAnno.txt'), row.names=F, col.names=F, sep="\t", quote=F)

}



