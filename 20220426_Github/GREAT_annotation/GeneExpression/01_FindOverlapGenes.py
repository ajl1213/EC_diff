#!/home/ajl1213/anaconda2/bin/python



import os
import numpy
import scipy.stats as stats



DIR=os.getcwd()
DegInfo='/home/ajl1213/Projects/Endothelial/data/RNA/Processing/DEGs/RNA.DegLabel.txt'
AnnoGeneDIR='/home/ajl1213/Projects/Endothelial/Analysis2/AnnoGREAT/1_AnnoGREAT/GREAT_GeneAnno'
dataTypeList=['EC_cRE:hESC-Mesoderm-EarlyEC-MidEC-LateEC-FullEC', 'nonEC_cRE:hESC-Mesoderm-EarlyNonEC-MidNonEC-LateNonEC-FullNonEC']
cellTypeList=['hESC','Mesoderm','EarlyEC','MidEC','LateEC','FullEC','EarlyNonEC','MidNonEC','LateNonEC','FullNonEC']
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'

os.system('mkdir '+DIR+'/OverlapGeneList')
os.system('mkdir '+DIR+'/CountMatrix')



def getGeneInfo():
    geneDict={}
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]
	geneDict[ensembleID]=geneID
    input1.close()

    return geneDict


def LoadDEGs():
    geneDict=getGeneInfo()

    DegDict={}
    input1=open(DegInfo,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=geneDict[ensembleID]
	stageType=each[3]
	if DegDict.has_key(stageType):
	    tmp=DegDict[stageType]
	    tmp.append(geneID)
	    DegDict[stageType]=tmp
	else:
	    DegDict[stageType]=[geneID]
    input1.close()

    return DegDict


def LoadAnnoGenes():
    annoGeneDict={}
    for cellType in cellTypeList:
	input1=open(AnnoGeneDIR+'/'+cellType+'.geneList.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    geneID=each[0]
	    if annoGeneDict.has_key(cellType):
		tmp=annoGeneDict[cellType]
		tmp.append(geneID)
		annoGeneDict[cellType]=tmp
	    else:
		annoGeneDict[cellType]=[geneID]
	input1.close()

    return annoGeneDict


def findOverlapGenes():
    DegDict=LoadDEGs()
    annoGeneDict=LoadAnnoGenes()

    for cellType in cellTypeList:
	output1=open(DIR+'/OverlapGeneList/'+cellType+'.overlapGeneList.txt','w')
	overlapGeneList=list(set(DegDict[cellType]) & set(annoGeneDict[cellType]))
	output1.write('\n'.join(overlapGeneList))
	output1.close()



def makeCountMatrix():
    geneDict=getGeneInfo()
    DegDict=LoadDEGs()
    annoGeneDict=LoadAnnoGenes()

    for i in dataTypeList:
	dataType=i.split(':')[0]
	cellTypeList=i.split(':')[1].split('-')

	## Count matrix of overlap genes
	output1=open(DIR+'/CountMatrix/'+dataType+'.overlapCount.txt','w')
	output1.write('CellType.PeakAnno\t'+'\t'.join(cellTypeList)+'\n')
	for PeakAnnoCellType in cellTypeList:
	    tmpList=[]
	    for DegCellType in cellTypeList:
		nOverlapGene=str(len(list(set(annoGeneDict[PeakAnnoCellType]) & set(DegDict[DegCellType]))))
		tmpList.append(nOverlapGene)

	    output1.write(PeakAnnoCellType+'\t'+'\t'.join(tmpList)+'\n')
	output1.close()


	## P-value matrix (Fisher's exact test of overlap genes)
	# get nTotalGene
	nGene={}
	for ensembleID in geneDict:
	    geneID=geneDict[ensembleID]
	    nGene[geneID]=1
	nTotalGene= len(nGene)

	# get P-val matrix	
	output1=open(DIR+'/CountMatrix/'+dataType+'.overlapPval.txt','w')
	output1.write('CellType.PeakAnno\t'+'\t'.join(cellTypeList)+'\n')
	for PeakAnnoCellType in cellTypeList:
	    tmpList=[]
	    for DegCellType in cellTypeList:

		nAnnoGene=len(set(annoGeneDict[PeakAnnoCellType]))
		nDEG=len(set(DegDict[DegCellType]))
		nOverlapGene=len(list(set(annoGeneDict[PeakAnnoCellType]) & set(DegDict[DegCellType])))

		oddsratio, pval = stats.fisher_exact([[nTotalGene - nDEG, nDEG], [nAnnoGene-nOverlapGene, nOverlapGene]], alternative='greater')
		logPval=str(-numpy.log10(pval))

		tmpList.append(logPval)

	    output1.write(PeakAnnoCellType+'\t'+'\t'.join(tmpList)+'\n')
	output1.close()




findOverlapGenes()
makeCountMatrix()




