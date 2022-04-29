#!/home/ajl1213/anaconda2/bin/python



import os
import numpy
import operator



DIR=os.getcwd()
inputDIR='/home/ajl1213/Projects/Endothelial/Analysis2/AnnoGREAT/1_AnnoGREAT/GREAT_GeneAnno'
overlapGeneDIR=DIR+'/OverlapGeneList'
humanTfCatalog='/home/ajl1213/genome.info/HumanTF_catalog/DatabaseExtract_v_1.01.txt'
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
cellTypeList=['hESC', 'Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC', 'EarlyNonEC', 'MidNonEC', 'LateNonEC', 'FullNonEC']
cellTypeDict={
'hESC':['hESC'],
'Mesoderm':['ME'],
'EarlyEC':['EC01','EC02','EC04'],
'MidEC':['EC06','EC08','EC12'],
'LateEC':['EC24'],
'FullEC':['EC48'],
'EarlyNonEC':['nonEC01','nonEC02','nonEC04'],
'MidNonEC':['nonEC06','nonEC08','nonEC12'],
'LateNonEC':['nonEC24'],
'FullNonEC':['nonEC48']
}

distCutoff= 100000



def getGeneInfo():
    geneDict={}
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        geneDict[geneID]=ensembleID
    input1.close()

    return geneDict


def loadPeakInfo(cellType):
    geneDict=getGeneInfo()

    nPeakDict={}
    PeakDict={}
    input1=open(inputDIR+'/'+cellType+'.geneList.txt', 'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
	geneID=each[0]
	peakList=each[1]

	newPeakList=[]
	for i in peakList.split(','):
	    if geneDict.has_key(geneID):
		peakID=i.split(';')[0]
		peakDist=i.split(';')[1]

		if abs(int(peakDist)) < distCutoff:
		    newPeakList.append(peakID+'('+peakDist+')')

	nPeak=len(newPeakList)

	if nPeak > 0:
	    nPeakDict[geneID]=nPeak
	    PeakDict[geneID]=newPeakList

        input1.close()

    return nPeakDict, PeakDict


def loadTfCatalog():
    tfDict={}
    input1=open(humanTfCatalog,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
	ensembleID=each[1]
        geneID=each[2]
        query=each[4]
        if query=='Yes':
            tfDict[geneID]=1
    input1.close()

    return tfDict


def loadOverlapGenes():
    overlapGeneDict={}
    for cellType in cellTypeList:
        input1=open(overlapGeneDIR+'/'+cellType+'.overlapGeneList.txt','r')
        all_input1=input1.readlines()
        geneList=[]
        for line in all_input1:
            each=line.strip().split('\t')
            geneID=each[0]
            geneList.append(geneID)
            overlapGeneDict[cellType]=geneList

        input1.close()

    return overlapGeneDict


def makeData():
    tfDict=loadTfCatalog()
    overlapGeneDict=loadOverlapGenes()

    output1=open(DIR+'/CandidateList.txt','w')
    output1.write('DEG_type\tGene_ID\tTF_Label\tnPeak\tAnnotatedPeak_GREAT\n')
    for cellType in cellTypeList:
	if not cellType in ['hESC', 'Mesoderm', 'LateEC', 'FullEC', 'EarlyNonEC', 'MidNonEC', 'LateNonEC', 'FullNonEC']:
	    nPeakDict, PeakDict = loadPeakInfo(cellType)

	    for key in sorted(nPeakDict.items(), key=operator.itemgetter(1), reverse=True):
		geneID=key[0]
		if geneID in overlapGeneDict[cellType]:
		    if tfDict.has_key(geneID):
			tfLabel='humanTF'
		    else:
			tfLabel='non-TF'
    
		    nPeak=key[1]
		    peakList=PeakDict[geneID]
		    newLine=[cellType, geneID, tfLabel, str(nPeak), ', '.join(peakList)]

		    output1.write('\t'.join(newLine)+'\n')

    output1.close()




makeData()
 


