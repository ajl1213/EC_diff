#!/home/ajl1213/anaconda2/bin/python



import os
import numpy



DIR=os.getcwd()
inputDIR='/home/ajl1213/Projects/Endothelial/data/RNA/Mapping'
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
sampleList=[
'hESC_rep1',
'hESC_rep2',
'ME_rep1',
'ME_rep2',
'EC01_rep1',
'EC01_rep2',
'EC02_rep1',
'EC02_rep2',
'EC04_rep1',
'EC04_rep2',
'EC06_rep1',
'EC06_rep2',
'EC08_rep1',
'EC08_rep2',
'EC12_rep1',
'EC12_rep2',
'EC24_rep1',
'EC24_rep2',
'EC48_rep1',
'EC48_rep2',
'nonEC01_rep1',
'nonEC01_rep2',
'nonEC02_rep1',
'nonEC02_rep2',
'nonEC04_rep1',
'nonEC04_rep2',
'nonEC06_rep1',
'nonEC06_rep2',
'nonEC08_rep1',
'nonEC08_rep2',
'nonEC12_rep1',
'nonEC12_rep2',
'nonEC24_rep1',
'nonEC24_rep2',
'nonEC48_rep1',
'nonEC48_rep2'
]

print len(sampleList)

os.system('mkdir '+DIR+'/CountMatrix')



def makeDataInfo():
    output1=open(DIR+'/DataInfo.txt','w')
    output1.write('sampleID\tstageLabel\tdataType\tcellType\n')
    for sampleID in sampleList:
        stageLabel=sampleID.split('_')[0]

        if stageLabel in ['hESC']:
            dataType='hESC'
            cellType='hESC'
        if stageLabel in ['ME']:
            dataType='Mesoderm'
            cellType='Mesoderm'
        if stageLabel in ['EC01','EC02','EC04']:
            dataType='EarlyEC'
            cellType='EC'
        if stageLabel in ['EC06','EC08','EC12']:
            dataType='MidEC'
            cellType='EC'
        if stageLabel in ['EC24']:
            dataType='LateEC'
            cellType='EC'
        if stageLabel in ['EC48']:
            dataType='FullEC'
            cellType='EC'
        if stageLabel in ['nonEC01','nonEC02','nonEC04']:
            dataType='EarlyNonEC'
            cellType='nonEC'
        if stageLabel in ['nonEC06','nonEC08','nonEC12']:
            dataType='MidNonEC'
            cellType='nonEC'
        if stageLabel in ['nonEC24']:
            dataType='LateNonEC'
            cellType='nonEC'
        if stageLabel in ['nonEC48']:
            dataType='FullNonEC'
            cellType='nonEC'

        newLine=[sampleID, stageLabel, dataType, cellType]
        output1.write('\t'.join(newLine)+'\n')
    output1.close()


def getGeneSet():
    GeneDict={}
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]
	GeneDict[ensembleID]=geneID
    input1.close()

    return GeneDict


def getStat():
    output1=open(DIR+'/MapStat.txt','w')
    output1.write('sampleID\tReadTotal\tUniqMapPercent\tUniqMapRead\n')
    for sampleID in sampleList:
	input1=open(inputDIR+'/AlignedBam/'+sampleID+'/Log.final.out','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    if line.count('Number of input reads')>0:
		inputReadNum=each[1]
	    if line.count('Uniquely mapped reads %')>0:
		mapPercent=each[1]
	    if line.count('Uniquely mapped reads number')>0:
		uniqReadNum=each[1]
	output1.write(sampleID+'\t'+inputReadNum+'\t'+mapPercent+'\t'+uniqReadNum+'\n')
	input1.close()
    output1.close()


def getCountInfo():
    GeneDict=getGeneSet()

    ValDict={}
    for sampleID in sampleList:
	tmpDict={}
	input1=open(inputDIR+'/CountInfo/'+sampleID+'.genes.results','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0].split('_')[0]
	    transcriptID=each[1]
	    count_val=each[4]
	    tpm_val=each[5]
	    fpkm_val=each[6]

	    tmpDict[ensembleID]=count_val
	input1.close()

	ValDict[sampleID]=tmpDict

    return ValDict


def makeMatrix():
    GeneDict=getGeneSet()
    ValDict=getCountInfo()

    ## rawCountMatrix
    output1=open(DIR+'/CountMatrix/RNA.rawCount.txt','w')
    output1.write('ensembleID\tgeneID\t'+'\t'.join(sampleList)+'\n')
    for ensembleID in GeneDict:
	geneID=GeneDict[ensembleID]
	output1.write(ensembleID+'\t'+geneID)
	for sampleID in sampleList:
	    count_val=ValDict[sampleID][ensembleID]
	    output1.write('\t'+count_val)
	output1.write('\n')
    output1.close()




makeDataInfo()
getStat()
makeMatrix()


