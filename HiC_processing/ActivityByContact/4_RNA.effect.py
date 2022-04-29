#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
resolution='10kb'
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
RnaCountFile='/home/ajl1213/Projects/Endothelial/data/RNA/Processing/CountMatrix/RNA.TPM.qq.txt'
AbcScoreDIR=DIR+'/ABC_score'
cellTypeList=['Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC']

os.system('mkdir '+DIR+'/RnaEffect')



def loadGeneInfo():
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


def loadRnaCount():
    valDict={}
    input1=open(RnaCountFile,'r')
    all_input1=input1.readlines()
    sampleList=all_input1[0].strip().split('\t')[2:]
    for line in all_input1[1:]:
        each=line.strip().split('\t')
	ensembleID=each[0]
        geneID=each[1]
        vals=numpy.array(each[2:], dtype='float')
        tmpDict={}
        index=0
        for val in vals:
            tmpDict[sampleList[index]]=val
            index+=1
	valDict[ensembleID]=tmpDict
    input1.close()

    #
    rnaCountDict={}
    for ensembleID in valDict:
	tmpDict={}
	for cellType in cellTypeList:
	    valList=[]

            for sampleID in sampleList:
                if sampleID.split('_')[0] in ['ME']:
                    dataType='Mesoderm'
                elif sampleID.split('_')[0] in ['EC01','EC02','EC04']:
                    dataType='EarlyEC'
                elif sampleID.split('_')[0] in ['EC06','EC08','EC12']:
                    dataType='MidEC'
                elif sampleID.split('_')[0] in ['EC24']:
                    dataType='LateEC'
                elif sampleID.split('_')[0] in ['EC48']:
                    dataType='FullEC'
                else:
                    dataType=''

                if dataType==cellType:
                    val=valDict[ensembleID][sampleID]
                    valList.append(val)
            meanVal=numpy.mean(valList)
            tmpDict[cellType]=meanVal

	rnaCountDict[ensembleID]=tmpDict

    return rnaCountDict


def loadInterVal(cellType):
    creSumDict={}
    interSumDict={}
    contactSumDict={}
    input1=open(AbcScoreDIR+'/'+cellType+'.abcScore.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')

	ensembleID=each[0]
	geneID=each[1]
	peakCount=float(each[6])
	interCount=float(each[7])
	contactVal=float(each[8])

	if creSumDict.has_key(ensembleID):
	    creSumDict[ensembleID]+=peakCount
	else:
	    creSumDict[ensembleID]=peakCount
	
	if interSumDict.has_key(ensembleID):
	    interSumDict[ensembleID]+=interCount
	else:
	    interSumDict[ensembleID]=interCount

	if contactSumDict.has_key(ensembleID):
	    contactSumDict[ensembleID]+=contactVal
	else:
	    contactSumDict[ensembleID]=contactVal	
    input1.close()

    return creSumDict, interSumDict, contactSumDict


def makeData():
    geneDict=loadGeneInfo()
    rnaCountDict=loadRnaCount()

    for cellType in cellTypeList:
	print cellType

	creSumDict, interSumDict, contactSumDict=loadInterVal(cellType)

	output1=open(DIR+'/RnaEffect/'+cellType+'.ContactValuePerGene.txt','w')
	output1.write('ensembleID\tgeneID\trnaCount\tinterSum\tcreSum\ttotalContactVal\n')
	for ensembleID in geneDict:
	    geneID=geneDict[ensembleID]

	    if rnaCountDict.has_key(ensembleID):
		rnaCount=str(rnaCountDict[ensembleID][cellType])
		if interSumDict.has_key(ensembleID):
		    interSum=str(interSumDict[ensembleID])
		    creSum=str(creSumDict[ensembleID])
		    contactSum=str(contactSumDict[ensembleID])

		    newLine=[ensembleID, geneID, rnaCount, interSum, creSum, contactSum]
		    output1.write('\t'.join(newLine)+'\n')
	output1.close()


makeData()



