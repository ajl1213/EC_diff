#!/home/ajl1213/anaconda2/bin/python



import os
import numpy



DIR=os.getcwd()
DataInfo=DIR+'/DataInfo.txt'
inputMatrix=DIR+'/CountMatrix/RNA.TPM.qq.txt'
stageList=['hESC', 'ME', 'EC01', 'EC02', 'EC04', 'EC06', 'EC08', 'EC12', 'EC24', 'EC48', 'nonEC01', 'nonEC02', 'nonEC04', 'nonEC06', 'nonEC08', 'nonEC12', 'nonEC24', 'nonEC48']
dataTypeList=['hESC', 'Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC', 'EarlyNonEC', 'MidNonEC', 'LateNonEC', 'FullNonEC']
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'


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


def loadDataInfo():
    dataDict={}
    input1=open(DataInfo,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	sampleID=each[0]
	stageLabel=each[1]
	dataType=each[2]
	dataDict[sampleID]=[stageLabel, dataType]
    input1.close()

    return dataDict


def loadCountData():
    valDict={}
    input1=open(inputMatrix,'r')
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

    return sampleList, valDict


def loadDEGs():
    degList=[]
    for dataType in dataTypeList:
	input1=open(DIR+'/DEGs/'+dataType+'.valTable.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    degList.append(ensembleID)
	input1.close()

    degList=set(degList)

    return degList


def makeLabel():
    geneDict=getGeneInfo()
    dataDict = loadDataInfo()
    sampleList, valDict= loadCountData()
    degList = loadDEGs()

    output1=open(DIR+'/DEGs/RNA.DegLabel.txt', 'w')
    output1.write('ensembleID\tgeneID\tmaxStage\tmaxDataType\n')
    for ensembleID in degList:
	geneID=geneDict[ensembleID]
	tmpDict={}
	for stageLabel in stageList:
	    tmpList=[]
	    for sampleID in sampleList:
		if dataDict[sampleID][0]==stageLabel:
		    val=valDict[ensembleID][sampleID]
		    tmpList.append(val)
	    meanVal=numpy.mean(tmpList)
	    tmpDict[stageLabel]=meanVal

	maxStage=max(tmpDict, key=tmpDict.get)

        if maxStage in ['hESC']:
            maxDataType='hESC'
        if maxStage in ['ME']:
            maxDataType='Mesoderm'
        if maxStage in ['EC01','EC02','EC04']:
            maxDataType='EarlyEC'
        if maxStage in ['EC06','EC08','EC12']:
            maxDataType='MidEC'
        if maxStage in ['EC24']:
            maxDataType='LateEC'
        if maxStage in ['EC48']:
            maxDataType='FullEC'
        if maxStage in ['nonEC01','nonEC02','nonEC04']:
            maxDataType='EarlyNonEC'
        if maxStage in ['nonEC06','nonEC08','nonEC12']:
            maxDataType='MidNonEC'
        if maxStage in ['nonEC24']:
            maxDataType='LateNonEC'
        if maxStage in ['nonEC48']:
            maxDataType='FullNonEC'

	newLine=[ensembleID, geneID, maxStage, maxDataType]
	output1.write('\t'.join(newLine)+'\n')

    output1.close()


def makeAllGeneMatrix():
    geneDict=getGeneInfo()
    dataDict = loadDataInfo()
    sampleList, valDict= loadCountData()

    meanDict={}
    stdDict={}
    for ensembleID in valDict:
	tmpList=[]
	for sampleID in sampleList:
	    val=valDict[ensembleID][sampleID]
	    tmpList.append(val)
	meanVal=numpy.mean(tmpList)
	stdVal=numpy.std(tmpList)

	meanDict[ensembleID]=meanVal
	stdDict[ensembleID]=stdVal

    ## Make Zscore Matrix
    output1=open(DIR+'/CountMatrix/RNA.TPM.qq.zscore.txt', 'w')
    output1.write('ensembleID\tgeneID\t'+'\t'.join(sampleList)+'\n')
    for ensembleID in valDict:
	geneID=geneDict[ensembleID]
	output1.write(ensembleID+'\t'+geneID)
	for sampleID in sampleList:
	    val=valDict[ensembleID][sampleID]
	    zscore=(val-meanDict[ensembleID])/stdDict[ensembleID] if not stdDict[ensembleID]==0 else 0
	    output1.write('\t'+str(zscore))
	output1.write('\n')
    output1.close()

    ## Merged TPM Matrix
    output1=open(DIR+'/CountMatrix/RNA.TPM.qq.stage.txt','w')
    output1.write('ensembleID\tgeneID\t'+'\t'.join(stageList)+'\n')
    for ensembleID in valDict:
	geneID=geneDict[ensembleID]
	meanList=[]
	for stageLabel in stageList:
	    tmpList=[]
	    for sampleID in sampleList:
		if dataDict[sampleID][0]==stageLabel:
		    val=valDict[ensembleID][sampleID]
		    tmpList.append(val)
	    meanVal=numpy.mean(tmpList)
	    meanList.append(str(meanVal))
	output1.write(ensembleID+'\t'+geneID+'\t'+'\t'.join(meanList)+'\n')
    output1.close()


def makeDiffMatrix():
    geneDict=getGeneInfo()
    dataDict = loadDataInfo()
    sampleList, valDict= loadCountData()

    ## Load DEGs
    degDict={}
    input1=open(DIR+'/DEGs/RNA.DegLabel.txt', 'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]
	maxStage=each[2]
	if degDict.has_key(maxStage):
	    tmp=degDict[maxStage]
	    tmp.append(ensembleID)
	    degDict[maxStage]=tmp
	else:
	    degDict[maxStage]=[ensembleID]
    input1.close()

    ## Make Zscore Matrix
    output1=open(DIR+'/DEGs/RNA.DegLabel.zscore.txt', 'w')
    output1.write('ensembleID\tgeneID\tmaxStage\t'+'\t'.join(stageList)+'\n')
    for i in stageList:
	for ensembleID in degDict[i]:
	    geneID=geneDict[ensembleID]
	    tmpDict={}
	    for stageLabel in stageList:
		tmpList=[]
		for sampleID in sampleList:
		    if dataDict[sampleID][0]==stageLabel:
			val=valDict[ensembleID][sampleID]
			tmpList.append(val)
		meanVal=numpy.mean(tmpList)
		tmpDict[stageLabel]=meanVal


	    newLine=[ensembleID, geneID, i]

	    valList=[]
	    for stageLabel in stageList:
		val=tmpDict[stageLabel]
		valList.append(val)

	    meanVal=numpy.mean(valList)
	    stdVal=numpy.std(valList)

	    for stageLabel in stageList:
		val=tmpDict[stageLabel]
		zscore=(val-meanVal)/stdVal if not stdVal==0 else 0
		newLine.append(str(zscore))
		
	    output1.write('\t'.join(newLine)+'\n')

    output1.close()




makeLabel()
makeAllGeneMatrix()
makeDiffMatrix()



