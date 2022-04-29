#!/home/ajl1213/anaconda2/bin/python



import os
import numpy



DIR=os.getcwd()
StageSpecificGeneInfo='/home/ajl1213/Projects/Endothelial/data/RNA/Processing/DEGs/RNA.DegLabel.txt'
GeneExpressionInfo='/home/ajl1213/Projects/Endothelial/data/RNA/Processing/CountMatrix/RNA.TPM.qq.txt'
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

os.system('mkdir '+DIR+'/StageSpecificTF')



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


def LoadStageSpecificGenes():
    geneDict=getGeneInfo()

    DegDict={}
    input1=open(StageSpecificGeneInfo, 'r')
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


def getStageSpecificTfList():
    DegDict=LoadStageSpecificGenes()

    output1=open(DIR+'/StageSpecificTF/StageSpecificTfList.txt','w')
    output1.write('OrigCellType\tMotifID\tTfID\tPval\tTargetFrac\n')
    for cellType in cellTypeList:
        input1=open(DIR+'/RESULT/'+cellType+'.TfList.final.txt','r')
        all_input1=input1.readlines()
        for line in all_input1[1:]:
            each=line.strip().split('\t')
            motifID=each[0]
            pval=each[1]
            targetFrac=each[2]
            tfList=each[3]
	    for tf in tfList.split(','):
		if tf in DegDict[cellType]:
		    newLine=[cellType, motifID, tf, pval, targetFrac]
		    output1.write('\t'.join(newLine)+'\n')
	input1.close()
    output1.close()


def makePlotData():
    geneDict=getGeneInfo()

    ## Load P-values
    pvalDict={}
    for cellType in cellTypeList:
	input1=open(DIR+'/RESULT/'+cellType+'.TfList.final.txt','r')
	all_input1=input1.readlines()
	tmpDict={}
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    if float(each[1]) == 0:
		logPval=300
	    else:
		logPval=-numpy.log10(float(each[1]))

	    tfList=each[3]
	    for tf in tfList.split(','):
		tmpDict[tf]=logPval
	pvalDict[cellType]=tmpDict
	input1.close()


    ## Load RNA log2TPM values
    valDict={}
    input1=open(GeneExpressionInfo,'r')
    all_input1=input1.readlines()
    sampleList=all_input1[0].strip().split('\t')[2:]
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=geneDict[ensembleID]
	vals=numpy.array(each[2:], dtype='float')
	index=0
	tmpDict={}
	for val in vals:
	    tmpDict[sampleList[index]]=numpy.log2(val+1)
	    index+=1
	valDict[geneID]=tmpDict
    input1.close()

    finalDict={}
    for geneID in valDict:
	tmpDict={}
	for cellType in cellTypeList:
	    tmpList=[]
	    for sampleID in sampleList:
		if sampleID.split('_')[0] in cellTypeDict[cellType]:
		    val=valDict[geneID][sampleID]
		    tmpList.append(val)
	    meanVal=numpy.mean(tmpList)
	    tmpDict[cellType]=meanVal
	finalDict[geneID]=tmpDict


    ## Make data
    input1=open(DIR+'/StageSpecificTF/StageSpecificTfList.txt','r')
    output1=open(DIR+'/StageSpecificTF/StageSpecificTfList.ValMat.txt','w')
    output1.write('OrigCellType\tMotifID\tTfID\tCellType\tlogPval\tlog2RNA\n')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	origCellType=each[0]
	motifID=each[1]
	tfID=each[2]
	for cellType in cellTypeList:
	    if pvalDict[cellType].has_key(tfID):
		logPval=pvalDict[cellType][tfID]
	    else:
		logPval='0.0'

	    log2RNA=finalDict[tfID][cellType]

	    newLine=[origCellType, motifID, tfID, cellType, str(logPval), str(log2RNA)]
	    output1.write('\t'.join(newLine)+'\n')
    input1.close()
    output1.close()


    ## Make interval
    input1=open(DIR+'/StageSpecificTF/StageSpecificTfList.ValMat.txt','r')
    output1=open(DIR+'/StageSpecificTF/StageSpecificTfList.ValMat.final.txt','w')
    output1.write('OrigCellType\tMotifID\tTfID\tCellType\tlogPval\tlog2RNA\tPvalLabel\tRnaLabel\n')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[origCellType, motifID, tfID, cellType, logPval, log2RNA]=line.strip().split('\t')

	logPval=float(logPval)
	log2RNA=float(log2RNA)

	if (logPval >= 0 and logPval < 10):
	    PvalLabel='i0-10'
	if (logPval >= 10 and logPval < 30):
	    PvalLabel='i10-30'
	if (logPval >= 30 and logPval < 60):
	    PvalLabel='i30-60'
	if (logPval >= 60 and logPval < 100):
	    PvalLabel='i60-100'
	if logPval >= 100:
	    PvalLabel='above100'

	if (log2RNA >= 0 and log2RNA < 1):
	    RnaLabel='i0-1'
	if (log2RNA >= 1 and log2RNA < 3):
	    RnaLabel='i1-3'
	if (log2RNA >= 3 and log2RNA < 5):
	    RnaLabel='i3-5'
	if (log2RNA >= 5 and log2RNA < 7):
	    RnaLabel='i5-7'
	if log2RNA >= 7:
	    RnaLabel='above7'

	newLine=[origCellType, motifID, tfID, cellType, str(logPval), str(log2RNA), PvalLabel, RnaLabel]
	output1.write('\t'.join(newLine)+'\n')

    input1.close()
    output1.close()




getStageSpecificTfList()
makePlotData()



