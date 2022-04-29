#!/home/ajl1213/anaconda2/bin/python



import os
import numpy


DIR=os.getcwd()
inputDIR=DIR+'/DEGs'
outFile=DIR+'/DEGs/RNA.STable.txt'
cellTypeList=[
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
]


def makeSTable():
    # ordered gene list
    geneList=[]
    input1=open(DIR+'/DEGs/RNA.DegLabel.zscore.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneList.append(ensembleID)
    input1.close()

    #
    pvalDict={}
    qvalDict={}
    for cellType in cellTypeList:
	tmpDict1={}
	tmpDict2={}
	input1=open(inputDIR+'/'+cellType+'.valTable.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    logFC=each[1]
	    logCPM=each[2]
	    F=each[3]
	    pval=float(each[4])
	    qval=float(each[5])

	    tmpDict1[ensembleID]=pval
	    tmpDict2[ensembleID]=qval

	pvalDict[cellType]=tmpDict1
	qvalDict[cellType]=tmpDict2

	input1.close()

    #
    finalValDict={}
    for ensembleID in geneList:
	tmpList1=[]
	for cellType in cellTypeList:
	    if pvalDict[cellType].has_key(ensembleID):
		pval=pvalDict[cellType][ensembleID]
		tmpList1.append(pval)

	tmpList2=[]
	for cellType in cellTypeList:
	    if qvalDict[cellType].has_key(ensembleID):
		qval=qvalDict[cellType][ensembleID]
		tmpList2.append(qval)

	minPval=min(tmpList1)
	minQval=min(tmpList2)

	finalValDict[ensembleID]=[str(minPval), str(minQval)]

    #
    allDict={}
    input1=open(DIR+'/DEGs/RNA.DegLabel.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]
	maxStage=each[2]
	maxDataType=each[3]
	allDict[ensembleID]=[geneID, maxStage, maxDataType]
    input1.close()

    #
    output1=open(outFile,'w')
    output1.write('ensembleID\tgeneID\tmaxTimePoint\tDiffStage\tPval\tQval\n')
    for ensembleID in geneList:

	geneID, maxStage, maxDataType = allDict[ensembleID]
	pval, qval = finalValDict[ensembleID]

	newLine=[ensembleID, geneID, maxStage, maxDataType, pval, qval]	
	output1.write('\t'.join(newLine)+'\n')

    output1.close()

makeSTable()


