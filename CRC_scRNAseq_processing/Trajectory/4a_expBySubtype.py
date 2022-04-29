#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import gzip
import operator


DIR=os.getcwd()
subTypeInfo='/home/ajl1213/Projects/Endothelial/data/snRNA/ColonCancer_GSE132465/1_downLoad'
valMat=DIR+'/ValMat/zmat.txt'


os.system('mkdir '+DIR+'/BySubtype')

subtypeDict={
'SMC06-T':'CMS1',
'SMC08-T':'CMS1',
'SMC15-T':'CMS1',
'SMC10-T':'CMS1',
'SMC03-T':'CMS1',
'SMC22-T':'CMS2',
'SMC18-T':'CMS2',
'SMC21-T':'CMS2',
'SMC09-T':'CMS2',
'SMC23-T':'CMS2',
'SMC25-T':'CMS2',
'SMC11-T':'CMS2',
'SMC07-T':'CMS2',
'SMC16-T':'CMS3',
'SMC01-T':'CMS3',
'SMC19-T':'CMS3',
'SMC05-T':'CMS3',
'SMC20-T':'CMS4',
'SMC14-T':'CMS4',
'SMC17-T':'CMS4',
'SMC04-T':'CMS4',
'SMC02-T':'CMS4',
'SMC24-T':'CMS4',
'SMC09-N':'Normal',
'SMC04-N':'Normal',
'SMC02-N':'Normal',
'SMC05-N':'Normal',
'SMC01-N':'Normal',
'SMC10-N':'Normal',
'SMC07-N':'Normal',
'SMC03-N':'Normal',
'SMC06-N':'Normal',
'SMC08-N':'Normal'
}

def loadValue():
    valDict={}
    input1=open(valMat,'r')
    all_input1=input1.readlines()
    cellList=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	vals=numpy.array(each[1:], dtype='float')
	tmpDict={}
	idx=0
	for val in vals:
	    tmpDict[cellList[idx]]=val
	    idx+=1

	valDict[geneID]=tmpDict
    input1.close()
 
    return cellList, valDict


def loadGeneInfo():

    geneLabelDict={}
    input1=open(DIR+'/ValMat/targetGene.log2fc.order.txt','r')
    all_input1=input1.readlines()
    barcodeList=all_input1[0].strip().split('\t')[5:]
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	log2fc=float(each[4])
	if log2fc > 0:
	    geneLabel='up'
	else:
	    geneLabel='down'
	geneLabelDict[geneID]=geneLabel

    input1.close()

    return barcodeList, geneLabelDict


def makeData():
    cellList, valDict=loadValue()
    barcodeList, geneLabelDict=loadGeneInfo()

    output1=open(DIR+'/BySubtype/zval.bysubtype.txt','w')
    output1.write('geneID\tgeneLabel\tbarcodeID\tcms\tzval\n')
    for geneID in valDict:
	geneLabel=geneLabelDict[geneID]

	for cell in cellList:
	    if cell in barcodeList:
		barcodeID=cell.replace('.','-')
		cmsID=subtypeDict[cell.split('_')[0].replace('.','-')]
		val=valDict[geneID][cell]

		newLine=[geneID, geneLabel, barcodeID, cmsID, str(val)] 

		output1.write('\t'.join(newLine)+'\n')

    output1.close()
    



makeData()



