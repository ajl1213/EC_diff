#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
nGroup=100
cellTypeList=['Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC']


def combinePCC():
    corDict={}
    for idx in range(nGroup):
	idx=idx+1

	input1=open(DIR+'/PCC/PCC.'+str(idx)+'.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    geneID=each[1]
	    creID=each[2]
	    corVal=each[3]
	    if not corVal=='nan':
		corDict[ensembleID+';'+creID]=corVal
	input1.close()

    return corDict


def makeData():
    corDict=combinePCC()

    for cellType in cellTypeList:
	print cellType

	input1=open(DIR+'/ABC_score/'+cellType+'.abcScore.txt','r')
	output1=open(DIR+'/ABC_score/'+cellType+'.abcScore.cor.txt','w')
	output1.write('ensembleID\tgeneID\ttssBin\tcreID\tcreBin\tdist\tpeakCount\tinterCount\tcontactVal\tabcScore\tPCC\n')
	
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    [ensembleID, geneID, tssBin, creID, creBin, dist, peakCount, interCount, contactVal, abcScore]=line.strip().split('\t')
	    if corDict.has_key(ensembleID+';'+creID):
		corVal=corDict[ensembleID+';'+creID]
		newLine=[ensembleID, geneID, tssBin, creID, creBin, dist, peakCount, interCount, contactVal, abcScore, corVal]
		output1.write('\t'.join(newLine)+'\n')

	input1.close()
	output1.close()

	

makeData()



