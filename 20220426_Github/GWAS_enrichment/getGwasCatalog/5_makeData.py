#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import random



DIR=os.getcwd()
LdSnpDIR=DIR+'/LdFinal'
PeakDIR=DIR+'/InputPeakBED'
cellType='MergedEC'


gwasList=[]
input1=open(DIR+'/GwasList.txt','r')
all_input1=input1.readlines()
for line in all_input1[1:]:
    each=line.strip().split('\t')
    gwasID=each[0]
    gwasList.append(gwasID)
input1.close()

os.system('mkdir '+DIR+'/ObservedData')
os.system('mkdir '+DIR+'/RESULT')



def getActualData():

    for gwasID in gwasList:

	runScript='intersectBed -a '+PeakDIR+'/'+cellType+'.bed -b '+LdSnpDIR+'/'+gwasID+'.bed -wa -wb | cut -f1,2,3 | sort -u | wc -l > '+DIR+'/ObservedData/'+gwasID+'.'+cellType+'.nPeakIntersect.txt'
	print runScript
	os.system(runScript)

	runScript='intersectBed -a '+PeakDIR+'/'+cellType+'.bed -b '+LdSnpDIR+'/'+gwasID+'.bed -wa -wb | cut -f1,2,3 | sort -u > '+DIR+'/ObservedData/'+gwasID+'.'+cellType+'.Peak.intersect.bed'
	print runScript
	os.system(runScript)

	runScript='intersectBed -a '+PeakDIR+'/'+cellType+'.bed -b '+LdSnpDIR+'/'+gwasID+'.bed -wa -wb | sort -u > '+DIR+'/ObservedData/'+gwasID+'.'+cellType+'.Peak.intersect.full.bed'
	print runScript
	os.system(runScript)


def makeData():

    output1=open(DIR+'/RESULT/'+cellType+'.FcOverExpectation.MergedRecord.txt','w')
    output1.write('gwasID\tdataLabel\tfcVal\n')
    output2=open(DIR+'/RESULT/'+cellType+'.DataList.txt','w')
    output2.write('gwasID\tnPeak\tfcVal\tempiricalPval\n')

    for gwasID in gwasList:

	## Load observed values
	input1=open(DIR+'/ObservedData/'+gwasID+'.'+cellType+'.nPeakIntersect.txt','r')
	all_input1=input1.readlines()
	nPeak=all_input1[0].strip().split('\t')[0]
	observedVal=float(nPeak)

	## Load random values
	randValList=[]
	input1=open(DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.RandPeakPermuteRecord.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    nPeak=float(each[1])
	    randValList.append(nPeak)
	input1.close()

	## Compute empirical pval
	tmpList=[]
	for randVal in randValList:
	    if randVal > observedVal:
		tmpList.append(randVal)

	pval=float(len(tmpList))/float(len(randValList))

	## Fold change (over expectation)
	tmpList=[]
	for randVal in randValList:
	    tmpList.append(randVal)
	meanVal=numpy.mean(tmpList)
	
	for randVal in randValList:
	    fcVal=float(randVal) / float(meanVal) if not meanVal == 0 else str(0)
	    newLine=[gwasID, 'RandomPermute', str(fcVal)]
	    output1.write('\t'.join(newLine)+'\n')

	fcVal=float(observedVal) / float(meanVal) if not meanVal == 0 else str(0)

	newLine=[gwasID, 'Observed', str(fcVal)]
	output1.write('\t'.join(newLine)+'\n')

	##
	newLine=[gwasID, str(observedVal), str(fcVal), str(pval)]
	output2.write('\t'.join(newLine)+'\n')

    output1.close()
    output2.close()





getActualData()
makeData()



