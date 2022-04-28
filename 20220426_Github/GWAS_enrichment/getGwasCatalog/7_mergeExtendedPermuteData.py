#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
cellType='MergedEC'
inFile=DIR+'/RESULT/'+cellType+'.DataList.txt'
LdSnpDIR=DIR+'/LdFinal'
fai='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
maxCopy=10



def selectGWAS(): 
    gwasList=[]
    input1=open(inFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[gwasID, nPeak, fcVal, empiricalPval]=line.strip().split('\t')
	if float(empiricalPval) == 0:
	    gwasList.append(gwasID)
    input1.close()

    return gwasList



gwasList=selectGWAS()

for gwasID in gwasList:

    output1=open(DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.RandPeakPermuteRecord.txt','w')
    output1.write('nPermute\tnPeak\n')

    for j in range(maxCopy):
	nCopy=str(j+1)

	input1=open(DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.'+nCopy+'.RandPeakPermuteRecord.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    nPermute=nCopy+'-'+each[0]
	    nPeak=each[1]
	    newLine=[nPermute, nPeak]
	    output1.write('\t'.join(newLine)+'\n')
	input1.close()
    output1.close()
	





