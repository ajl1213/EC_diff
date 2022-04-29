#!/home/ajl1213/anaconda2/bin/python


import os
import sys
import numpy
import random



DIR=os.getcwd()
gwasID=sys.argv[1]
cellType=sys.argv[2]
fai=sys.argv[3]
LdSnpDIR=sys.argv[4]
cellTypePeakDIR=sys.argv[5]
maxIter=int(sys.argv[6])



def loadSigPeaks():
    PeakDict={}
    input1=open(cellTypePeakDIR+'/'+cellType+'.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	chrID=each[0]
	pt1=each[1]
	pt2=each[2]
	peakSize=int(pt2)-int(pt1)

	if PeakDict.has_key(chrID):
	    tmp=PeakDict[chrID]
	    tmp.append(peakSize)
	    PeakDict[chrID]=tmp
	else:
	    PeakDict[chrID]=[peakSize]
    input1.close()

    return PeakDict


## Number of peaks matched in each chr
def genRandPeak():
    PeakDict = loadSigPeaks()

    genomeDict={}
    input1=open(fai,'r')
    all_input1=input1.readlines()
    for line in all_input1:
        each=line.strip().split('\t')
        chrID=each[0]
        chrSize=each[1]
	genomeDict[chrID]=chrSize
    input1.close()

    output1=open(DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.RandPeak.bed','w')
    for chrID in PeakDict:
        for peakSize in PeakDict[chrID]:
            pt1=random.randrange(0, int(genomeDict[chrID])-peakSize)
            pt2=pt1+peakSize
            output1.write(chrID+'\t'+str(pt1)+'\t'+str(pt2)+'\n')
    output1.close()



def permuteRandPeaks():
    output1=open(DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.RandPeakPermuteRecord.txt','w')
    output1.write('nPermute\tnPeak\n')
    for nIter in range(0, maxIter):
	print str(nIter+1)+'/'+str(maxIter)

	genRandPeak()

	##
	runScript='intersectBed -a '+DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.RandPeak.bed -b '+LdSnpDIR+'/'+gwasID+'.bed -wa -wb | cut -f1,2,3 | sort -u | wc -l > '+DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.RandPeak.intersect.bed'
	#print runScript
	os.system(runScript)

	##
	input1=open(DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.RandPeak.intersect.bed', 'r')
	all_input1=input1.readlines()
	nPeak=all_input1[0].strip().split('\t')[0]
	input1.close()

	newLine=[str(nIter+1), nPeak]
	output1.write('\t'.join(newLine)+'\n')

    output1.close()




permuteRandPeaks()


