#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
DataInfo=DIR+'/DataInfo.txt'
inputMatrix=DIR+'/CountMatrix/H3K27ac.RPM.qq.txt'
stageList=['hESC', 'ME', 'EC01', 'EC02', 'EC04', 'EC06', 'EC08', 'EC12', 'EC24', 'EC48', 'nonEC01', 'nonEC02', 'nonEC04', 'nonEC06', 'nonEC08', 'nonEC12', 'nonEC24', 'nonEC48']
dataTypeList=['hESC', 'Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC', 'EarlyNonEC', 'MidNonEC', 'LateNonEC', 'FullNonEC']

os.system('mkdir '+DIR+'/PeakBED')



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
    sampleList=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	peakID=each[0]
	vals=numpy.array(each[1:], dtype='float')
	tmpDict={}
	index=0
	for val in vals:
	    tmpDict[sampleList[index]]=val
	    index+=1
	valDict[peakID]=tmpDict
    input1.close()

    return sampleList, valDict


def loadPeak():
    peakList=[]
    for dataType in dataTypeList:
	input1=open(DIR+'/DiffPeaks/'+dataType+'.valTable.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    peakID=each[0]
	    peakList.append(peakID)
	input1.close()

    peakList=set(peakList)

    return peakList


def makeLabel():
    dataDict = loadDataInfo()
    sampleList, valDict= loadCountData()
    peakList = loadPeak()

    output1=open(DIR+'/DiffPeaks/H3K27ac.DiffPeakLabel.txt', 'w')
    output1.write('peakID\tmaxStage\tmaxDataType\n')

    for peakID in peakList:
	tmpDict={}
	for stageLabel in stageList:
	    tmpList=[]
	    for sampleID in sampleList:
		if dataDict[sampleID][0]==stageLabel:
		    val=valDict[peakID][sampleID]
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

	newLine=[peakID, maxStage, maxDataType]
	output1.write('\t'.join(newLine)+'\n')

    output1.close()


def makeAllPeakMatrix():
    sampleList, valDict= loadCountData()

    meanDict={}
    stdDict={}
    for peakID in valDict:
	tmpList=[]
	for sampleID in sampleList:
	    val=valDict[peakID][sampleID]
	    tmpList.append(val)
	meanVal=numpy.mean(tmpList)
	stdVal=numpy.std(tmpList)

	meanDict[peakID]=meanVal
	stdDict[peakID]=stdVal

    output1=open(DIR+'/CountMatrix/H3K27ac.readCount.qq.zscore.txt', 'w')
    output1.write('peakID\t'+'\t'.join(sampleList)+'\n')
    for peakID in valDict:
	output1.write(peakID)
	for sampleID in sampleList:
	    val=valDict[peakID][sampleID]
	    zscore=str((val-meanDict[peakID])/stdDict[peakID])
	    output1.write('\t'+zscore)
	output1.write('\n')
    output1.close()


def makeDiffMatrix():
    dataDict = loadDataInfo()
    sampleList, valDict= loadCountData()

    ## Load DiffPeak
    peakDict={}
    input1=open(DIR+'/DiffPeaks/H3K27ac.DiffPeakLabel.txt', 'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	peakID=each[0]
	maxStage=each[1]
	if peakDict.has_key(maxStage):
	    tmp=peakDict[maxStage]
	    tmp.append(peakID)
	    peakDict[maxStage]=tmp
	else:
	    peakDict[maxStage]=[peakID]
    input1.close()

    ## Make Zscore
    output1=open(DIR+'/DiffPeaks/H3K27ac.DiffPeakLabel.zscore.txt', 'w')
    output1.write('peakID\tmaxStage\t'+'\t'.join(stageList)+'\n')
    for i in stageList:
	for peakID in peakDict[i]:
	    tmpDict={}
	    for stageLabel in stageList:
		tmpList=[]
		for sampleID in sampleList:
		    if dataDict[sampleID][0]==stageLabel:
			val=valDict[peakID][sampleID]
			tmpList.append(val)
		meanVal=numpy.mean(tmpList)
		tmpDict[stageLabel]=meanVal


	    newLine=[peakID, i]

	    valList=[]
	    for stageLabel in stageList:
		val=tmpDict[stageLabel]
		valList.append(val)

	    meanVal=numpy.mean(valList)
	    stdVal=numpy.std(valList)

	    for stageLabel in stageList:
		val=tmpDict[stageLabel]
		zscore=str((val-meanVal)/stdVal)
		newLine.append(zscore)
		
	    output1.write('\t'.join(newLine)+'\n')
    output1.close()


def separateBED():
    ##
    for dataType in dataTypeList:
        input1=open(DIR+'/DiffPeaks/H3K27ac.DiffPeakLabel.txt','r')
        output1=open(DIR+'/PeakBED/'+dataType+'.bed','w')
        all_input1=input1.readlines()
        for line in all_input1[1:]:
            each=line.strip().split('\t')
            peakID=each[0]
            if dataType==each[2]:
                chrID=peakID.split(':')[0]
                pt1=peakID.split(':')[1].split('-')[0]
                pt2=peakID.split(':')[1].split('-')[1]

                newLine=[chrID, pt1, pt2, peakID]
                output1.write('\t'.join(newLine)+'\n')
        input1.close()
        output1.close()

    ##
    command1=['cat ']
    for dataType in ['EarlyEC', 'MidEC', 'LateEC', 'FullEC']:
        command1.append(DIR+'/PeakBED/'+dataType+'.bed ')

    command1.append('| cut -f1,2,3 | sortBed -i stdin | mergeBed -i stdin | sort -V -k1,1 -k2,2n - > '+DIR+'/PeakBED/AllEC.bed')
    command1=''.join(command1)
    print command1
    os.system(command1)

    ##
    for stage in stageList:
        input1=open(DIR+'/DiffPeaks/H3K27ac.DiffPeakLabel.txt','r')
        output1=open(DIR+'/PeakBED/'+stage+'.bed','w')
        all_input1=input1.readlines()
        for line in all_input1[1:]:
            each=line.strip().split('\t')
            peakID=each[0]
            if stage == each[1]:
                chrID=peakID.split(':')[0]
                pt1=peakID.split(':')[1].split('-')[0]
                pt2=peakID.split(':')[1].split('-')[1]

                newLine=[chrID, pt1, pt2, peakID]
                output1.write('\t'.join(newLine)+'\n')
        input1.close()
        output1.close()




#makeLabel()
#makeAllPeakMatrix()
#makeDiffMatrix()
separateBED()



