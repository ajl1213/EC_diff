#!/home/ajl1213/anaconda2/bin/python



import os
import numpy



DIR=os.getcwd()
inputDIR='/home/ajl1213/Projects/Endothelial/data/ChIP/Mapping/PeakCoverage'
histMark='H3K27ac'

##
sampleList=[
'hESC_rep1',
'hESC_rep2',
'ME_rep1',
'ME_rep2',
'EC01_rep1',
'EC01_rep2',
'EC01_rep3',
'EC02_rep1',
'EC02_rep2',
'EC02_rep3',
'EC04_rep1',
'EC04_rep2',
'EC06_rep1',
'EC06_rep2',
'EC08_rep1',
'EC08_rep2',
'EC12_rep1',
'EC12_rep2',
'EC24_rep1',
'EC24_rep2',
'EC48_rep1',
'EC48_rep2',
'EC48_rep3',
'nonEC01_rep1',
'nonEC01_rep2',
'nonEC02_rep1',
'nonEC02_rep2',
'nonEC02_rep3',
'nonEC04_rep1',
'nonEC04_rep2',
'nonEC06_rep1',
'nonEC06_rep2',
'nonEC08_rep1',
'nonEC08_rep2',
'nonEC12_rep1',
'nonEC12_rep2',
'nonEC24_rep1',
'nonEC24_rep2',
'nonEC48_rep1',
'nonEC48_rep2'
]

print len(sampleList)

os.system('mkdir ' +DIR+'/CountMatrix')



def makeData():
    output1=open(DIR+'/DataInfo.txt','w')
    output1.write('sampleID\tstageLabel\tdataType\tcellType\n')
    for sampleID in sampleList:
	stageLabel=sampleID.split('_')[0]

	if stageLabel in ['hESC']:
	    dataType='hESC'
	    cellType='hESC'
	if stageLabel in ['ME']:
	    dataType='Mesoderm'
	    cellType='Mesoderm'
	if stageLabel in ['EC01','EC02','EC04']:
	    dataType='EarlyEC'
	    cellType='EC'
	if stageLabel in ['EC06','EC08','EC12']:
	    dataType='MidEC'
	    cellType='EC'
	if stageLabel in ['EC24']:
	    dataType='LateEC'
	    cellType='EC'
	if stageLabel in ['EC48']:
	    dataType='FullEC'
	    cellType='EC'
	if stageLabel in ['nonEC01','nonEC02','nonEC04']:
	    dataType='EarlyNonEC'
	    cellType='nonEC'
	if stageLabel in ['nonEC06','nonEC08','nonEC12']:
	    dataType='MidNonEC'
	    cellType='nonEC'
	if stageLabel in ['nonEC24']:
	    dataType='LateNonEC'
	    cellType='nonEC'
	if stageLabel in ['nonEC48']:
	    dataType='FullNonEC'
	    cellType='nonEC'

	newLine=[sampleID, stageLabel, dataType, cellType]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()



def makeCountMatrix():
    countDict={}
    for sampleID in sampleList:
        input1=open(inputDIR+'/'+sampleID+'.'+histMark+'.peakCoverage.txt','r')
        all_input1=input1.readlines()
	tmpDict={}
        for line in all_input1:
            [chrID, pos1, pos2, readCount]=line.strip().split('\t')
	    coordID=chrID+':'+pos1+'-'+pos2	
	    tmpDict[coordID]=str(float(readCount))
	countDict[sampleID]=tmpDict
        input1.close()

    output1=open(DIR+'/CountMatrix/'+histMark+'.readCount.txt','w')
    output1.write('Coordinate'+'\t'+'\t'.join(sampleList)+'\n')
    for coordID in countDict[sampleList[0]]:
	output1.write(coordID)
	for sampleID in sampleList:
	    val=countDict[sampleID][coordID]
	    output1.write('\t'+str(val))
	output1.write('\n')




makeData()
makeCountMatrix()



