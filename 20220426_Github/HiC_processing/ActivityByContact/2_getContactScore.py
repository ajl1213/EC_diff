#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import gzip



DIR=os.getcwd()
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
InterCountDIR='/home/ajl1213/Projects/Endothelial/data/HiC/CovNorm/FeatureVec'
PeakCountFile='/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/CountMatrix/H3K27ac.RPM.qq.txt'
chrList=['chr1', 'chr2', 'chr3','chr4', 'chr5', 'chr6','chr7', 'chr8', 'chr9','chr10', 'chr11', 'chr12','chr13', 'chr14', 'chr15','chr16', 'chr17', 'chr18','chr19', 'chr20', 'chr21','chr22', 'chrX']
inputFile=DIR+'/CreListPerGene/CreListPerGene.txt'
resolution=['10kb',10000]
cellTypeList=['Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC']

os.system('mkdir '+DIR+'/ABC_score')


def loadGeneInfo():
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


def loadInterCount():
    interCountDict={}
    print 'loading HiC count..'
    for chrID in chrList:
	print chrID
	input1=gzip.open(InterCountDIR+'/MergedEC.'+chrID+'.'+resolution[0]+'.count.gz','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    frag1=each[0]
	    frag2=each[1]
	    freqVal=each[4]
	    interCountDict[frag1+'-'+frag2]=float(freqVal)
	    interCountDict[frag2+'-'+frag1]=float(freqVal)
	input1.close()
    print 'done'

    return interCountDict


def loadPeakCount():
    print 'loading Peak count..'
    valDict={}
    input1=open(PeakCountFile,'r')
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

    peakCountDict={}
    for peakID in valDict:
	tmpDict={}
	for cellType in cellTypeList:
	    valList=[]
	    for sampleID in sampleList:
		if sampleID.split('_')[0] in ['ME']:
		    dataType='Mesoderm'
		elif sampleID.split('_')[0] in ['EC01','EC02','EC04']:
		    dataType='EarlyEC'
		elif sampleID.split('_')[0] in ['EC06','EC08','EC12']:
		    dataType='MidEC'
		elif sampleID.split('_')[0] in ['EC24']:
		    dataType='LateEC'
		elif sampleID.split('_')[0] in ['EC48']:
		    dataType='FullEC'
		else:
		    dataType=''

		if dataType==cellType:
		    val=valDict[peakID][sampleID]
		    valList.append(val)
	    meanVal=numpy.mean(valList)
	    tmpDict[cellType]=meanVal

        peakCountDict[peakID]=tmpDict

    print 'done'

    return peakCountDict


def computeContactScore():
    geneDict=loadGeneInfo()
    interCountDict=loadInterCount()
    peakCountDict=loadPeakCount()

    print 'computing ABC score..'
    for cellType in cellTypeList:
	print cellType

	totalScoreDict={}
	input1=open(inputFile,'r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    tssBin=each[2]
	    creID=each[3]
	    creBin=each[4]

	    peakCount=peakCountDict[creID][cellType]

	    interID=tssBin+'-'+creBin
	    dist=int(each[5])

	    if interCountDict.has_key(interID):
		if dist > 0:
		    hicCount=float(interCountDict[interID])/2
		else:
		    frag1=tssBin.split('.')[0]+'.'+str(int(tssBin.split('.')[1])-resolution[1])+'.'+str(int(tssBin.split('.')[2])-resolution[1])
		    frag2=tssBin.split('.')[0]+'.'+str(int(tssBin.split('.')[1])+resolution[1])+'.'+str(int(tssBin.split('.')[2])+resolution[1])
		    interID1=frag1+'-'+tssBin
		    interID2=frag2+'-'+tssBin

		    count1=float(interCountDict[interID1])/2
		    count2=float(interCountDict[interID2])/2

		    hicCount=float(count1+count2)/2  ###

		contactVal=hicCount*peakCount
		if totalScoreDict.has_key(ensembleID):
		    totalScoreDict[ensembleID]+=contactVal
		else:
		    totalScoreDict[ensembleID]=contactVal

	    else:
		print 'no Hi-C count : '+ensembleID+' '+tssBin+' '+creID+' '+creBin+' '+str(peakCount)+' '+str(dist)

	input1.close()

	##
	output1=open(DIR+'/ABC_score/'+cellType+'.SumContactPerGene.txt','w')
	output1.write('ensembleID\tgeneID\tsumContactVal\n')
	for ensembleID in totalScoreDict:
	    geneID=geneDict[ensembleID]
	    sumContactVal=str(totalScoreDict[ensembleID])
	    newLine=[ensembleID, geneID, sumContactVal]
	    output1.write('\t'.join(newLine)+'\n')
	output1.close()

	##
	output1=open(DIR+'/ABC_score/'+cellType+'.abcScore.txt','w')
	output1.write('ensembleID\tgeneID\ttssBin\tcreID\tcreBin\tdist\tpeakCount\tinterCount\tcontactVal\tabcScore\n')
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    geneID=each[1]
	    tssBin=each[2]
	    creID=each[3]
	    creBin=each[4]

	    peakCount=peakCountDict[creID][cellType]

	    interID=tssBin+'-'+creBin
	    dist=int(each[5])

	    if interCountDict.has_key(interID):
		if dist > 0:
		    hicCount=float(interCountDict[interID])/2
		else:
		    frag1=tssBin.split('.')[0]+'.'+str(int(tssBin.split('.')[1])-resolution[1])+'.'+str(int(tssBin.split('.')[2])-resolution[1])
		    frag2=tssBin.split('.')[0]+'.'+str(int(tssBin.split('.')[1])+resolution[1])+'.'+str(int(tssBin.split('.')[2])+resolution[1])
		    interID1=frag1+'-'+tssBin
		    interID2=frag2+'-'+tssBin

		    count1=float(interCountDict[interID1])/2
		    count2=float(interCountDict[interID2])/2

		    hicCount=float(count1+count2)/2  ###

		contactVal=hicCount*peakCount
		sumScore=totalScoreDict[ensembleID]
		activityByContactScore= contactVal/sumScore*100 if not sumScore==0 else 0

	    newLine=[ensembleID, geneID, tssBin, creID, creBin, str(dist), str(peakCount), str(hicCount), str(contactVal), str(activityByContactScore)]
	    output1.write('\t'.join(newLine)+'\n')

	output1.close()



computeContactScore()



