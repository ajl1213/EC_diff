#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import operator


DIR=os.getcwd()
GwasTargetInfo='/home/ajl1213/Projects/Endothelial/Analysis2/GwasTarget/02_GetTargetGenes/TargetGeneList/Colorectal_cancer.MergedTarget.abcScore.txt'
ScDegFile='/home/ajl1213/Projects/Endothelial/data/snRNA/ColonCancer_GSE132465/2_Analysis.harmony/CellTypeDEGs/snRNA.cellTypeDEGs.txt'
PseudoTimeInfo=DIR+'/pseudotime.metadata.txt'
ZscoreMat=DIR+'/ValMat/zmat.txt'


def loadLog2fc():
    log2fcDict={}
    input1=open(ScDegFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[p_val, avg_log2FC, pct1, pct2, p_val_adj, diff, celltype, geneID]=line.strip().split('\t')
	if celltype=='EndothelialCells':
	    log2fcDict[geneID]=float(avg_log2FC)
    input1.close()

    return log2fcDict


def loadZscore():
    zscoreDict={}
    input1=open(ZscoreMat,'r')
    all_input1=input1.readlines()
    cellList=all_input1[0].strip().split('\t')[1:] 
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	vals=each[1:]
	tmpDict={}
	idx=0
	for val in vals:
	    tmpDict[cellList[idx]]=val
	    idx+=1

	zscoreDict[geneID]=tmpDict

    input1.close()

    return zscoreDict


def loadPseudoTime():
    barcodeList=[]
    input1=open(PseudoTimeInfo,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	barcodeID=each[0]
	barcodeList.append(barcodeID)
    input1.close()

    return barcodeList


def loadTargetGeneInfo():
    geneSet=[]

    #
    degDict={}
    input1=open(DIR+'/GenesOverlap/DEG_intersect.txt', 'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	label=each[1]
	degDict[geneID]=label

	geneSet.append(geneID)
    input1.close()

    #
    trajDict={}
    input1=open(DIR+'/GenesOverlap/Traj_intersect.txt', 'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	label=each[1]
	trajDict[geneID]=label

	geneSet.append(geneID)
    input1.close()

    #
    gwasDict={}
    input1=open(GwasTargetInfo,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[peakID, peakLabel, ensembleID, geneID, assoType, abcScore]=line.strip().split('\t')
	if gwasDict.has_key(geneID):
	    tmp=gwasDict[geneID]
	    if not peakLabel in tmp:
		tmp.append(peakLabel)
	    gwasDict[geneID]=tmp
	else:
	    gwasDict[geneID]=[peakLabel]
    input1.close()

    output1=open(DIR+'/GenesOverlap/MergedList.final.txt','w')
    output1.write('geneID\tpeakLabel\tDEG_label\tTraj_label\n')
    for geneID in set(geneSet):

	if len(gwasDict[geneID]) > 1:
#	    peakLabel='multiple'	    
	    peakLabel=';'.join(gwasDict[geneID])
	else:
	    peakLabel=gwasDict[geneID][0]

	if degDict.has_key(geneID):
	    deg_label=degDict[geneID]
	else:
	    deg_label='none'
    
	if trajDict.has_key(geneID):
	    traj_label=trajDict[geneID]
	else:
	    traj_label='none'

	newLine=[geneID, peakLabel, deg_label, traj_label]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()

    infoDict={}
    input1=open(DIR+'/GenesOverlap/MergedList.final.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[geneID, peakLabel, deg_label, traj_label]=line.strip().split('\t')
	infoDict[geneID]=[peakLabel, deg_label, traj_label]

    input1.close()

    return infoDict


def makeData():
    log2fcDict = loadLog2fc()
    barcodeList = loadPseudoTime()
    zscoreDict = loadZscore()
    infoDict = loadTargetGeneInfo()

    # 
    output1=open(DIR+'/ValMat/targetGene.log2fc.order.txt','w')
    output1.write('geneID\tpeakLabel\tDEG\ttrajGene\tlog2fc\t'+'\t'.join(barcodeList)+'\n')
    for key in sorted(log2fcDict.items(), key=operator.itemgetter(1), reverse=True):
	geneID=key[0]
	log2fc=key[1]
	if infoDict.has_key(geneID):
	    valList=[]
	    for cellID in barcodeList:
		val=zscoreDict[geneID][cellID]
		valList.append(val)		    
	    [peakLabel, deg_label, traj_label]=infoDict[geneID]

	    output1.write(geneID+'\t'+peakLabel+'\t'+deg_label+'\t'+traj_label+'\t'+str(log2fc)+'\t'+'\t'.join(valList)+'\n')
    output1.close()


def smoothing():
    nInterval=10

    input1=open(DIR+'/ValMat/targetGene.log2fc.order.txt','r')
    all_input1=input1.readlines()
    nVal=int(float(len(all_input1[0].strip().split('\t')[5:]))/nInterval)
    output1=open(DIR+'/ValMat/targetGene.log2fc.order.smooth.txt','w')
    output1.write('geneID\tpeakLabel\tDEG\ttrajGene\tlog2fc')
    for i in range(nVal+1):
	i=i+1
	output1.write('\tcell_'+str(i))
    output1.write('\n')

    valDict={}
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	peakLabel=each[1]
	deg_label=each[2]
	traj_label=each[3]
	log2fc=each[4]
	vals=each[5:]
	tmpDict={}
	idx=1
	for val in vals:
	    if tmpDict.has_key(idx):
		if len(tmpDict[idx]) < nInterval:
		    tmp=tmpDict[idx]
		    tmp.append(val)
		    tmpDict[idx]=tmp
		else:
		    idx+=1
		    tmpDict[idx]=[val]

	    else:
		tmpDict[idx]=[val]

	valList=[]
	for key in sorted(tmpDict.items(), key=operator.itemgetter(0)):
	    meanVal=numpy.mean(numpy.array(tmpDict[key[0]], dtype='float'))
	    valList.append(str(meanVal))
	output1.write(geneID+'\t'+peakLabel+'\t'+deg_label+'\t'+traj_label+'\t'+str(log2fc)+'\t'+'\t'.join(valList)+'\n')

    input1.close()
    output1.close()


makeData()
smoothing()


