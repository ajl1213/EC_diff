#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
rawSnpDIR='/home/ajl1213/Projects/Endothelial/Analysis2/GwasTarget/01_MergeSNP/MergedRawSNP'
peakIntersectDIR='/home/ajl1213/Projects/Endothelial/Analysis2/GwasTarget/01_MergeSNP/ObservedData'
peakLabelFile='/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/DiffPeaks/H3K27ac.DiffPeakLabel.txt'

gwasList=[
'Crohns_disease',
'Asthma',
'Inflammatory_bowel_disease',
'Multiple_sclerosis',
'Eczema',
'Atrial_fibrillation',
'Atopic_asthma',
'Coronary_artery_disease',
'Varicose_veins',
'Alzheimers_disease_late_onset',
'Colorectal_cancer',
'Ulcerative_colitis',
'Rheumatoid_arthritis'
]

os.system('mkdir '+DIR+'/RawSNP')
os.system('mkdir '+DIR+'/RESULT')



def makeData():

    for gwasID in gwasList:
	print gwasID

	# load SNP label
	snpLabelDict={}
	tagSnpList=[]
	input1=open(rawSnpDIR+'/'+gwasID+'.rawSNP.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    chrID=each[0]
	    pt=each[1]
	    snpID=chrID+'.'+pt
	    jList=each[2]
	    snpLabelDict[snpID]=jList
	    tagSnpList.append(snpID)
	input1.close()

	# load peak label
	peakLabelDict={}
	input1=open(peakLabelFile,'r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    peakID=each[0]
	    peakLabel=each[2]
	    peakLabelDict[peakID]=peakLabel
	input1.close()

	# load intersect info
	intersectDict={}
	input1=open(peakIntersectDIR+'/'+gwasID+'.MergedEC.Peak.intersect.full.bed','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    peakID=each[3]
	    ldSNP=each[4]+'.'+each[5]
	    tagSNP=each[7]
	    if intersectDict.has_key(tagSNP):
		tmp=intersectDict[tagSNP]
		if not peakID in tmp:
		    tmp.append(peakID)
		intersectDict[tagSNP]=tmp
	    else:
		intersectDict[tagSNP]=[peakID]
	input1.close()

	# 
	peakCountDict={}
	for tagSNP in intersectDict:
	    for peakID in intersectDict[tagSNP]:
		if peakCountDict.has_key(peakID):
		    peakCountDict[peakID]+=1
		else:
		    peakCountDict[peakID]=1

	# make data
	output1=open(DIR+'/RESULT/'+gwasID+'.SnpToPeak.txt','w')
	output1.write('tagSNP\tsnpLabel\tpeakID\tpeakLabel\tfrac\n')
	for tagSNP in tagSnpList:
	    snpLabel=snpLabelDict[tagSNP]

	    if intersectDict.has_key(tagSNP):
		for peakID in intersectDict[tagSNP]:
		    peakLabel=peakLabelDict[peakID]
		    frac=float(1)/float(peakCountDict[peakID])

		    newLine=[tagSNP, snpLabel, peakID, peakLabel, str(frac)]
		    output1.write('\t'.join(newLine)+'\n')

	output1.close()



makeData()



