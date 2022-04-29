#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
inputDIR=DIR+'/TargetGeneList'
AbcScoreDIR='/home/ajl1213/Projects/Endothelial/data/HiC/ActivityByContact/ABC_score'
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

cellTypeList=['Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC']
minAbcScore=0



def loadAbcScore():
    abcScoreDict={}

    print 'Loading ABC score..'
    for cellType in cellTypeList:
        print cellType

        tmpDict={}
        input1=open(AbcScoreDIR+'/'+cellType+'.abcScore.txt','r')
        all_input1=input1.readlines()
        for line in all_input1[1:]:
            each=line.strip().split('\t')

            ensembleID=each[0]
            peakID=each[3]
            contactVal=float(each[8])
            abcScore=float(each[9])

            keyID=ensembleID+';'+peakID
            tmpDict[keyID]=abcScore
        input1.close()

        abcScoreDict[cellType]=tmpDict
    print 'Loading complete'

    return abcScoreDict


def attachAbcScore():
    abcScoreDict=loadAbcScore()

    for gwasID in gwasList:
	input1=open(inputDIR+'/'+gwasID+'.MergedTarget.txt','r')
	all_input1=input1.readlines()
	output1=open(inputDIR+'/'+gwasID+'.MergedTarget.abcScore.txt','w')
	output1.write('peakID\tpeakLabel\tensembleID\tgeneID\tassoType\tabcScore\n')
	for line in all_input1[1:]:
	    [peakID, peakLabel, ensembleID, geneID, assoType]=line.strip().split('\t')
	    if abcScoreDict[peakLabel].has_key(ensembleID+';'+peakID):
		abcScore=abcScoreDict[peakLabel][ensembleID+';'+peakID]

		newLine=[peakID, peakLabel, ensembleID, geneID, assoType, str(abcScore)]
		output1.write('\t'.join(newLine)+'\n')
	    else:
		print peakLabel, ensembleID, geneID, peakID

	input1.close()


def makeFinalSet():

    for gwasID in gwasList:
	peakList=[]

	#
	tmpDict={}	
	input1=open(inputDIR+'/'+gwasID+'.MergedTarget.abcScore.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    [peakID, peakLabel, ensembleID, geneID, assoType, abcScore]=line.strip().split('\t')
	    if float(abcScore) > minAbcScore:
		peakList.append(peakID)
		if tmpDict.has_key(peakID):
		    tmpDict[peakID]+=1
		else:
		    tmpDict[peakID]=1
	input1.close()

	output1=open(DIR+'/RESULT/'+gwasID+'.PeakToGene.final.txt','w')
	output1.write('peakID\tpeakLabel\tensembleID\tgeneID\tassoType\tabcScore\tfrac\n')
	for line in all_input1[1:]:
	    [peakID, peakLabel, ensembleID, geneID, assoType, abcScore]=line.strip().split('\t')
	    if float(abcScore) > minAbcScore:
		frac=float(1)/float(tmpDict[peakID])
		newLine=[peakID, peakLabel, ensembleID, geneID, assoType, abcScore, str(frac)]
		output1.write('\t'.join(newLine)+'\n')
	output1.close()


	#
	peakCountDict={}
	input1=open(DIR+'/RESULT/'+gwasID+'.SnpToPeak.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    tagSNP=each[0]
	    snpLabel=each[1]
	    peakID=each[2]
	    peakLabel=each[3]
	    if peakID in peakList:
		if peakCountDict.has_key(peakID):
		    peakCountDict[peakID]+=1
		else:
		    peakCountDict[peakID]=1
	input1.close()


	# 
	input1=open(DIR+'/RESULT/'+gwasID+'.SnpToPeak.txt','r')
	all_input1=input1.readlines()
	output1=open(DIR+'/RESULT/'+gwasID+'.SnpToPeak.final.txt','w')
	output1.write('tagSNP\tsnpLabel\tpeakID\tpeakLabel\tfrac\n')
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    tagSNP=each[0]
	    snpLabel=each[1]
	    peakID=each[2]
	    peakLabel=each[3]
	    if peakID in peakList:
		frac=float(1)/float(peakCountDict[peakID])
		newLine=[tagSNP, snpLabel, peakID, peakLabel, str(frac)]
		output1.write('\t'.join(newLine)+'\n')
	input1.close()
	output1.close()



attachAbcScore()
makeFinalSet()


