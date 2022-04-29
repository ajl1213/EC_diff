#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import operator


DIR=os.getcwd()
inputDIR=DIR+'/RESULT'
outputDIR=DIR+'/ValMat'
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

cellTypeList=['EarlyEC', 'MidEC', 'LateEC', 'FullEC']

minAbcScore=5
maxGeneCount=10

os.system('mkdir '+outputDIR)



def makeValMat():

    for gwasID in gwasList:
	# load value
	valDict={}
	geneList=[]

	for cellType in cellTypeList:
	    input1=open(inputDIR+'/'+gwasID+'.PeakToGene.final.txt', 'r')
	    all_input1=input1.readlines()
	    tmpDict={}
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		peakID=each[0]
		peakLabel=each[1]
		ensembleID=each[2]
		geneID=each[3]
		abcScore=float(each[5])

		if peakLabel==cellType:
		    if tmpDict.has_key(geneID):
			tmpDict[geneID]+=abcScore
		    else:
			tmpDict[geneID]=abcScore

		    if not geneID in geneList:
			geneList.append(geneID)

		valDict[cellType]=tmpDict
	    input1.close()

        # filter
        filteredList=[]
        for geneID in geneList:
            tmpList=[]
            for cellType in cellTypeList:
		if valDict[cellType].has_key(geneID):
		    val=valDict[cellType][geneID]
		else:
		    val=0
		    valDict[cellType][geneID]=0

                tmpList.append(val)
            if max(tmpList) > minAbcScore:
                filteredList.append(geneID)

        # order
        assignDict={}
        for geneID in filteredList:
            finalCellType='none'
            cur_val=0
            for cellType in cellTypeList:
                val=float(valDict[cellType][geneID])
                if val > cur_val:
                    finalCellType=cellType
                    cur_val=val
            assignDict[geneID]=finalCellType

        output1=open(outputDIR+'/'+gwasID+'.txt','w')
        output1.write('geneID\tmaxCell\t'+'\t'.join(cellTypeList)+'\n')
        for cellType in cellTypeList:
            tmpDict={}
            for geneID in filteredList:
                if cellType == assignDict[geneID]:
                    tmpDict[geneID]=float(valDict[cellType][geneID])

            for key in sorted(tmpDict.items(), key=operator.itemgetter(1), reverse=True)[:maxGeneCount]:
                geneID=key[0]

                valList=[]
                for cellType in cellTypeList:
                    val=str(valDict[cellType][geneID])
                    valList.append(val)
                output1.write(geneID+'\t'+assignDict[geneID]+'\t'+'\t'.join(valList)+'\n')
        output1.close()



makeValMat()


