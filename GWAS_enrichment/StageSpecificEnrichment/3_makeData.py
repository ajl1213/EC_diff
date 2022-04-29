#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import random



DIR=os.getcwd()
LdSnpDIR='/home/ajl1213/Projects/Endothelial/Analysis2/GwasEnrich/1_getGwasCatalog/LdFinal'
PeakDIR='/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/PeakBED'

cellTypeList=['hESC','Mesoderm','EarlyEC','MidEC','LateEC','FullEC']
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

os.system('mkdir '+DIR+'/ObservedData')
os.system('mkdir '+DIR+'/RESULT')



def getActualData():
    for cellType in cellTypeList:
	for gwasID in gwasList:
	    ##  
	    runScript='intersectBed -a '+PeakDIR+'/'+cellType+'.bed -b '+LdSnpDIR+'/'+gwasID+'.bed -wa -wb | cut -f1,2,3 | sort -u | wc -l > '+DIR+'/ObservedData/'+gwasID+'.'+cellType+'.nPeakIntersect.txt'
	    print runScript
	    os.system(runScript)

	    runScript='intersectBed -a '+PeakDIR+'/'+cellType+'.bed -b '+LdSnpDIR+'/'+gwasID+'.bed -wa -wb | cut -f1,2,3 | sort -u > '+DIR+'/ObservedData/'+gwasID+'.'+cellType+'.Peak.intersect.bed'
	    print runScript
	    os.system(runScript)


def makeData():

    for cellType in cellTypeList:
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

	    if not len(tmpList) == 0 :
		pval=float(len(tmpList))/float(len(randValList))
	    else:
		pval=1e-6

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


    #   
    output1=open(DIR+'/RESULT/MergedDataList.txt','w')
    output1.write('gwasID\tcellType\tnPeak\tfcVal\tempiricalPval\n')
    for cellType in cellTypeList:
        input1=open(DIR+'/RESULT/'+cellType+'.DataList.txt','r')
        all_input1=input1.readlines()
        for line in all_input1[1:]:
            each=line.strip().split('\t')
            gwasID=each[0]
            nPeak=each[1]
            fcVal=each[2]
            empiricalPval=each[3]
            newLine=[gwasID, cellType, nPeak, fcVal, empiricalPval]
            output1.write('\t'.join(newLine)+'\n')
        input1.close()
    output1.close()



getActualData()
makeData()



