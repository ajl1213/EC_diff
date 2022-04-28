#!/python2


import os
import numpy
import operator



DIR=os.getcwd()
dataTypeList=['Mesoderm','EarlyEC','MidEC','LateEC','FullEC']
nLimit=20


def combineGO():
    for dataType in dataTypeList:
	tmpDict={}
	input1=open(DIR+'/GREAT_GOBP/'+dataType+'.GOBP.txt','r')
	output1=open(DIR+'/GREAT_GOBP/'+dataType+'.GOBP.n'+str(nLimit)+'.txt','w')
	output1.write('BPID\tDescription\tBinomPval(-log10)\tRegionFoldEnrich\tGeneFoldEnrich\tObsGenes\tTotalGenes\n')
	all_input1=input1.readlines()
	idx=0
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    traitID=each[1]
	    termID=each[2]
	    termID=termID.replace(' ','_')
	    logPval=-(numpy.log2(float(each[3])))
	    regionFC=each[4]
	    geneFC=each[5]
	    nObsGenes=each[6]
	    nTotalGenes=each[7]
	    if idx < nLimit:
		newLine=[traitID, termID, str(logPval), regionFC, geneFC, nObsGenes, nTotalGenes]	
		output1.write('\t'.join(newLine)+'\n')
	    idx+=1
	input1.close()
	output1.close()


def makeData():
    selectedList=[
    'skeletal_system_development',
    'bone_morphogenesis',
    'Wnt_signaling_pathway',
    'limb_development',
    'limb_morphogenesis',
    'cardiac_right_ventricle_morphogenesis',
    'artery_development',
    'cardiac_septum_development',
    'apoptotic_process_involved_in_development',
    'heart_morphogenesis',
    'tube_development',
    'epithelium_development',
    'positive_regulation_of_signal_transduction',
    'tissue_morphogenesis',
    'morphogenesis_of_an_epithelium',
    'vasculature_development',
    'blood_vessel_development',
    'blood_vessel_morphogenesis',
    'angiogenesis',
    'circulatory_system_development',
    'vasculature_development',
    'blood_vessel_development',
    'regulation_of_angiogenesis',
    'positive_regulation_of_vasculature_development',
    'blood_vessel_morphogenesis'
    ]

    termList=[]
    for i in selectedList:
	if not i in termList:
	    termList.append(i)

    allList=[]
    for dataType in dataTypeList:
	tmpDict={}
	input1=open(DIR+'/GREAT_GOBP/'+dataType+'.GOBP.n'+str(nLimit)+'.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    
	    termID=each[1]
	    logPval=each[2]

	    tmpDict[termID]=logPval

	allList.append(tmpDict)

	input1.close()

    output1=open(DIR+'/GREAT_GOBP/EC_cRE.GOTABLE.txt','w')
    output1.write('GOBP\t'+'\t'.join(dataTypeList)+'\n')
    for termID in termList:
	output1.write(termID)
	for i in allList:
	    if i.has_key(termID):
		output1.write('\t'+str(i[termID]))
	    else:
		output1.write('\t0.000001')
	output1.write('\n')
    output1.close()	



combineGO()
makeData()




