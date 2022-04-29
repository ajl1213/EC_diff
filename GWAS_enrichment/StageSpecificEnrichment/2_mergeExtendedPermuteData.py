#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
fai='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
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

maxCopy=10

for cellType in cellTypeList:
    for gwasID in gwasList:

	output1=open(DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.RandPeakPermuteRecord.txt','w')
	output1.write('nPermute\tnPeak\n')

	for j in range(maxCopy):
	    nCopy=str(j+1)

	    input1=open(DIR+'/RandPeakPermute/'+gwasID+'.'+cellType+'.'+nCopy+'.RandPeakPermuteRecord.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		nPermute=nCopy+'-'+each[0]
		nPeak=each[1]
		newLine=[nPermute, nPeak]
		output1.write('\t'.join(newLine)+'\n')
	    input1.close()
	output1.close()
	    





