#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
inputDIR='/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/PeakBED'
LdSnpDIR=DIR+'/LdFinal'
PeakDIR=DIR+'/InputPeakBED'
cellType='MergedEC'
gwasList=[
'Colorectal_cancer',
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
'Ulcerative_colitis',
'Rheumatoid_arthritis'
]


os.system('mkdir '+DIR+'/InputPeakBED')
os.system('mkdir '+DIR+'/ObservedData')



def makeInputPeakBED():
    runLine = ['cat ']
    for i in ['EarlyEC','MidEC','LateEC','FullEC']:
        runLine.append(inputDIR+'/'+i+'.bed ')
    runLine.append(' | sortBed -i stdin > '+DIR+'/InputPeakBED/'+cellType+'.bed')
    os.system(''.join(runLine))


def getActualData():

    for gwasID in gwasList:

	runScript='intersectBed -a '+PeakDIR+'/'+cellType+'.bed -b '+LdSnpDIR+'/'+gwasID+'.bed -wa -wb | cut -f1,2,3 | sort -u | wc -l > '+DIR+'/ObservedData/'+gwasID+'.'+cellType+'.nPeakIntersect.txt'
	print runScript
	os.system(runScript)

	runScript='intersectBed -a '+PeakDIR+'/'+cellType+'.bed -b '+LdSnpDIR+'/'+gwasID+'.bed -wa -wb | cut -f1,2,3 | sort -u > '+DIR+'/ObservedData/'+gwasID+'.'+cellType+'.Peak.intersect.bed'
	print runScript
	os.system(runScript)

	runScript='intersectBed -a '+PeakDIR+'/'+cellType+'.bed -b '+LdSnpDIR+'/'+gwasID+'.bed -wa -wb | sort -u > '+DIR+'/ObservedData/'+gwasID+'.'+cellType+'.Peak.intersect.full.bed'
	print runScript
	os.system(runScript)



makeInputPeakBED()
getActualData()



