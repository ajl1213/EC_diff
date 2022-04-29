#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
ncpu=1
fai='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
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

maxCopy=10
maxIter=100000

os.system('mkdir '+DIR+'/InputPeakBED')
os.system('mkdir '+DIR+'/RandPeakPermute')



def makeScript():
    for cellType in cellTypeList:
        for gwasID in gwasList:
	    for j in range(maxCopy):
		nCopy=str(j+1)

		output1=open('PBS_randPeakPermute.'+gwasID+'.'+cellType+'.'+nCopy+'.pbs','w')
		output1.write('#PBS -N randPeakPermute.'+gwasID+'.'+cellType+'.'+nCopy+'\n')
		output1.write('#PBS -q workq\n')
		output1.write('#PBS -l nodes=1:ppn='+str(ncpu)+'\n')
		output1.write('#PBS -j oe\n')
		output1.write('\n')
		output1.write('# go workdir\n')
		output1.write('cd $PBS_O_WORKDIR\n')
		output1.write('\n')
		output1.write('# run command \n')
		output1.write('sleep 5\n')
		output1.write('\n')
		output1.write('echo -n \"I am on: \"\n')
		output1.write('hostname;\n')
		output1.write('echo finding ssh-accessible nodes:\n')
		output1.write('echo -n \"running on: \"\n')
		output1.write('\n')

		output1.write('python '+DIR+'/1b_permuteDiffPeak.py '+gwasID+' '+cellType+' '+fai+' '+LdSnpDIR+' '+PeakDIR+' '+str(maxIter)+' '+nCopy+'\n')

		output1.write('\n')
		output1.write('sleep 30\n')
		output1.write('exit 0')

		output1.close()

		os.system('qsub PBS_randPeakPermute.'+gwasID+'.'+cellType+'.'+nCopy+'.pbs')


makeScript()




