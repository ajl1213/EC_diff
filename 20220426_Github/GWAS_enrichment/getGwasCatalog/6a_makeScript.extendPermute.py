#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
ncpu=1
cellType='MergedEC'
inFile=DIR+'/RESULT/'+cellType+'.DataList.txt'
LdSnpDIR=DIR+'/LdFinal'
fai='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
maxCopy=10
maxIter=100000



def selectGWAS(): 
    gwasList=[]
    input1=open(inFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[gwasID, nPeak, fcVal, empiricalPval]=line.strip().split('\t')
	if float(empiricalPval) == 0:
	    gwasList.append(gwasID)
    input1.close()

    return gwasList


def makeScript(gwasID, cellType, fai, LdSnpDIR, maxIter, nCopy):
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

    output1.write('python '+DIR+'/6b_permuteDiffPeak.py '+gwasID+' '+cellType+' '+fai+' '+LdSnpDIR+' '+DIR+'/InputPeakBED '+str(maxIter)+' '+nCopy+'\n')

    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')

    output1.close()


gwasList=selectGWAS()


for gwasID in gwasList:
    for j in range(maxCopy):
	nCopy=str(j+1)

	makeScript(gwasID, cellType, fai, LdSnpDIR, maxIter, nCopy)
	os.system('qsub PBS_randPeakPermute.'+gwasID+'.'+cellType+'.'+nCopy+'.pbs')




