#!/home/ajl1213/anaconda2/bin/python


import os
import sys 
import numpy
import random



DIR=os.getcwd()
ncpu=1
fai='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
LdSnpDIR=DIR+'/LdFinal'
PeakDIR='/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/PeakBED'
cellType='MergedEC'
maxIter=10000

gwasList=[]
input1=open(DIR+'/GwasList.txt','r')
all_input1=input1.readlines()
for line in all_input1[1:]:
    each=line.strip().split('\t')
    gwasID=each[0]
    gwasList.append(gwasID)
input1.close()


os.system('mkdir '+DIR+'/InputPeakBED')
os.system('mkdir '+DIR+'/RandPeakPermute')



def makeInputPeakBED():
    runLine = ['cat ']
    for i in ['EarlyEC','MidEC','LateEC','FullEC']:
	runLine.append(PeakDIR+'/'+i+'.bed ')
    runLine.append(' | sortBed -i stdin > '+DIR+'/InputPeakBED/'+cellType+'.bed')
    os.system(''.join(runLine))


def makeScript(gwasID, cellType, fai, LdSnpDIR, maxIter):
    output1=open('PBS_randPeakPermute.'+gwasID+'.'+cellType+'.pbs','w')
    output1.write('#PBS -N randPeakPermute.'+gwasID+'.'+cellType+'\n')
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

    output1.write('python '+DIR+'/4b_permuteDiffPeak.py '+gwasID+' '+cellType+' '+fai+' '+LdSnpDIR+' '+DIR+'/InputPeakBED '+str(maxIter)+'\n')

    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')

    output1.close()



makeInputPeakBED()

for gwasID in gwasList:
    makeScript(gwasID, cellType, fai, LdSnpDIR, maxIter)
    os.system('qsub PBS_randPeakPermute.'+gwasID+'.'+cellType+'.pbs')




