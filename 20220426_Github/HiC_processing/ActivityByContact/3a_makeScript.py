#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
GeneCountFile='/home/ajl1213/Projects/Endothelial/data/RNA/Processing/CountMatrix/RNA.TPM.qq.txt'
H3K27acCountFile='/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/CountMatrix/H3K27ac.RPM.qq.txt'
inputFile=DIR+'/CreListPerGene/CreListPerGene.txt'
nGroup=100

os.system('mkdir '+DIR+'/PCC')


def loadGeneInfo():
    geneDict={}
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        geneDict[ensembleID]=geneID
    input1.close()

    return geneDict


def separateAsso():
    geneDict=loadGeneInfo()

    #
    elementList=[]
    input1=open(inputFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        creID=each[3]

        assoID=ensembleID+';'+creID
        elementList.append(assoID)
    input1.close()

    #

    sepList=numpy.array_split(elementList, nGroup)
    idx=1
    for i in sepList:
	print str(idx)
	output1=open(DIR+'/PCC/AssoList.'+str(idx)+'.txt','w')
	output1.write('ensembleID\tgeneID\tcreID\n')
	idx+=1
	for j in i:
	    ensembleID=j.split(';')[0]
	    geneID=geneDict[ensembleID]
	    creID=j.split(';')[1]
    
	    newLine=[ensembleID, geneID, creID]
	    output1.write('\t'.join(newLine)+'\n')

	output1.close()


def makeScript():
    ncpu=1

    for idx in range(nGroup):
	idx=idx+1

	output1=open('PBS_getPCC_'+str(idx)+'.pbs','w')
	output1.write('#PBS -N getPCC_'+str(idx)+'\n')
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

	output1.write('python '+DIR+'/3b_getPCC.py '+gtf+' '+GeneCountFile+' '+H3K27acCountFile+' ' +DIR+'/PCC/AssoList.'+str(idx)+'.txt '+str(idx)+'\n')

	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')

	output1.close()

	os.system('qsub PBS_getPCC_'+str(idx)+'.pbs')



separateAsso()
makeScript()




