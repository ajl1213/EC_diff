#!/home/ajl1213/anaconda2/bin/python


import os



ncpu=2
DIR=os.getcwd()
inputDIR=DIR+'/rawGWAS'
LD_DIR='/home/ajl1213/genome.info/LD/LD_naoki'

gwasList=[]
input1=open(DIR+'/GwasList.txt','r')
all_input1=input1.readlines()
for line in all_input1[1:]:
    each=line.strip().split('\t')
    gwasID=each[0]
    gwasList.append(gwasID)
input1.close()

os.system('mkdir '+DIR+'/LdExpand')



def makeScript(inputDIR, LD_DIR, gwasID, chrID):
    output1=open('PBS_ld_'+gwasID+'_'+chrID+'.pbs','w')
    output1.write('#PBS -N ld_'+gwasID+'_'+chrID+'\n')
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

    output1.write('python '+DIR+'/2b_LdExpand.py '+inputDIR+' '+LD_DIR+' '+gwasID+' '+chrID+'\n')

    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')
    output1.close()



for gwasID in gwasList:
    for i in range(23):
	if i==22:
	    chrID='chrX'
	else:
	    chrID='chr'+str(i+1)
	makeScript(inputDIR, LD_DIR, gwasID, chrID)
	os.system('qsub PBS_ld_'+gwasID+'_'+chrID+'.pbs')



