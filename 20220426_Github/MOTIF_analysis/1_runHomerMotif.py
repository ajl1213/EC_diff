#!/home/ajl1213/anaconda2/bin/python


import os



DIR=os.getcwd()
genomeVersion='hg19'
inputDIR='/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/PeakBED'
sampleList=['hESC', 'Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC', 'EarlyNonEC', 'MidNonEC', 'LateNonEC', 'FullNonEC']
os.system('mkdir '+DIR+'/HOMER')
cpu=8



def makeScript():
    for sampleID in sampleList:
	output1=open('PBS_Motif.'+sampleID+'.pbs','w')
	output1.write('#PBS -N Motif.'+sampleID+'\n')
	output1.write('#PBS -q workq\n')
	output1.write('#PBS -l nodes=1:ppn='+str(cpu)+'\n')
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
	output1.write('findMotifsGenome.pl '+inputDIR+'/'+sampleID+'.bed '+genomeVersion+' '+DIR+'/HOMER/'+sampleID+' -size given\n')
	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()
	os.system('qsub PBS_Motif.'+sampleID+'.pbs')


makeScript()




