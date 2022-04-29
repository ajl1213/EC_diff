#!/home/ajl1213/anaconda2/bin/python


import os


ncpu=4
DIR=os.getcwd()
inputDIR='/home/ajl1213/Projects/Endothelial/data/HiC/Mapping/SortedNodupBam'
FAI='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'

sampleList=[
'ME_rep1.HiC',
'ME_rep2.HiC',
'ME',
'EC12_rep1.HiC',
'EC12_rep2.HiC',
'EC12',
'EC48_rep1.HiC',
'EC48_rep2.HiC',
'EC48',
'MergedEC'
]

chrList=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',  'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16','chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'] 
#resolution='5kb'
resolution='10kb'
maxDist=1000000

os.system('mkdir '+DIR+'/CisOnlyBam')
os.system('mkdir '+DIR+'/Coverage')
os.system('mkdir '+DIR+'/FeatureVec')
os.system('mkdir '+DIR+'/RScript')
os.system('mkdir '+DIR+'/Figs')
os.system('mkdir '+DIR+'/CovNorm')
os.system('mkdir '+DIR+'/DistNorm')



def makeScript(sampleID, chrID):
    output1=open('PBS_covNorm.'+sampleID+'.'+chrID+'.'+resolution+'.pbs','w')
    output1.write('#PBS -N covNorm.'+sampleID+'.'+chrID+'.'+resolution+'\n')
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

    output1.write('samtools view '+inputDIR+'/'+sampleID+'.sorted.nodup.bam '+chrID+' | python 01b_filterCisRead.py | samtools view -Sb -t '+FAI+' -o - - > '+DIR+'/CisOnlyBam/'+sampleID+'.'+chrID+'.bam\n')
    output1.write('coverageBed -counts -sorted -a '+DIR+'/Bin/'+chrID+'.'+resolution+'.bin  -b '+DIR+'/CisOnlyBam/'+sampleID+'.'+chrID+'.bam > '+DIR+'/Coverage/'+sampleID+'.'+chrID+'.'+resolution+'.coverage\n')
    output1.write('samtools view '+DIR+'/CisOnlyBam/'+sampleID+'.'+chrID+'.bam | python 01c_makeFeature.py '+resolution+' '+str(maxDist)+' ' +DIR+'/Bin/'+chrID+'.'+resolution+'.bin '+DIR+'/Coverage/'+sampleID+'.'+chrID+'.'+resolution+'.coverage '+DIR+'/FeatureVec/'+sampleID+'.'+chrID+'.'+resolution+'.count.gz\n')
    R_argument1=' \'--args '+DIR+'/FeatureVec/'+sampleID+'.'+chrID+'.'+resolution+'.count.gz '+sampleID+' '+chrID+' '+resolution+' '+str(maxDist)+'\''
    output1.write('source activate r_env\n')
    output1.write('R CMD BATCH --no-save --no-restore '+R_argument1+' 01d_CovNorm.R '+DIR+'/RScript/covNorm.'+sampleID+'.'+chrID+'.'+resolution+'.Rout\n')

    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')
    output1.close()



for sampleID in sampleList:
    for chrID in chrList: 
	makeScript(sampleID, chrID)
	os.system('qsub PBS_covNorm.'+sampleID+'.'+chrID+'.'+resolution+'.pbs')





