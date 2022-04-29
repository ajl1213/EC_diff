#!/home/ajl1213/anaconda2/bin/python


import os



DIR=os.getcwd()
chromSizeFile='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
chrList=['chr1', 'chr2', 'chr3','chr4', 'chr5', 'chr6','chr7', 'chr8', 'chr9','chr10', 'chr11', 'chr12','chr13', 'chr14', 'chr15','chr16', 'chr17', 'chr18','chr19', 'chr20', 'chr21','chr22', 'chrX']

##
chromSizeDict = {}    
file1 = open(chromSizeFile, 'r')
for rawline in file1:
    splitted  = rawline.split("\t")
    chromSizeDict[splitted[0]] = int(splitted[1])
file1.close()

os.system('mkdir '+DIR+'/Bin')



def makeBin_5kb():
    outputall=open(DIR+'/Bin/allchr.5kb.bin','w')
    for chrID in chrList:
	output1=open(DIR+'/Bin/'+chrID+'.5kb.bin','w')
	i=0
	size=int(chromSizeDict[chrID])
	while i<size:
	    output1.write(chrID+'\t'+str(i)+'\t'+str(i+5000)+'\n')
	    outputall.write(chrID+'\t'+str(i)+'\t'+str(i+5000)+'\n')
	    i+=5000
	output1.close()
    outputall.close()


def makeBin_10kb():
    outputall=open(DIR+'/Bin/allchr.10kb.bin','w')
    for chrID in chrList:
	output1=open(DIR+'/Bin/'+chrID+'.10kb.bin','w')
	i=0 
	size=int(chromSizeDict[chrID])
	while i<size:
	    output1.write(chrID+'\t'+str(i)+'\t'+str(i+10000)+'\n')
	    outputall.write(chrID+'\t'+str(i)+'\t'+str(i+10000)+'\n')
	    i+=10000
	output1.close()
    outputall.close()


def makeBin_20kb():
    outputall=open(DIR+'/Bin/allchr.20kb.bin','w')
    for chrID in chrList:
	output1=open(DIR+'/Bin/'+chrID+'.20kb.bin','w')
	i=0 
	size=int(chromSizeDict[chrID])
	while i<size:
	    output1.write(chrID+'\t'+str(i)+'\t'+str(i+20000)+'\n')
	    outputall.write(chrID+'\t'+str(i)+'\t'+str(i+20000)+'\n')
	    i+=20000
	output1.close()
    outputall.close()


def makeBin_40kb():
    outputall=open(DIR+'/Bin/allchr.40kb.bin','w')
    for chrID in chrList:
	output1=open(DIR+'/Bin/'+chrID+'.40kb.bin','w')
	i=0 
	size=int(chromSizeDict[chrID])
	while i<size:
	    output1.write(chrID+'\t'+str(i)+'\t'+str(i+40000)+'\n')
	    outputall.write(chrID+'\t'+str(i)+'\t'+str(i+40000)+'\n')
	    i+=40000
	output1.close()
    outputall.close()




makeBin_5kb()
makeBin_10kb()
#makeBin_20kb()
#makeBin_40kb()



