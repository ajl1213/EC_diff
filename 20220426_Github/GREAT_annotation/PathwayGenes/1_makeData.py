#!/home/ajl1213/anaconda2/bin/python



import os
import glob



DIR=os.getcwd()
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
DegInfo='/home/ajl1213/Projects/Endothelial/data/RNA/Processing/DEGs/RNA.DegLabel.txt'
GOBP_DIR='/home/ajl1213/Projects/Endothelial/Analysis2/AnnoGREAT/1_AnnoGREAT/GREAT_GOBP'

queryList=[
'EarlyEC;artery_development',
'EarlyEC;apoptotic_process_involved_in_development',
'EarlyEC;heart_morphogenesis',

'MidEC;epithelium_development',
'MidEC;morphogenesis_of_an_epithelium',

'LateEC;vasculature_development',
'LateEC;angiogenesis',

'FullEC;vasculature_development',
'FullEC;angiogenesis',
'FullEC;regulation_of_angiogenesis'
]




def getGeneInfo():
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


def LoadDEGs():
    geneDict=getGeneInfo()

    DegDict={}
    input1=open(DegInfo,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=geneDict[ensembleID]
        stageType=each[3]
        if DegDict.has_key(stageType):
            tmp=DegDict[stageType]
            tmp.append(geneID)
            DegDict[stageType]=tmp
        else:
            DegDict[stageType]=[geneID]
    input1.close()

    return DegDict


def LoadPathwayGenes():
    pathwayDict={}
    for i in queryList:
	dataType=i.split(';')[0]
	traitID=i.split(';')[1]

	input1=open(GOBP_DIR+'/'+dataType+'.GOBP.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    geneList=each[9]
	    if each[2].replace(' ','_')==traitID:
		pathwayDict[i]=geneList

	input1.close()


    return pathwayDict


def makeData():
    DegDict = LoadDEGs()
    pathwayDict = LoadPathwayGenes()

    # matched gene list
    output1=open(DIR+'/MatchedGeneList.DEGs.GO.txt','w')
    output1.write('stageType\ttraitID\tmatchedGeneList\n')
    for i in queryList:
	dataType=i.split(';')[0]
	traitID=i.split(';')[1]

	geneList=pathwayDict[i].split(',')
	matchedList=[]
	for geneID in geneList:
	    if geneID in DegDict[dataType]:
		matchedList.append(geneID)

	output1.write(dataType+'\t'+traitID+'\t'+','.join(matchedList)+'\n')
    output1.close()

    # unique gene list
    tmpList=[]
    input1=open(DIR+'/MatchedGeneList.DEGs.GO.txt','r')
    all_input1=input1.readlines()
    
    output1=open(DIR+'/MatchedGeneList.unique.DEGs.GO.txt','w')
    output1.write('stageType\ttraitID\tuniqGeneList\n')

    for line in all_input1[1:]:
	each=line.strip().split('\t')
	stageType=each[0]
	traitID=each[1]
	geneList=each[2]

	matchedList=[]
	for geneID in geneList.split(','):
	    if not geneID in tmpList:
		matchedList.append(geneID)
		tmpList.append(geneID)

	output1.write(stageType+'\t'+traitID+'\t'+','.join(matchedList)+'\n')
    output1.close()





makeData()



