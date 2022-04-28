#!/python2



import os
import numpy



DIR=os.getcwd()
dataTypeList=['EC_cRE:EarlyEC-MidEC-LateEC-FullEC', 'nonEC_cRE:EarlyNonEC-MidNonEC-LateNonEC-FullNonEC']
inputDIR=DIR+'/GREAT_GeneAnno'

os.system('mkdir '+DIR+'/geneCategory')



def makeGeneList():
    for cellType in ['hESC', 'Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC', 'EarlyNonEC', 'MidNonEC', 'LateNonEC', 'FullNonEC']:
	tmpDict={}	
        input1=open(inputDIR+'/'+cellType+'.GeneAnno.txt', 'r')
        all_input1=input1.readlines()
        for line in all_input1:
            each=line.strip().split('\t')
	    chrID=each[0]
	    pt1=each[1]
	    pt2=each[2]
	    geneID=each[3]
	    dist=each[4]

	    peakID=chrID+':'+pt1+'-'+pt2
	    valID=peakID+';'+dist 

	    if tmpDict.has_key(geneID):
		tmp=tmpDict[geneID]
		tmp.append(valID)
		tmpDict[geneID]=tmp
	    else:
		tmpDict[geneID]=[valID]
	input1.close()

        output1=open(inputDIR+'/'+cellType+'.geneList.txt', 'w')
	output1.write('geneID\tpeakList\n')
	for geneID in tmpDict:
	    output1.write(geneID+'\t'+','.join(tmpDict[geneID])+'\n')
	output1.close()


def categorize():
    for i in dataTypeList:
        #
        dataType=i.split(':')[0]
        cellTypeList=i.split(':')[1]
        geneListDict={}
        for cellType in cellTypeList.split('-'):
            input1=open(inputDIR+'/'+cellType+'.geneList.txt', 'r')
            all_input1=input1.readlines()
            for line in all_input1[1:]:
                each=line.strip().split('\t')
		geneID=each[0]
		peakList=each[1]
		if geneListDict.has_key(geneID):
		    tmp=geneListDict[geneID]
		    tmp.append(cellType)
		    geneListDict[geneID]=tmp
		else:
		    geneListDict[geneID]=[cellType]
            input1.close()

        #
        output1=open(DIR+'/geneCategory/'+dataType+'.GreatAnno.txt','w')
        output1.write('geneID\tdataType\tnAsso\n')
        for geneID in geneListDict:
	    nAsso=len(geneListDict[geneID])
            output1.write(geneID+'\t'+';'.join(geneListDict[geneID])+'\t'+str(nAsso)+'\n')
        output1.close()



makeGeneList()
categorize()




