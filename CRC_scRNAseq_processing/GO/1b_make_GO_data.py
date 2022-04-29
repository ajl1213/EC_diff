#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import operator



DIR=os.getcwd()
inputDIR=DIR+'/GOBP'
dataTypeList=['snDEGs','trajGenes']
cellTypeList=['EarlyEC','MidEC','LateEC','FullEC']

selectedTerms=[
'peptide_biosynthetic_process',
'regulation_of_apoptotic_process',
'DNA_damage_response,_signal_transduction_by_p53_class_mediator',
'extracellular_matrix_organization',
'neutrophil_mediated_immunity',
'endoderm_formation',
'positive_regulation_of_endothelial_cell_proliferation',
'extracellular_matrix_disassembly',
'regulation_of_angiogenesis',
'endothelial_cell_migration',
'positive_regulation_of_inflammatory_response'
]


def load_CRC_GOBP():
    crcAllDict={}
    crcPvalDict={}
    crcTermList=[]
    for dataType in dataTypeList:
	tmp1Dict={}
	tmp2Dict={}
	input1=open(DIR+'/GOBP/'+dataType+'.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    termID=each[0].replace('\'','')
	    termID='_'.join(termID.split(' ')[:-1])
	    if not termID in crcTermList:
		crcTermList.append(termID)

	    goID=each[0].replace('\'','').split(' ')[-1]
	    goID=goID.replace('(','')
	    goID=goID.replace(')','')

	    pval=each[2]
	    adjPval=each[3]
	    oddRatio=each[6]
	    combinedScore=each[7]
	    genes=each[8]

	    tmp1Dict[termID]=[pval, oddRatio, combinedScore, genes]
	    tmp2Dict[termID]=-numpy.log10(float(pval))

	input1.close()

	crcAllDict[dataType]=tmp1Dict
	crcPvalDict[dataType]=tmp2Dict

    return crcTermList, crcAllDict, crcPvalDict


def load_stageDEG_GOBP():
    stageAllDict={}
    stagePvalDict={}
    stageTermList=[]
    for cellType in cellTypeList:
	tmp1Dict={}
	tmp2Dict={}
	input1=open(DIR+'/GOBP/'+cellType+'.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    termID=each[0].replace('\'','')
	    termID='_'.join(termID.split(' ')[:-1])
	    if not termID in stageTermList:
		stageTermList.append(termID)

	    goID=each[0].replace('\'','').split(' ')[-1]
	    goID=goID.replace('(','')
	    goID=goID.replace(')','')

	    pval=each[2]
	    adjPval=each[3]
	    oddRatio=each[6]
	    combinedScore=each[7]
	    genes=each[8]

	    tmp1Dict[termID]=[pval, oddRatio, combinedScore, genes]
	    tmp2Dict[termID]=-numpy.log10(float(pval))

	input1.close()

	stageAllDict[cellType]=tmp1Dict
	stagePvalDict[cellType]=tmp2Dict

    return stageTermList, stageAllDict, stagePvalDict


def filterGOBP():
    crcTermList, crcAllDict, crcPvalDict = load_CRC_GOBP()
    stageTermList, stageAllDict, stagePvalDict = load_stageDEG_GOBP()

    output1=open(DIR+'/StageGOBP.merged.txt','w')
    output1.write('CellStage\tGOBP\t'+'\t'.join(dataTypeList)+'\t'+'\t'.join(cellTypeList)+'\tcommonGeneSet\n')
    for cellType in cellTypeList:
	print cellType
	for key in sorted(stagePvalDict[cellType].items(), key=operator.itemgetter(1), reverse=True):
	    termID=key[0]
	    logPval=key[1]
	    if logPval > 3:
		if (crcPvalDict['snDEGs'].has_key(termID) and crcPvalDict['trajGenes'].has_key(termID)):
		    commonGeneSet=set(stageAllDict[cellType][termID][3].split(';')) & set(crcAllDict['snDEGs'][termID][3].split(';')) & set(crcAllDict['trajGenes'][termID][3].split(';'))
		    if len(commonGeneSet) > 0:

			pvalList=[]
			pvalList.append(str(crcPvalDict['snDEGs'][termID]))
			pvalList.append(str(crcPvalDict['trajGenes'][termID]))
			for cellID in cellTypeList:
			    if stagePvalDict[cellID].has_key(termID):
				pvalList.append(str(stagePvalDict[cellID][termID]))
			    else:
				pvalList.append('0.0')
			output1.write(cellType+'\t'+termID+'\t'+'\t'.join(pvalList)+'\t'+';'.join(commonGeneSet)+'\n')
    output1.close()


def makeData():

    input1=open(DIR+'/StageGOBP.merged.txt','r')
    output1=open(DIR+'/StageGOBP.merged.selected.txt','w')
    output1.write('GOBP\t'+'\t'.join(dataTypeList)+'\t'+'\t'.join(cellTypeList)+'\n')
    all_input1=input1.readlines()
    tmpList=[]
    geneDict={}
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	stage=each[0]
	termID=each[1]
	geneList=each[-1]
	if termID in selectedTerms:
	    #
	    if not termID in tmpList:
		output1.write('\t'.join(each[1:-1])+'\n')
		tmpList.append(termID)

	    #
	    if stage=='LateEC':
		for geneID in geneList.split(';'):
		    if geneDict.has_key(termID):
			tmp=geneDict[termID]
			if not geneID in tmp:
			    tmp.append(geneID)
			geneDict[termID]=tmp
		    else:
			geneDict[termID]=[geneID]

    input1.close()
    output1.close()

    return geneDict


def abc():
    geneDict=makeData()

    specificGeneExpression=[
    'ADAMTS4',
    'ADAMTS9',
    'APLN',
    'APLNR',
    'COL4A1',
    'COL12A1',
    'ITGAV',
    'NRP2',
    'P4HB',
    'PHLDA1',
    'PLAUR',
    'SERPINH1',
    'SOX17',
    'SPHK1',
    'TIMP1',
    'TNFRSF10D'
    ]

    #
    finalDict={}
    for termID in geneDict:
	tmpList=[]
	for geneID in geneDict[termID]:
	    if geneID in specificGeneExpression:
		tmpList.append(geneID)
	finalDict[termID]=tmpList	


    #
    for termID in finalDict:
	print termID
	print finalDict[termID]



filterGOBP()
makeData()
abc()




