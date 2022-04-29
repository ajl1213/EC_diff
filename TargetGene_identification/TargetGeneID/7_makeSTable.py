#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
rawSnpDIR='/home/ajl1213/Projects/Endothelial/Analysis2/GwasTarget/01_MergeSNP/RawGWAS'
sumStatDIR='/home/ajl1213/Projects/Endothelial/Analysis2/GwasEnrich/3_LDSC/GWAS_sumstat'

inputDict1={
'Colorectal_cancer':['Colorectal_cancer','Colorectal_adenoma_advanced','Colorectal_cancer_or_advanced_adenoma']
}
inputDict2={
'CRC_Jiang_2021':['Nat_Genet','2021'],
'CRC_Rashkin_2020':['Nat Commun', '2020'],
'CRC_Zhou_2018':['Nat_Genet', '2018']
}

gwasList=[
'Crohns_disease',
'Asthma',
'Inflammatory_bowel_disease',
'Multiple_sclerosis',
'Eczema',
'Atrial_fibrillation',
'Atopic_asthma',
'Coronary_artery_disease',
'Varicose_veins',
'Alzheimers_disease_late_onset',
'Ulcerative_colitis',
'Rheumatoid_arthritis'
]


pvalCut=5e-8



def makeSTable():

    output1=open(DIR+'/STable.GwasTarget.txt','w')
    output1.write('GWAS_ID\tPeakID\tPeakLabel\tEnsembleID\tGeneID\tInter_type\ttagSNP\tSNP_info\tABC_score\n')

    ## CRC data
    snpDict={}
    for gwasID in inputDict1:
	for i in inputDict1[gwasID]:
	    input1=open(rawSnpDIR+'/'+i+'.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		rsID=each[0]
		chrID=each[1]
		pos=each[2]
		pval=each[3]
		trait=each[4]
		journal=each[5].replace(' ','_')
		date=each[6]

		snpID=chrID+'.'+pos
		valID=rsID+','+pval+','+journal+','+date

		if snpDict.has_key(snpID):
		    tmp=snpDict[snpID]
		    tmp.append(valID)
		    snpDict[snpID]=tmp
		else:
		    snpDict[snpID]=[valID]
	    input1.close()

    for gwasID in inputDict2:
	input1=open(sumStatDIR+'/'+gwasID+'/sumStat.txt', 'r')
	all_input1=input1.readlines()

	idx=0
	for j in all_input1[0].strip().split('\t'):
	    if j == 'P':
		pval_idx=idx
	    if j == 'SNP':
		rsID_idx=idx
	    idx+=1

	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    chrID='chr'+each[0]
	    snpPos=chrID+'.'+each[1]
	    pval=each[pval_idx]
	    rsID=each[rsID_idx]
	    journal=inputDict2[gwasID][0]
	    year=inputDict2[gwasID][1]

	    valID=rsID+','+pval+','+journal+','+year

	    if float(pval) < float(pvalCut):
		if snpDict.has_key(snpPos):
		    tmp=snpDict[snpPos]
		    tmp.append(valID)
		    snpDict[snpPos]=tmp
		else:
		    snpDict[snpPos]=[valID]
	input1.close()

    peakSnpDict={}
    input1=open(DIR+'/RESULT/Colorectal_cancer.SnpToPeak.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	tagSNP=each[0]
	snpLabel=each[1]
	peakID=each[2]
	peakLabel=each[3]
	if peakSnpDict.has_key(peakID):
	    tmp=peakSnpDict[peakID]
	    tmp.append(tagSNP)
	    peakSnpDict[peakID]=tmp
	else:
	    peakSnpDict[peakID]=[tagSNP]
    input1.close()

    #
    input1=open(DIR+'/TargetGeneList/Colorectal_cancer.MergedTarget.abcScore.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[peakID, peakLabel, ensembleID, geneID, assoType, abcScore]=line.strip().split('\t')

	snpInfoList=[]
	for tagSNP in peakSnpDict[peakID]:
	    for i in snpDict[tagSNP]:
		snpInfoList.append(i)

	newLine=['Colorectal_cancer', peakID, peakLabel, ensembleID, geneID, assoType, ';'.join(set(peakSnpDict[peakID])), ';'.join(set(snpInfoList)), abcScore]
	output1.write('\t'.join(newLine)+'\n')

    input1.close()

    ## Rest of GWAS diseases
    for gwasID in gwasList:
	print gwasID
	# load SNP info
	snpDict={}
	input1=open(rawSnpDIR+'/'+gwasID+'.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    rsID=each[0]
	    chrID=each[1]
	    pos=each[2]
	    pval=each[3]
	    trait=each[4]
	    journal=each[5].replace(' ','_')
	    date=each[6]

	    snpID=chrID+'.'+pos
	    valID=rsID+','+pval+','+journal+','+date

	    if snpDict.has_key(snpID):
		tmp=snpDict[snpID]
		tmp.append(valID)
		snpDict[snpID]=tmp
	    else:
		snpDict[snpID]=[valID]
	input1.close()

	# load SNP-cRE info
	peakSnpDict={}
	input1=open(DIR+'/RESULT/'+gwasID+'.SnpToPeak.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    tagSNP=each[0]
	    snpLabel=each[1]
	    peakID=each[2]
	    peakLabel=each[3]
	    if peakSnpDict.has_key(peakID):
		tmp=peakSnpDict[peakID]
		tmp.append(tagSNP)
		peakSnpDict[peakID]=tmp
	    else:
		peakSnpDict[peakID]=[tagSNP]
	input1.close()

	#
	input1=open(DIR+'/TargetGeneList/'+gwasID+'.MergedTarget.abcScore.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    [peakID, peakLabel, ensembleID, geneID, assoType, abcScore]=line.strip().split('\t')

	    snpInfoList=[]
	    for tagSNP in peakSnpDict[peakID]:
		for i in snpDict[tagSNP]:
		    snpInfoList.append(i)

	    newLine=[gwasID, peakID, peakLabel, ensembleID, geneID, assoType, ';'.join(set(peakSnpDict[peakID])), ';'.join(set(snpInfoList)), abcScore]
	    output1.write('\t'.join(newLine)+'\n')

	input1.close()
    output1.close()


def reportStat():

    #
    allDict={}
    input1=open(DIR+'/STable.GwasTarget.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	gwasID=each[0]
	ensembleID=each[3]
	geneID=each[4]
	if allDict.has_key(gwasID):
	    tmp=allDict[gwasID]
	    if not geneID in tmp:
		tmp.append(geneID)
	    allDict[gwasID]=tmp
	else:
	    allDict[gwasID]=[geneID]
    input1.close()


    nTotal=0
    for gwasID in allDict:
	nTarget=len(allDict[gwasID])
	print gwasID, str(nTarget)

	nTotal+=nTarget

    print 'Total target gene #: '+str(nTotal)

    #
    totalGeneCount=0
    diffGeneCount=0
    for gwasID in allDict:
	tmpDict={}
	input1=open(DIR+'/STable.GwasTarget.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    if gwasID ==each[0]:
		gwasID=each[0]
		peakLabel=each[2]
		ensembleID=each[3]
		geneID=each[4]
		if tmpDict.has_key(geneID):
		    tmp=tmpDict[geneID]
		    if not peakLabel in tmp:
			tmp.append(peakLabel)
		    tmpDict[geneID]=tmp
		else:
		    tmpDict[geneID]=[peakLabel]
	input1.close()

	for geneID in tmpDict:
	    if not (len(tmpDict[geneID])==1 and tmpDict[geneID][0] == 'FullEC'):
		diffGeneCount+=1
	    totalGeneCount+=1

    print 'totalGeneCount: ' +str(totalGeneCount)
    print 'differentiationGeneCount: ' +str(diffGeneCount)
    print 'diffGeneFrac:' + str(float(diffGeneCount)/float(totalGeneCount))




#makeSTable()
reportStat()




