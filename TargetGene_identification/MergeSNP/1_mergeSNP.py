#!/home/ajl1213/anaconda2/bin/python
# merge SNPs from GWAS catalogue and individual GWAS summary statistics


import os
import numpy


DIR=os.getcwd()
inputGwasCatalog='/home/ajl1213/genome.info/GWAS_catalog_downloaded_2022.01.13/GWAS_catalog.truncated.hg19.liftover.final.txt'
SumStatDIR='/home/ajl1213/Projects/Endothelial/Analysis2/GwasEnrich/3_LDSC/GWAS_sumstat'

inputDict1={
'Colorectal_cancer':['Colorectal_cancer','Colorectal_adenoma_advanced','Colorectal_cancer_or_advanced_adenoma']
}
inputDict2={
'Colorectal_cancer':['CRC_Jiang_2021', 'CRC_Rashkin_2020', 'CRC_Zhou_2018']
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

os.system('mkdir '+DIR+'/RawGWAS')
os.system('mkdir '+DIR+'/MergedRawSNP')



def getColorectalCancerSNPs():
    #
    input1=open(inputGwasCatalog,'r')
    all_input1=input1.readlines()
    for i in inputDict1:
	for gwasID in inputDict1[i]:
	    output1=open(DIR+'/RawGWAS/'+gwasID+'.txt','w')
	    output1.write('rsID\tchrID\tpos\tpval\ttrait\tjournal\tdate\n')
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		rsID=each[0]
		chrID=each[1]
		pt=each[2]
		pval=each[3]
		trait=each[4].replace(' ','_')
		trait=trait.replace('/','_')
		trait=trait.replace('\'','')
		trait=trait.replace('(','')
		trait=trait.replace(')','')
		journal=each[5]
		date=each[6]
		if gwasID == trait:
		    newLine=[rsID, chrID, pt, pval, trait, journal, date]
		    output1.write('\t'.join(newLine)+'\n')
	    output1.close()
    input1.close()

    allDict={}
    # Catalog
    for i in inputDict1:
	for gwasID in inputDict1[i]:
	    input1=open(DIR+'/RawGWAS/'+gwasID+'.txt', 'r')
	    all_input1=input1.readlines()
	    for line in all_input1[1:]:
		[rsID, chrID, pos, pval, trait, journal, date]=line.strip().split('\t')

		journal=journal.replace(' ','_')
		valID=gwasID+'_'+journal+'_'+date
		snpPos=chrID+'.'+pos
		if allDict.has_key(snpPos):
		    tmp=allDict[snpPos]
		    tmp.append(valID)
		    allDict[snpPos]=tmp
		else:
		    allDict[snpPos]=[valID]
	    input1.close()

    # Summary statistics
    for i in inputDict2:
	for gwasID in inputDict2[i]:
	    input1=open(SumStatDIR+'/'+gwasID+'/sumStat.txt', 'r')
	    all_input1=input1.readlines()

	    idx=0
	    for j in all_input1[0].strip().split('\t'):
		if j == 'P':
		    pval_idx=idx
		idx+=1

	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		chrID='chr'+each[0]
		snpPos=chrID+'.'+each[1]
		pval=each[pval_idx]

		if float(pval) < float(pvalCut):
		    if allDict.has_key(snpPos):
			tmp=allDict[snpPos]
			tmp.append(gwasID)
			allDict[snpPos]=tmp
		    else:
			allDict[snpPos]=[gwasID]
	    input1.close()

    # make list
    for i in inputDict1:
	output1=open(DIR+'/MergedRawSNP/'+i+'.rawSNP.txt','w')
	output1.write('chrID\tpos\tjournalList\n')
	for snpPos in allDict:
	    chrID=snpPos.split('.')[0]
	    pt=snpPos.split('.')[1]
	    jList=set(allDict[snpPos])
	    newLine=[chrID, pt, ';'.join(jList)]
	    output1.write('\t'.join(newLine)+'\n')
	output1.close() 


def getOtherSNPs():
    input1=open(inputGwasCatalog,'r')
    all_input1=input1.readlines()
    for gwasID in gwasList:

        output1=open(DIR+'/RawGWAS/'+gwasID+'.txt','w')
        output1.write('rsID\tchrID\tpos\tpval\ttrait\tjournal\tdate\n')
        for line in all_input1[1:]:
            each=line.strip().split('\t')
            rsID=each[0]
            chrID=each[1]
            pt=each[2]
            pval=each[3]
            trait=each[4].replace(' ','_')
            trait=trait.replace('/','_')
            trait=trait.replace('\'','')
            trait=trait.replace('(','')
            trait=trait.replace(')','')
            journal=each[5]
            date=each[6]
            if gwasID == trait:
                newLine=[rsID, chrID, pt, pval, trait, journal, date]
                output1.write('\t'.join(newLine)+'\n')
        output1.close()

    # Catalog
    for gwasID in gwasList:
	allDict={}
	input1=open(DIR+'/RawGWAS/'+gwasID+'.txt', 'r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    [rsID, chrID, pos, pval, trait, journal, date]=line.strip().split('\t')

	    journal=journal.replace(' ','_')
	    valID=gwasID+'_'+journal+'_'+date
	    snpPos=chrID+'.'+pos
	    if allDict.has_key(snpPos):
		tmp=allDict[snpPos]
		tmp.append(valID)
		allDict[snpPos]=tmp
	    else:
		allDict[snpPos]=[valID]
	input1.close()

	output1=open(DIR+'/MergedRawSNP/'+gwasID+'.rawSNP.txt','w')
	output1.write('chrID\tpos\tjournalList\n')
	for snpPos in allDict:
	    chrID=snpPos.split('.')[0]
	    pt=snpPos.split('.')[1]
	    jList=set(allDict[snpPos])
	    newLine=[chrID, pt, ';'.join(jList)]
	    output1.write('\t'.join(newLine)+'\n')
	output1.close() 




getColorectalCancerSNPs()
getOtherSNPs()



