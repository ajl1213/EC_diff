#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
inputDIR=DIR+'/LdExpand'

gwasList=[
'Colorectal_cancer',
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

os.system('mkdir '+DIR+'/LdFinal')


def PostProcess():

    for gwasID in gwasList:
	print gwasID

	snpDict={}
	for i in range(23):
	    if i==22:
		chrID='chrX'
	    else:
		chrID='chr'+str(i+1)

	    input1=open(DIR+'/LdExpand/'+gwasID+'.'+chrID+'.LD.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1:
		each=line.strip().split('\t')
		rsID=each[0].replace('.','_')
		chrID=each[1]
		parentalPos=each[2]
		ldSnpList=each[3]
		snpDict[chrID+'.'+parentalPos]=[rsID+'.'+chrID+'.'+parentalPos]

		if not ldSnpList=='no':
		    for i in ldSnpList.split(',')[:-1]:
			if int(i.split(':')[1]) > 1:   ## Get 2, 3, 4 or 5
			    ldSnp=chrID+'.'+i.split(':')[0]
			    if snpDict.has_key(ldSnp):
				tmp=snpDict[ldSnp]
				tmp.append(rsID+'.'+chrID+'.'+parentalPos)
				snpDict[ldSnp]=tmp
			    else:
				snpDict[ldSnp]=[rsID+'.'+chrID+'.'+parentalPos]
	    input1.close()
		    
	output1=open(DIR+'/LdFinal/'+gwasID+'.allchr.txt','w')
	output1.write('LdSnpPos\tParentalSnpPos\tParantalSnpID\n')
	for snpID in snpDict:
	    for i in snpDict[snpID]:
		chrID=i.split('.')[1]
		parentalSnpID=i.split('.')[0]
		parentalSnpPos=i.split('.')[2]
		parentalSNP=chrID+'.'+parentalSnpPos

		newLine=[snpID, parentalSNP, parentalSnpID]
		output1.write('\t'.join(newLine)+'\n')
	output1.close()
	    

def getBedFiles():
    for gwasID in gwasList:
	output1=open(DIR+'/LdFinal/'+gwasID+'.bed','w')
	input1=open(DIR+'/LdFinal/'+gwasID+'.allchr.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    if len(each)==3:
		LdSnpPos = each[0]
		ParentalSnpPos = each[1]
		ParentalSnpID = each[2]

		chrID=LdSnpPos.split('.')[0]
		pt1=LdSnpPos.split('.')[1]
		pt2=str(int(pt1)+1)
		newLine=[chrID, pt1, pt2, ParentalSnpPos, ParentalSnpID]
		output1.write('\t'.join(newLine)+'\n')
	input1.close()
	output1.close()




PostProcess()
getBedFiles()




