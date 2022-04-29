#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
inputGwasCatalog='/home/ajl1213/genome.info/GWAS_catalog_downloaded_2022.01.13/GWAS_catalog.truncated.hg19.liftover.final.txt'
nLimit=100
pvalCut=5e-8

os.system('mkdir '+DIR+'/rawGWAS')



def getSNPs():

    sizeDict={}
    input1=open(inputGwasCatalog,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	rsID=each[0]
        chrID=each[1]
        pt=each[2]
	try:
	    pval=float(each[3])
	    if pval < pvalCut:
		trait=each[4].replace(' ','_')
		trait=trait.replace('/','_')
		trait=trait.replace('\'','')
		trait=trait.replace('(','')
		trait=trait.replace(')','')
		journal=each[5]
		date=each[6]
		if sizeDict.has_key(trait):
		    sizeDict[trait]+=1
		else:
		    sizeDict[trait]=1
	except:
	    print each[4], rsID, chrID, pt, pval

    input1.close()

    gwasList=[]
    output1=open(DIR+'/GwasList.txt','w')
    output1.write('gwasID\tnSNP\n')
    for trait in sizeDict:
	if sizeDict[trait] > nLimit:
	    gwasList.append(trait)
	    newLine=[trait, str(sizeDict[trait])]
	    output1.write('\t'.join(newLine)+'\n')
    output1.close()

    print len(gwasList)
    idx=1
    for gwasID in gwasList:
	print float(idx)/float(len(gwasList))*100
	idx+=1

	output1=open(DIR+'/rawGWAS/'+gwasID+'.txt','w')
	output1.write('rsID\tchrID\tpos\tpval\ttrait\tjournal\tdate\n')
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    rsID=each[0]
	    chrID=each[1]
	    pt=each[2]
#	    try:
	    pval=each[3]
#		if float(pval) < pvalCut:
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
#	    except:
#		pass

	output1.close()



getSNPs()



