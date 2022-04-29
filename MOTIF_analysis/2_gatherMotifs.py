#!/home/ajl1213/anaconda2/bin/python



import os
import operator



DIR=os.getcwd()
sampleList=['hESC', 'Mesoderm', 'EarlyEC', 'MidEC', 'LateEC', 'FullEC', 'EarlyNonEC', 'MidNonEC', 'LateNonEC', 'FullNonEC']

os.system('mkdir '+DIR+'/RESULT')



def getMotifCount():
    countDict={}
    for sampleID in sampleList:
	input1=open(DIR+'/HOMER/'+sampleID+'/homerResults.html', 'r')
	all_input1=input1.readlines()
	motifCount=0
	for line in all_input1[1:]:
	    if line.count('homerResults')>0:
		line=line.replace('>','').replace('<','')
		each=line.strip().split('TD')
		pval=float(each[2].replace('/',''))
		motifCount+=1
		countDict[sampleID]=motifCount
	input1.close()

    return countDict


def getMotifInfo():
    countDict=getMotifCount()

    scoreDict={}
    for sampleID in sampleList:
	pvalDict={}
	targetFracDict={}
	tfDict={}
	output1=open(DIR+'/RESULT/'+sampleID+'.TfList.txt','w')
	output1.write('MotifID\tPval\tTargetFrac\tTfList\n')
	for index in range(1, countDict[sampleID]+1):
	    keyID='Motif'+str(index)
	    input1=open(DIR+'/HOMER/'+sampleID+'/homerResults/motif'+str(index)+'.info.html','r')
	    all_input1=input1.readlines()
	    TfList=[]
	    nTF=0
	    for line in all_input1:
		if (line.count('p-value')>0 and line.count('log')==0):
		    pval=line.strip().split('<TD>')[2].split('<')[0]
		    pvalDict[keyID]=float(pval)

		if line.count('Percentage of Target Sequences')>0:
		    targetFrac=line.strip().split('TD>')[3].split('<')[0]
		    targetFracDict[keyID]=targetFrac

		if line.count('<H4>')>0:
		    TF=line.strip().split('<H4>')[1].split('/')[0].upper()
		    if TF.count('(')>0:
			TF=TF.split('(')[0]
		    TfList.append(TF)

		if line.count('Score')>0:
		    score=float(line.strip().split('<TD>')[2].split('<')[0])
		    if score>0.70:  ##
			scoreDict[sampleID+';'+TF]=score
			nTF+=1
	    input1.close()

	    TfFinal=[]
	    for TF in TfList[:nTF]:
		if not TF in TfFinal:
		    TfFinal.append(TF)

	    tfDict[keyID]=TfFinal
	    
	for key in sorted(pvalDict.items(), key=operator.itemgetter(1)):
	    motifID=key[0]
	    pval=str(key[1])
	    targetFrac=targetFracDict[motifID]
	    tfList=tfDict[motifID]
	    if len(tfList)>0:
		output1.write(motifID+'\t'+pval+'\t'+targetFrac+'\t'+','.join(tfList)+'\n')
	output1.close()

    ##
    totalList=[]
    for sampleID in sampleList:
	input1=open(DIR+'/RESULT/'+sampleID+'.TfList.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')

	    pval=each[1]
	    targetFrac=each[2]
	    tf_list=each[3].split(',')

	    for tf in tf_list:
		valID=[sampleID, tf, pval, str(scoreDict[sampleID+';'+tf]), targetFrac]
		totalList.append(valID)
	input1.close()

    output1=open(DIR+'/tflist.unique.txt','w')
    output1.write('dataType\tTF\tPval\tMotifScore\tTargetFrac\n')
    for i in totalList:
	output1.write('\t'.join(i)+'\n')
    output1.close()



def modifyName():

    ## Manual curation of TF ID
    conversionDict={}
    input1=open(DIR+'/STable.TF_conversion.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	homer_tf=each[2]
	geneID=each[3]
	conversionDict[homer_tf]=geneID
    input1.close()

    ##
    for sampleID in sampleList:
	input1=open(DIR+'/RESULT/'+sampleID+'.TfList.txt','r')
	all_input1=input1.readlines()
	output1=open(DIR+'/RESULT/'+sampleID+'.TfList.final.txt','w')
	output1.write('MotifID\tPval\tTargetFrac\tTfList\n')
	uniqList=[]
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    motifID=each[0]
	    pval=each[1]
	    targetFrac=each[2]
	    tfList=each[3]

	    newList=[]
	    for i in tfList.split(','):
                if conversionDict.has_key(i):
                    tf=conversionDict[i]
		
		    if not tf in uniqList:
			newList.append(tf)
			uniqList.append(tf)

	    if len(newList) > 0:
		newLine=[motifID, pval, targetFrac, ','.join(newList)]
		output1.write('\t'.join(newLine)+'\n')

	input1.close()
	output1.close()



getMotifInfo()
modifyName()




