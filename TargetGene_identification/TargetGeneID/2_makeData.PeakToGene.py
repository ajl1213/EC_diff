#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
NearbyWindow=[5000,'5kb']
GTF='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
peakLabelFile='/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/DiffPeaks/H3K27ac.DiffPeakLabel.txt'
SigInterFile='/home/ajl1213/Projects/Endothelial/data/HiC/FitHiC/SigInter/MergedEC.10kb.1e-1.FitHiC.ProCre.txt'
SnpPeakDIR=DIR+'/RESULT'

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
'Colorectal_cancer',
'Ulcerative_colitis',
'Rheumatoid_arthritis'
]

os.system('mkdir '+DIR+'/BedFiles')
os.system('mkdir '+DIR+'/TargetGeneList')


def loadGeneInfo():
    geneDict={}
    input1=open(GTF,'r')
    all_input1=input1.readlines()
    output1=open(DIR+'/BedFiles/GeneTSS.bed','w')
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        chrID=each[2]
        tssPos=each[3]
        if each[6]=='protein_coding':
            geneDict[ensembleID]=geneID
            output1.write(chrID+'\t'+tssPos+'\t'+str(int(tssPos)+1)+'\t'+ensembleID+'\t'+geneID+'\n')
    input1.close()

    return geneDict


def loadGwasPeak():
    gwasPeakDict={}
    for gwasID in gwasList:
	tmpDict={}
        input1=open(SnpPeakDIR+'/'+gwasID+'.SnpToPeak.txt','r')
        all_input1=input1.readlines()
        for line in all_input1[1:]:
            each=line.strip().split('\t')
	    tagSNP=each[0]
	    snpLabel=each[1]
	    peakID=each[2]
	    peakLabel=each[3]
            if tmpDict.has_key(peakID):
                tmp=tmpDict[peakID]
		tmp.append(tagSNP)
                tmpDict[peakID]=tmp
            else:
                tmpDict[peakID]=[tagSNP]

	gwasPeakDict[gwasID]=tmpDict

        input1.close()

    return gwasPeakDict


def loadPeakLabel():
    peakLabelDict={}
    input1=open(peakLabelFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	peakID=each[0]
	peakLabel=each[2]
	peakLabelDict[peakID]=peakLabel
    input1.close()

    return peakLabelDict


def loadSigInter():
    interDict={}
    input1=open(SigInterFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        interLabel=each[0]
        if interLabel=='PC':
            frag1=each[1]
            geneList=each[2]
            frag2=each[3]
            peakList=each[4]

            for ensembleID in geneList.split(';'):
                for peakID in peakList.split(';'):
                    if interDict.has_key(peakID):
                        tmp=interDict[peakID]
                        tmp.append(ensembleID)
                        interDict[peakID]=tmp
                    else:
                        interDict[peakID]=[ensembleID]
    input1.close()

    return interDict


def getTarget():
    geneDict=loadGeneInfo()
    gwasPeakDict=loadGwasPeak()
    peakLabelDict=loadPeakLabel()
    interDict=loadSigInter()

    #
    nearbyDict={}
    for gwasID in gwasList:
        output1=open(DIR+'/BedFiles/'+gwasID+'.'+NearbyWindow[1]+'.bed','w')
        for peakID in gwasPeakDict[gwasID]:
            chrID=peakID.split(':')[0]
            pt1=int(peakID.split(':')[1].split('-')[0])
            pt2=int(peakID.split(':')[1].split('-')[1])
            newLine=[chrID, str(pt1-NearbyWindow[0]), str(pt2+NearbyWindow[0]), peakID]
            output1.write('\t'.join(newLine)+'\n')
        output1.close()

        runLine='intersectBed -a '+DIR+'/BedFiles/'+gwasID+'.'+NearbyWindow[1]+'.bed -b '+DIR+'/BedFiles/GeneTSS.bed -wa -wb > '+DIR+'/BedFiles/'+gwasID+'.TSS_intersect.'+NearbyWindow[1]+'.bed'
        os.system(runLine)

        tmpDict={}
        input1=open(DIR+'/BedFiles/'+gwasID+'.TSS_intersect.'+NearbyWindow[1]+'.bed','r')
        all_input1=input1.readlines()
        for line in all_input1:
            each=line.strip().split('\t')
            peakID=each[3]
            ensembleID=each[7]
            if tmpDict.has_key(peakID):
                tmp=tmpDict[peakID]
		if not ensembleID in tmp:
		    tmp.append(ensembleID)
                tmpDict[peakID]=tmp
            else:
                tmpDict[peakID]=[ensembleID]
        input1.close()
	nearbyDict[gwasID]=tmpDict

    #
    distantDict={}
    for gwasID in gwasList:
        tmpDict={}
        for peakID in gwasPeakDict[gwasID]:
            if interDict.has_key(peakID):
                for ensembleID in interDict[peakID]:

                    if tmpDict.has_key(peakID):
                        tmp=tmpDict[peakID]
			if not ensembleID in tmp:
			    tmp.append(ensembleID)
                        tmpDict[peakID]=tmp
                    else:
                        tmpDict[peakID]=[ensembleID]
	distantDict[gwasID]=tmpDict

    #
    for gwasID in gwasList:
	output1=open(DIR+'/TargetGeneList/'+gwasID+'.MergedTarget.txt','w')
	output1.write('peakID\tpeakLabel\tensembleID\tgeneID\tassoType\n')

        for peakID in gwasPeakDict[gwasID]:

	    peakLabel=peakLabelDict[peakID]

	    nearbyList=[]
	    if nearbyDict[gwasID].has_key(peakID):
		for ensembleID in nearbyDict[gwasID][peakID]:
		    nearbyList.append(ensembleID)

	    distantList=[]
	    if distantDict[gwasID].has_key(peakID):
		for ensembleID in distantDict[gwasID][peakID]:
		    if not ensembleID in nearbyList:
			distantList.append(ensembleID)

	    for ensembleID in nearbyList:
		if geneDict.has_key(ensembleID):
		    geneID=geneDict[ensembleID]
		    assoType='nearby'
	    
		    newLine=[peakID, peakLabel, ensembleID, geneID, assoType]
		    output1.write('\t'.join(newLine)+'\n')

	    for ensembleID in distantList:
		if geneDict.has_key(ensembleID):
		    geneID=geneDict[ensembleID]
		    assoType='longRange'
	    
		    newLine=[peakID, peakLabel, ensembleID, geneID, assoType]
		    output1.write('\t'.join(newLine)+'\n')
		
	output1.close()	    



getTarget()



