#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import sys 


DIR=os.getcwd()
gtf=sys.argv[1]
GeneCountFile=sys.argv[2]
H3K27acCountFile=sys.argv[3]
inputFile=sys.argv[4]
index=str(sys.argv[5])


def loadGeneInfo():
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


def getCor(elementList):
    #   
    rnaCountDict={}
    input1=open(GeneCountFile,'r')
    all_input1=input1.readlines()
    rnaSamples=all_input1[0].strip().split('\t')[2:]
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        vals=numpy.array(each[2:], dtype='float')
        tmpDict={}
        idx=0
        for val in vals:
            tmpDict[rnaSamples[idx]]=val
            idx+=1
        rnaCountDict[ensembleID]=tmpDict
    input1.close()

    #   
    h3k27acCountDict={}
    input1=open(H3K27acCountFile,'r')
    all_input1=input1.readlines()
    h3k27acSamples=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        vals=numpy.array(each[1:], dtype='float')
        tmpDict={}
        idx=0
        for val in vals:
            tmpDict[h3k27acSamples[idx]]=val
            idx+=1
        h3k27acCountDict[peakID]=tmpDict
    input1.close()

    #
    tmpSamples = set(rnaSamples) & set(h3k27acSamples)
    commonSamples = []

    for sampleID in tmpSamples:
        if not sampleID.count('hESC') > 0:
            if not sampleID.count('nonEC') > 0:
                commonSamples.append(sampleID)

    print 'Number of samples.. '+str(len(commonSamples))
    print commonSamples
    print 'getting correlation values..'
    corDict={}
    idx=1
    for i in set(elementList):
        print str(idx)+' / '+str(len(set(elementList)))
        idx+=1

        ensembleID=i.split(';')[0]
        peakID=i.split(';')[1]

        if rnaCountDict.has_key(ensembleID):
            rnaValList=[]
            h3k27acValList=[]

            for sampleID in commonSamples:
                rnaVal= rnaCountDict[ensembleID][sampleID]
                h3k27acVal= h3k27acCountDict[peakID][sampleID]

                rnaValList.append(rnaVal)
                h3k27acValList.append(h3k27acVal)

            rnaZscoreList=[]
            meanRnaVal=numpy.mean(rnaValList)
            stdRnaVal=numpy.std(rnaValList)
            for j in rnaValList:
                zVal=(j-meanRnaVal)/stdRnaVal
                rnaZscoreList.append(zVal)

            h3k27acZscoreList=[]
            meanRnaVal=numpy.mean(h3k27acValList)
            stdRnaVal=numpy.std(h3k27acValList)
            for j in h3k27acValList:
                zVal=(j-meanRnaVal)/stdRnaVal
                h3k27acZscoreList.append(zVal)

            corVal=str(numpy.corrcoef(rnaZscoreList, h3k27acZscoreList)[0,1])
            keyID=ensembleID+';'+peakID
            corDict[keyID]=corVal

    print 'done'

    return corDict


def makeData():
    geneDict=loadGeneInfo()

    #
    elementList=[]
    input1=open(inputFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        creID=each[2]

        assoID=ensembleID+';'+creID
        elementList.append(assoID)
    input1.close()

    corDict=getCor(elementList)

    #
    output1=open(DIR+'/PCC/PCC.'+str(index)+'.txt','w')
    output1.write('ensembleID\tgeneID\tcreID\tPCC\n')
    for i in elementList:
        if corDict.has_key(i):

            ensembleID = i.split(';')[0]
            geneID=geneDict[ensembleID]
            peakID = i.split(';')[1]

            corVal=corDict[i]
            newLine=[ensembleID, geneID, peakID, corVal]
            output1.write('\t'.join(newLine)+'\n')
    output1.close()



makeData()




