#!/home/ajl1213/anaconda2/bin/python



import os
import numpy


DIR=os.getcwd()
inputDIR=DIR+'/DiffPeaks'
outFile=DIR+'/DiffPeaks/cRE.STable.txt'
cellTypeList=[
'hESC',
'Mesoderm',
'EarlyEC',
'MidEC',
'LateEC',
'FullEC',
'EarlyNonEC',
'MidNonEC',
'LateNonEC',
'FullNonEC'
]


def makeSTable():
    # ordered gene list
    peakList=[]
    input1=open(DIR+'/DiffPeaks/H3K27ac.DiffPeakLabel.zscore.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        peakList.append(peakID)
    input1.close()

    #   
    pvalDict={}
    qvalDict={}
    for cellType in cellTypeList:
        tmpDict1={}
        tmpDict2={}
        input1=open(inputDIR+'/'+cellType+'.valTable.txt','r')
        all_input1=input1.readlines()
        for line in all_input1[1:]:
            each=line.strip().split('\t')
            peakID=each[0]
            logFC=each[1]
            logCPM=each[2]
            F=each[3]
            pval=float(each[4])
            qval=float(each[5])

            tmpDict1[peakID]=pval
            tmpDict2[peakID]=qval

        pvalDict[cellType]=tmpDict1
        qvalDict[cellType]=tmpDict2

        input1.close()

    #
    finalValDict={}
    for peakID in peakList:
        tmpList1=[]
        for cellType in cellTypeList:
            if pvalDict[cellType].has_key(peakID):
                pval=pvalDict[cellType][peakID]
                tmpList1.append(pval)

        tmpList2=[]
        for cellType in cellTypeList:
            if qvalDict[cellType].has_key(peakID):
                qval=qvalDict[cellType][peakID]
                tmpList2.append(qval)

        minPval=min(tmpList1)
        minQval=min(tmpList2)

        finalValDict[peakID]=[str(minPval), str(minQval)]

    #
    allDict={}
    input1=open(DIR+'/DiffPeaks/H3K27ac.DiffPeakLabel.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        maxStage=each[1]
        maxDataType=each[2]
        allDict[peakID]=[maxStage, maxDataType]
    input1.close()

    #
    output1=open(outFile,'w')
    output1.write('peakID\tmaxTimePoint\tDiffStage\tPval\tQval\n')
    for peakID in peakList:

        maxStage, maxDataType = allDict[peakID]
        pval, qval = finalValDict[peakID]

        newLine=[peakID, maxStage, maxDataType, pval, qval]
        output1.write('\t'.join(newLine)+'\n')

    output1.close()



makeSTable()




