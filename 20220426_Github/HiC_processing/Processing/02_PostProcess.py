#!/home/ajl1213/anaconda2/bin/python


import os 
import gzip



DIR=os.getcwd()
sampleList=[
'ME_rep1.HiC',
'ME_rep2.HiC',
'ME',
'EC12_rep1.HiC',
'EC12_rep2.HiC',
'EC12',
'EC48_rep1.HiC',
'EC48_rep2.HiC',
'EC48',
'MergedEC'
]

chrList=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',  'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16','chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
#resolution='5kb'
resolution='10kb'

qvalCut=[0.05,'1e-1']
#qvalCut=[0.01,'1e-2']

os.system('mkdir '+DIR+'/OUTPUT')


def mergeInter():
    for sampleID in sampleList:
        print sampleID
        output1=open(DIR+'/OUTPUT/'+sampleID+'.'+resolution+'.'+str(qvalCut[1])+'.SigInter','w')
        output1.write('frag1Chr\tfrag1Mid\tfrag2Chr\tfrag2Mid\tcontactCount\tp_value\tq_value\n')
        i=1
        while i<=23:
            if i==23:
                chrID='chrX'
            else:
                chrID='chr'+str(i)

            print chrID
            input1=gzip.open(DIR+'/DistNorm/'+sampleID+'.'+chrID+'.'+resolution+'.DistNormFeat.gz','r')
            all_input1=input1.readlines()
            for j in all_input1[1:]:
                each=j.split()


		frag1Chr=each[0].split('.')[0]
		frag1Mid=str((int(each[0].split('.')[2]) + int(each[0].split('.')[1]))/2)
		frag2Chr=each[1].split('.')[0]
		frag2Mid=str((int(each[1].split('.')[2]) + int(each[1].split('.')[1]))/2)
		contactCount=each[4]
		p_value=each[11]
		q_value=each[12]

                if float(q_value) < qvalCut[0]:
		    newLine=[frag1Chr, frag1Mid, frag2Chr, frag2Mid, contactCount, p_value, q_value]
		    output1.write('\t'.join(newLine)+'\n')
            i+=1
        output1.close()




mergeInter()



