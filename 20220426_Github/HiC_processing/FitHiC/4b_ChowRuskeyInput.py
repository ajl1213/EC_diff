#!/home/ajl1213/anaconda2/bin/python


import os



DIR=os.getcwd()
inputDIR=DIR+'/SigInter'
res='10kb'
qval='1e-1'

# 000 100 010 110 001 101 011 111 

def makeInput():

    a=[]
    b=[]
    c=[]

    input1=open(inputDIR+'/ME.'+res+'.'+qval+'.FitHiC.all.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
	frag1=each[0]
	frag2=each[1]
	interID=frag1+'-'+frag2
        a.append(interID)
    input1.close()

    input1=open(inputDIR+'/EC12.'+res+'.'+qval+'.FitHiC.all.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
	frag1=each[0]
	frag2=each[1]
	interID=frag1+'-'+frag2
        b.append(interID)
    input1.close()

    input1=open(inputDIR+'/EC48.'+res+'.'+qval+'.FitHiC.all.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
	frag1=each[0]
	frag2=each[1]
	interID=frag1+'-'+frag2
        c.append(interID)
    input1.close()

    a=set(a)
    b=set(b)
    c=set(c)

    ## compute combinations
    # 000
    c000=0
    # 100
    c100=len(a-b-c)
    # 010
    c010=len(b-a-c)
    # 110
    c110=len(a&b-c)
    # 001
    c001=len(c-a-b)
    # 101
    c101=len(a&c-b)
    # 011
    c011=len(b&c-a)
    # 111 
    c111=len(a&b&c)

    ##
    labelList=['c000','c100','c010','c110','c001','c101','c011','c111']
    valList=[c000,c100,c010,c110,c001,c101,c011,c111]

    output1=open(DIR+'/VennInput.txt','w')
    output1.write('comb\tnVal\n')
    for i in range(len(valList)):
        output1.write(str(labelList[i])+'\t'+str(valList[i])+'\n')

    output1.close()



makeInput()






