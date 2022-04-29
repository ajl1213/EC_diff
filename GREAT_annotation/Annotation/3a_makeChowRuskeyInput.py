#!/home/ajl1213/anaconda2/bin/python


import os



DIR=os.getcwd()
inputDIR=DIR+'/geneCategory'
dataTypeList=['EC_cRE', 'nonEC_cRE']


## a,b,c,d,
# 0000 1000 0100 1100 0010 1010 0110 1110 0001 1001 0101 1101 0011 1011 0111 1111

os.system('mkdir '+DIR+'/ChowInput')



def makeInput():
    for dataType in dataTypeList:
	a=[]
	b=[]
	c=[]
	d=[]
	input1=open(inputDIR+'/'+dataType+'.GreatAnno.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    geneID=each[0]
	    cellInfo=each[1]

	    if cellInfo.count('Early'):
		a.append(geneID)

	    if cellInfo.count('Mid'):
		b.append(geneID)

	    if cellInfo.count('Late'):
		c.append(geneID)

	    if cellInfo.count('Full'):
		d.append(geneID)

	input1.close()

	a=set(a)
	b=set(b)
	c=set(c)
	d=set(d)


	## compute combinations
	# 0000
	c0000=0
	# 1000
	c1000=len(a-b-c-d)
	# 0100
	c0100=len(b-a-c-d)
	# 1100
	c1100=len(a&b-c-d)
	# 0010
	c0010=len(c-a-b-d)
	# 1010
	c1010=len(a&c-b-d)
	# 0110
	c0110=len(b&c-a-d)
	# 1110
	c1110=len(a&b&c-d)
	# 0001
	c0001=len(d-a-b-c)
	# 1001
	c1001=len(a&d-b-c)
	# 0101
	c0101=len(b&d-a-c)
	# 1101
	c1101=len(a&b&d-c)
	# 0011
	c0011=len(c&d-a-b)
	# 1011
	c1011=len(a&c&d-b)
	# 0111 
	c0111=len(b&c&d-a)
	# 1111
	c1111=len(a&b&c&d)

	##
	labelList=['c0000', 'c1000', 'c0100', 'c1100', 'c0010', 'c1010', 'c0110', 'c1110', 'c0001', 'c1001', 'c0101', 'c1101', 'c0011', 'c1011', 'c0111', 'c1111']
	valList=[c0000,c1000,c0100,c1100,c0010,c1010,c0110,c1110,c0001,c1001,c0101,c1101,c0011,c1011,c0111,c1111]


	output1=open(DIR+'/ChowInput/'+dataType+'.VennInput.txt','w')
	output1.write('comb\tnVal\n')
	for i in range(len(valList)):
	    output1.write(str(labelList[i])+'\t'+str(valList[i])+'\n')

	output1.close()



makeInput()



