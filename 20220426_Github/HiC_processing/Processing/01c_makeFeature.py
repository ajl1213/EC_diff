#!/home/ajl1213/anaconda2/bin/python


import sys
import gzip



resolution=sys.argv[1]
if resolution == '5kb':
    winSize=5000
if resolution == '10kb':
    winSize=10000
if resolution == '20kb':
    winSize=20000
if resolution == '40kb':
    winSize=40000

maxDist=int(sys.argv[2])

bin_file=sys.argv[3]
input1=open(bin_file,'r')
all_input1=input1.readlines()
bin_list=[]
for i in all_input1:
	each=i.split()
	bin_list.append(each[0]+'.'+each[1]+'.'+each[2])
input1.close()

coverage_file=sys.argv[4]
input1=open(coverage_file,'r')
all_input1=input1.readlines()
coverage_dic={}
for i in all_input1:
	each=i.split()
	coverage_dic[each[0]+'.'+each[1]+'.'+each[2]]=each[3]
input1.close()
	
count_dic={}

output_file=sys.argv[5]

for line in sys.stdin:
        if 1==1:
                if line[0]!='@':
                        # HWI-ST1113:549:HKWJFADXX:2:1115:20976:94471   163     chr22   16052163        41      64M36S  chrY    16541071        0       TCAAAAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTACAGTATGTGAAGAAGCTTAGTTTTTTCGCTCTTTGCAATAAATCTTGCT    CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJIGIJDGJJJJIJJJJJJIIIIJGHIIJJHHGHHFFCEFFDDDDDDDDDDDDEDDDDDEDDDDDC    NM:i:0  MD:Z:64 AS:i:64 XS:i:51 SA:Z:chrY,16540953,+,62S38M,27,0;
                        [id1,flag1,chr_from1,loc_from1,mapq1,cigar1,chr_from2, loc_from2, dist, read1, read_qual1]=line.split('\t')[0:11]
                        pos1_chr=chr_from1
                        pos2_chr=chr_from2
                        pos1=int(loc_from1)
                        pos2=int(loc_from2)
                        dist=abs(int(dist))

                        id1=(int(int(pos1/winSize)*winSize))
                        id2=(int(int(pos2/winSize)*winSize))
			if id1<=id2:
				id=str(id1)+'.'+str(id2)
			else:
				id=str(id2)+'.'+str(id1)
			if count_dic.has_key(id):
				count_dic[id]+=1
			else:
				count_dic[id]=1


output1=gzip.open(output_file,'w')
output1.write('frag1\tfrag2\tcov_frag1\tcov_frag2\tfreq\tdist\n')

for i in bin_list:
	each_i=i.split('.')
	for j in bin_list:
		each_j=j.split('.')
		if int(each_i[1])<=int(each_j[1]):
			id=each_i[1]+'.'+each_j[1]	
			dist=int(each_j[1])-int(each_i[1])
			if dist <= maxDist:
			    if count_dic.has_key(id):
				    output1.write(i+'\t'+j+'\t'+coverage_dic[i]+'\t'+coverage_dic[j]+'\t'+str(count_dic[id])+'\t'+str(dist)+'\n')
			    else:
				    output1.write(i+'\t'+j+'\t'+coverage_dic[i]+'\t'+coverage_dic[j]+'\t'+str(0)+'\t'+str(dist)+'\n')




