#!/usr/bin/python

import fileinput
import sys


filterDist=5000

for line in sys.stdin:
	if 1==1:
		if line[0]!='@':
                	# HWI-ST1113:549:HKWJFADXX:2:1115:20976:94471	163	chr22	16052163	41	64M36S	chrY	16541071	0	TCAAAAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTACAGTATGTGAAGAAGCTTAGTTTTTTCGCTCTTTGCAATAAATCTTGCT	CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJIGIJDGJJJJIJJJJJJIIIIJGHIIJJHHGHHFFCEFFDDDDDDDDDDDDEDDDDDEDDDDDC	NM:i:0	MD:Z:64	AS:i:64	XS:i:51	SA:Z:chrY,16540953,+,62S38M,27,0;
			[id1,flag1,chr_from1,loc_from1,mapq1,cigar1,chr_from2, loc_from2, dist, read1, read_qual1]=line.split('\t')[0:11]
			pos1_chr=chr_from1
			pos2_chr=chr_from2
			pos1=int(loc_from1)
			pos2=int(loc_from2)
			dist=abs(int(dist))
			
			if pos2_chr=='=':
				if abs(pos1-pos2)>filterDist:
					print line.rstrip("\n")
