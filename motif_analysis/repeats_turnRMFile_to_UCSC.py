#!/usr/bin/env python
# encoding: utf-8
"""
repeats_turnRMFile_to_UCSC.py

Created by Edward Chuong on 2010-09-21.
Copyright (c) 2010 Baker Lab, Stanford University. All rights reserved.
"""

import sys
import getopt

help_message = ''' '''

'''
   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  687  17.4  0.0  0.0  chr1      3000002 3000156 (194195276) C  L1_Mur2        LINE/L1             (4310) 1567   1413      1
  917  21.4 11.4  4.5  chr1      3000238 3000733 (194194699) C  L1_Mur2        LINE/L1             (4488) 1389    913      1
  215   3.1  0.0  3.0  chr1      3000734 3000766 (194194666) +  (TTTG)n        Simple_repeat            2   33    (0)      2
  845  23.3  7.6 11.4  chr1      3000767 3000792 (194194640) C  L1_Mur2        LINE/L1             (6816)  912    887      1
  621  25.0  6.5  3.7  chr1      3001288 3001583 (194193849) C  Lx9            LINE/L1             (1596) 6048   5742      3
 1320  19.7 33.2  0.7  chr1      3001723 3002005 (194193427) C  RLTR25A        LTR/ERVK               (0) 1028    625      4
  414   0.0  0.0  0.0  chr1      3002006 3002051 (194193381) +  (TA)n          Simple_repeat            1   46    (0)      5
baker-server:repeatmasker_new echuong$ head ../repeatmasker
repeatmasker/     repeatmasker_new/ 
baker-server:repeatmasker_new echuong$ head ../repeatmasker/
.DS_Store        11/              14/              17/              2/               5/               8/               Un/              chr_masked/      
1/               12/              15/              18/              3/               6/               9/               X/               chromOut.tar.gz  
10/              13/              16/              19/              4/               7/               M/               Y/               
baker-server:repeatmasker_new echuong$ head ../repeatmasker/1/chr1
chr1.fa.out         chr1_random.fa.out  
baker-server:repeatmasker_new echuong$ head ../repeatmasker/1/chr1.fa.out 
   SW  perc perc perc  query      position in query         matching        repeat          position in  repeat
score  div. del. ins.  sequence    begin    end   (left)   repeat          class/family     begin  end (left)  ID
0		1	2		3	4			5		6		7		8	9				10					11		12	13		14
  687  17.4  0.0  0.0  chr1      3000002 3000156 (194195276) C  L1_Mur2        LINE/L1             (4310) 1567   1413      1
  917  21.4 11.4  4.5  chr1      3000238 3000733 (194194699) C  L1_Mur2        LINE/L1             (4488) 1389    913      1
  215   3.1  0.0  3.0  chr1      3000734 3000766 (194194666) +  (TTTG)n        Simple_repeat            2   33    (0)      2
  845  23.3  7.6 11.4  chr1      3000767 3000792 (194194640) C  L1_Mur2        LINE/L1             (6816)  912    887      1
  621  25.0  6.5  3.7  chr1      3001288 3001583 (194193849) C  Lx9            LINE/L1             (1596) 6048   5742      3
 1320  19.7 33.2  0.7  chr1      3001723 3002005 (194193427) C  RLTR25A        LTR/ERVK               (0) 1028    625      4
  414   0.0  0.0  0.0  chr1      3002006 3002051 (194193381) +  (TA)n          Simple_repeat            1   46    (0)      5
baker-server:repeatmasker_new echuong$ head ../../
.DS_Store     hg18/         meug1/        meug1.1/      mm9/          monDelph/     repeat25.txt  uniprobe/     xenTro2/      
baker-server:repeatmasker_new echuong$ head ../../mm9/mm9.repeats.txt
#bin	swScore	milliDiv	milliDel	milliIns	genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEnd	repLeft	id
0		1		2			3			4			5			6			7		8			9		10		11			12			13			14		15
607	687	174	0	0	chr1	3000001	3000156	-194195276	-	L1_Mur2	LINE	L1	-4310	1567	1413	1
607	917	214	114	45	chr1	3000237	3000733	-194194699	-	L1_Mur2	LINE	L1	-4488	1389	913	1
607	215	31	0	30	chr1	3000733	3000766	-194194666	+	(TTTG)n	Simple_repeat	Simple_repeat	2	33	0	2
607	845	233	76	114	chr1	3000766	3000792	-194194640	-	L1_Mur2	LINE	L1	-6816	912	887	1
607	621	250	65	37	chr1	3001287	3001583	-194193849	-	Lx9	LINE	L1	-1596	6048	5742	3
607	1320	197	332	7	chr1	3001722	3002005	-194193427	-	RLTR25A	LTR	ERVK	0	1028	625	4
607	414	0	0	0	chr1	3002005	3002051	-194193381	+	(TA)n	Simple_repeat	Simple_repeat	1	46	0	5
607	2352	193	69	64	chr1	3002051	3002615	-194192817	-	RLTR25A	LTR	ERVK	-404	624	1	4
607	386	268	145	138	chr1	3003212	3003261	-194192171	-	L1M2	LINE	L1	-1343	4800	4769	6



'''





class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def main(argv=None):

	import optparse
	parser = optparse.OptionParser(description=help_message)
	required = optparse.OptionGroup(parser, "Required arguments")
	optional = optparse.OptionGroup(parser, "Optional arguments")
	
	required.add_option('-i', '--input', help="input", dest='input', action='store', default=None, metavar='[input]')
	optional.add_option('-n', '--number', help="some number [%default]", dest='someNumber', default=1000,  metavar='[# things]', type='int', action='store')
	optional.add_option('-e', '--ensembl', help="turn to ensembl mode (no chr) [%default]", dest='ensembl', default=False, action='store_true')
	optional.add_option('-d', '--debug', help="turns on debug mode [%default]", dest="debug", default=False, action='store_true')

	parser.add_option_group(required)
	parser.add_option_group(optional)
	(opts, args) = parser.parse_args()
	
	mandatories = ["input"]
	for m in mandatories:
		if not opts.__dict__[m]:
			print "\n*** Mandatory option \"%s\" missing!\n" % m
			print help_message
			sys.exit(1)
	

	f = open(opts.input)
	for line in f:
		if line.startswith('#') or len(line.strip())==0: continue
		l = line.strip().split()
		if l[0] == "SW":
			print "#bin	swScore	milliDiv	milliDel	milliIns	genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEnd	repLeft	id"
			continue
		if l[0] == "score": continue
		if len(l) == 0: continue
		
		bin = "607"
		swScore = l[0]
		milliDiv = str(float(l[1]) * 10)
		milliDel = str(float(l[2]) * 10)
		milliIns = str(float(l[3]) * 10)
		genoName = l[4]
		
		## for ensembl files
		if opts.ensembl:
			genoName = genoName.lstrip("chr")
		
		##
		genoStart = str(int(l[5]) - 1)
		genoEnd = l[6]
		genoLeft = "-" + l[7].strip("()")
		strand = "+"
		if l[8] == "C": strand = "-"
		repName = l[9]
		classAndFam = l[10].split("/")
		repClass = classAndFam[0]
		if len(classAndFam) == 2:
			repFamily = classAndFam[1]
		else:
			repFamily = classAndFam[0]
		if l[11].startswith("("):
			repStart = str(-int(l[11].strip("()")))
		else:
			repStart = l[11]
		repEnd = l[12]
		if l[13].startswith("("):
			repLeft = str(-int(l[13].strip("()")))
		else:
			repLeft = l[13]		
		id = ""
		if len(l) > 14: id = l[14]
		print "\t".join([bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id])
		
	
	

	f.close()
	

if __name__ == "__main__":
	sys.exit(main())
