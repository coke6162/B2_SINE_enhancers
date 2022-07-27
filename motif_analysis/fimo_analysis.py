#!/usr/bin/env python
# encoding: utf-8
"""
fimo_analysis.py

Created by Edward Chuong on 2009-06-08.
Copyright (c) 2009 Baker Lab, Stanford University. All rights reserved.
"""

import sys
import getopt
from ed import *

help_message = '''

-i fimo.txt <output from fimo --text --parse-genomic-coord --thresh 1e-4 motif.meme.txt genome.fa>
-o bed_output <folder of bed files>
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def main(argv=None):

	inFile = None 
	name = None
	outFolder = "beds"
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hi:vo:n:", ["help", "input=", "output="])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-i", "--input"):
				inFile = value
			if option in ("-o", "--output"):
				outFolder = value
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, help_message
		return 2

		
	if outFolder: createFolder(outFolder)
	f = open(inFile)
	inFileString = inFile.split("/")[-1].split(".")[0]
	
	motif2beds = dict()
	for line in f:
		if line.startswith('#') or len(line.strip())==0: continue
		l = line.strip().split()
		
		motif = l[0]
		motif_alt_id = l[1]
		chrm = l[2]
		hitStart = l[3]
		hitEnd = l[4]
		strand = l[5]
		score = l[6]
		#asdfas


		#bedPos = l[1] + ":" + hitStart + "-" + hitEnd
		bedPos = (chrm, hitStart, hitEnd, strand, score)
        
		if motif not in motif2beds:
			motif2beds[motif] = set()
		motif2beds[motif].add(bedPos)
		
		

	f.close()
	
	print "FIMO summary for %s" % inFile
	
	for motif in motif2beds:
		print "%s\t%d" % (motif, len(motif2beds[motif]))
	
	if outFolder:
		print "Outputting bed files into %s" % outFolder
		totalOutFileString = outFolder + "/all.bed"
		t = open(totalOutFileString, 'w')
		for motif in motif2beds:
			if not name:
				name = inFileString
			outFileString = outFolder + "/" + motif + ".bed"
			print "... %s" % outFileString
			o = open(outFileString, 'w')
			for bedPos in motif2beds[motif]:
				#if not bedPos.startswith("chr"): continue
				(chrm, hitStart, hitEnd, strand, score) = bedPos
				
				line = "\t".join([chrm, hitStart, hitEnd, motif, score, strand]) + "\n"
				t.write(line)
				o.write(line)
			o.close()
		print "... %s" % totalOutFileString
		t.close()


if __name__ == "__main__":
	sys.exit(main())
