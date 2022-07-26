#!/usr/bin/env python
# encoding: utf-8
"""
repeats_createExpandedRepeatFile.py

Input: UCSC-formatted repeatmasker out file
Output: a bed file of repeats with expanded coordinates based on "full-length" boundaries (used to align Deeptools Heatmaps)

Created by Edward Chuong on 2010-09-09.
"""

import sys
import getopt


import operator
help_message = '''

'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def main(argv=None):

	import optparse
	parser = optparse.OptionParser(description=help_message)
	required = optparse.OptionGroup(parser, "Required arguments")
	optional = optparse.OptionGroup(parser, "Optional arguments")
	
	required.add_option('-i', '--input', help="input repeatmasker file table", dest='input', action='store', default=None, metavar='[input]')

	parser.add_option_group(required)
	parser.add_option_group(optional)
	(opts, args) = parser.parse_args()
	
	mandatories = ["input"]
	for m in mandatories:
		if not opts.__dict__[m]:
			print "\n*** Mandatory option \"%s\" missing!\n" % m
			print help_message
			sys.exit(1)


	# begin program
	f = open(opts.input)
	
	locations = []
	starts = []
	ends = []
	allPoints = []
	
	found = 0
	for line in f:
		if line.startswith('#') or len(line.strip())==0: continue
		l = line.strip().split()
		score = l[1]
		repName = l[10]
		repFamily = l[12]
		repClass = l[11]
		contig = l[5]
		start = int(l[6])
		end = int(l[7])
	
					
		strand = l[9]
		repStart = int(l[13])
		repEnd = int(l[14])
		repLeft = int(l[15])
		if strand == "+":
			realStart = start - repStart
			realEnd = end - repLeft
		else:
			realStart = start + repStart
			realEnd = end + repLeft

	
		print "\t".join([contig, str(realStart), str(realEnd), repName, score, strand])
	
	f.close()
	
if __name__ == "__main__":
	sys.exit(main())
