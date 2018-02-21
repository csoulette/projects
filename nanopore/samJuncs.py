#!/usr/bin/env python3


########################################################################
# File: samJuncs.py
#  executable: samJuncs.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 07/26/2017 Created
#
# This program was written in emacs.
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import numpy as np
from multiprocessing import Pool
import re
import pysam

########################################################################
# CommandLine
########################################################################

class CommandLine(object) :
	'''
	Handle the command line, usage and help requests.
	CommandLine uses argparse, now standard in 2.7 and beyond. 
	it implements a standard command line argument parser with various argument options,
	and a standard usage and help,
	attributes:
	myCommandLine.args is a dictionary which includes each of the available command line arguments as
	myCommandLine.args['option'] 
	
	methods:
	
	'''
	
	def __init__(self, inOpts=None) :
		'''
		CommandLine constructor.
		Implements a parser to interpret the command line argv string using argparse.
		'''
		import argparse
		self.parser = argparse.ArgumentParser(description = 'samJuncs.py - lorem ipsium.',
											 epilog = 'Please feel free to forward any questions or concerns to /dev/null', 
											 add_help = True, #default is True 
											 prefix_chars = '-', 
											 usage = '%(prog)s -g gtf_file -f out_format ')
		# Add args
		self.parser.add_argument('-i', '--sam_file', type=str, action = 'store', required=True, help='Input SAM/BAM file.')
		self.parser.add_argument('--isSam', action = 'store_true',  default=False, help='File is sam')
		
				
		if inOpts is None :
			self.args = vars(self.parser.parse_args())
		else :
			self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# Sequence Alignment File
########################################################################

class SAM(object):
	'''
	Handles sequence alignment format file input and output.
	'''

	def __init__(self, inFile=None, isBam=False):
		# Attributes
		self.isBam = isBam
		self.inFile = inFile
		self.juncCounts = dict()

		# Start pysam object either as bam or sam object
		if self.isBam == False:
			readSam = "r"
		else:
			readSam = "rb"

		try:
			self.reader = pysam.AlignmentFile(self.inFile, readSam)
		except:
			#File does not exist.
			print("ERROR: Cannot find file %s. Exiting!" % self.inFile, file=sys.stderr)
			sys.exit(1)
	
		#print(self.reader.find_introns((read for read in self.reader.fetch() if read.is_reverse)))

		#sys.exit(1)

	def countJuncs(self):
		'''
		Reads self.reader and reports junction counts for each junction.
		Returns dictionary with junction counts per junction per chromosome.
		'''

		strandInfo = {0:'+', 16:'-'}
		
		try:
			firstRead = self.reader.fetch('chrMT')
		except:
			pysam.index()
			
		for read in self.reader.fetch():

			try:
				strand = strandInfo[read.flag]
			except:
				continue

			qName = read.query_name
			chromosome = read.reference_name
			
			refPos = read.pos
			refEnd = read.pos
			

			startPos = read.pos
			cigar = read.cigar
	
			for num, flagTuple in enumerate(cigar,1):
				flag, length = flagTuple 
				if flag not in [0,2,3,7,8]:
					continue
					
				if flag == 3:
					if chromosome not in self.juncCounts:
						self.juncCounts[chromosome] = dict()
					if (refEnd, refEnd+length) not in self.juncCounts[chromosome]:
						self.juncCounts[chromosome][(refEnd, refEnd+length)] = int()

					self.juncCounts[chromosome][(refEnd, refEnd+length)] += 1
					if refEnd+length == 74177847:
						print(qName)
				refPos = refEnd+length
				refEnd = refPos

		return self.juncCounts

def main():
	'''
	TDB
	'''
	myCommandLine = CommandLine()
	
	alignmentFile = myCommandLine.args['sam_file']
	isBam = myCommandLine.args['isSam']
	
	sObj = SAM(alignmentFile, isBam)

	d = sObj.countJuncs()
	for c,j in d.items():
		for i in j:
			print(c,"\t".join(str(x) for x in i[:2]), d[c][i], sep="\t")
	
########################################################################
# Main
# Here is the main program
# 
########################################################################

if __name__ == "__main__":
	main();        
