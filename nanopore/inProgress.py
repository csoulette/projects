#!/usr/bin/env python3


########################################################################
# File: 
#  executable: 
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 02/22/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import numpy as np
from multiprocessing import Pool
import re
from intervaltree import Interval, IntervalTree


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
		self.parser = argparse.ArgumentParser(description = ' TBD - TBD.',
											 epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
											 add_help = True, #default is True 
											 prefix_chars = '-', 
											 usage = '%(prog)s --gtf annotations.gtf -s aligned.sam --isBAM ')
		# Add args
		self.parser.add_argument('--gtf', action = 'store', required=True, help='Input GTF file.')
		self.parser.add_argument('-s', '--sam', action = 'store', required=True, help='Input alignment file.')
		self.parser.add_argument('--isBAM', action = 'store_true', default=False, required=False, help='Input alignment file is BAM formatted. [Default: SAM]')
				
		if inOpts is None :
			self.args = vars(self.parser.parse_args())
		else :
			self.args = vars(self.parser.parse_args(inOpts))



########################################################################
# GTF
########################################################################

class GTF(object):
	'''
	Reads GTF formatted file. Format is specific to gencode.
	'''
	def __init__(self, gtfFile=None): 

		try:
			self.gtfFile = open(gtfFile, 'r')
		except:
			print("File not found.", file=sys.stderr)
			sys.exit(1)

	def getLine(self):
		'''
		File reading generator. Calling on this returns line from self.gftFile
		'''

		with self.gtfFile as lines:
			for line in lines: 
				if line[0] == "#":
					continue
				yield line.rstrip()
		
	def buildTranscripts(self):
		'''
		Reads self.gftFile until EOF. Returns gene object with collection of exons.
		'''

		self.transcripts = dict()
		for line in self.getLine():

			cols = line.split("\t")
			feature = cols[2]

			if feature != "exon":
				continue

			chrom, caller, feature, start, end, phase, strand = cols[:7]
			dataCol = cols[-1]

			# Collect data from exon features.

			transcriptID = (re.search('transcript_id \"(\S+?)\"', dataCol)).group(1) 
			
			if transcriptID not in self.transcripts:
				self.transcripts[transcriptID] = Transcript(transcriptID, strand, chrom)

			start, end = int(start), int(end)
			exonID = (start, end)

			self.transcripts[transcriptID].exons.append(exonID)            

		return self.transcripts

########################################################################
# Gene
########################################################################

class Gene(object):
	'''
	Holds important transcript data.
	Transcipt is called by the GTF class object.
	This class can be used to retrieve various transcript data. 
	'''

	def __init__(self, name = None, strand = None):
		
		self.name = name
		self.strand = strand

		self.transcripts = dict()
		self.exons = set()

	def genomicRange(self):
		'''
		Returns min and max exon coordinates.
		'''
		self.start = min([x[0]] for x in list(self.exons))
		self.end = max([x[1]] for x in list(self.exons))

		return (self.start, self.end)


########################################################################
# Transcript
########################################################################

class Transcript(object):
	'''
	Holds important transcript data.
	Transcipt is called by the GTF class object.
	This class can be used to retrieve various transcript data. 
	'''

	def __init__(self, name = None, strand = None, chrom = None):
		
		self.name = name
		self.strand = strand
		self.chromosome = chrom
		self.exons = list()

	def genomicRange(self):
		'''
		Returns min and max exon coordinates.
		'''
		self.start = min([x[0]] for x in list(self.exons))
		self.end = max([x[1]] for x in list(self.exons))

		return (self.start, self.end)

	def makeIntrons(self):

		self.introns = list()
		self.junctions = list()

		if self.strand == "-":
			self.exons = self.exons[::-1]

		for num, exon in enumerate(self.exons, 0):

			if num == 0:
				continue

			current, null = exon
			null, previous = self.exons[num-1]

			current = current - 1

			intronID = (previous, current)
			self.introns.append(intronID)
			self.junctions.append(current)
			self.junctions.append(previous)

		return self.introns



########################################################################
# Functions
########################################################################



########################################################################
# Main
# Here is the main program
# 
########################################################################
def main():
	'''
	TDB
	'''
	myCommandLine = CommandLine()
	
	gtfFile = myCommandLine.args['gtf']
	alignmentFile = myCommandLine.args['sam']
	isBam = myCommandLine.args['isBAM']

	gtfObj = GTF(gtfFile)
	txns = gtfObj.buildTranscripts()
	
	# Build intron Tree
	intTree = dict()
	wiggle = 10
	allInt = dict()
	for txn, tObj in txns.items():
		print(txn, tObj.genomicRange(), tObj.makeIntrons())
		
		if tObj.chromosome not in allInt:
			allInt[tObj.chromosome] = list()

		allInt[tObj.chromosome].extend(tObj.junctions)


	for ch, ints in allInt.items():
		#print(ch, ints)
		ints = [(pos-wiggle,pos+wiggle,pos) for pos in set(list(ints))]
		intTree[ch] = IntervalTree.from_tuples(ints)
		
	for i in range(1000000):	
		l = intTree["chr1"][12228]


if __name__ == "__main__":
	main()

