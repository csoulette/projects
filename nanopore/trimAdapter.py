#!/usr/bin/env python3


########################################################################
# File: trimAdapter.py
#  executable: trimAdapter.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 02/15/2018 Created
#r
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################



import regex
import os, sys
import re
import numpy as np

from skbio import TabularMSA, DNA
from skbio.alignment import local_pairwise_align_ssw
from skbio.alignment import StripedSmithWaterman

import matplotlib.pyplot as plt
import pandas as pd

from fuzzywuzzy import fuzz
from fuzzywuzzy import process
########################################################
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
		self.parser = argparse.ArgumentParser(description = 'trimAdapter.py - lorem ipsium.',
											 epilog = 'Please feel free to forward any questions or concerns to /dev/null', 
											 add_help = True, #default is True 
											 prefix_chars = '-', 
											 usage = '%(prog)s -i file.fastq ')
		# Add args
		self.parser.add_argument('-i', '--input_file', type=str, action = 'store', required=True, help='Input FQ/FA file.')
		self.parser.add_argument('--isFasta', action = 'store_true',  default=False, help='File is sam (Not working yet)')
		self.parser.add_argument('-a','--adapter_fasta', type=str, action = 'store', required=False, default="ISPCR", help='Custom adapter fasta sequence files. ISPCR used by default')
		self.parser.add_argument('-s','--score_threshold', type=int, action = 'store', required=False, default=20, help='Match score threshold for matching adapter sequence (Default 40)')				
		self.parser.add_argument('--pA', action = 'store_true',  default=False, help='Look for polyA tail in read sequence (Not working yet)')
		self.parser.add_argument('--plot_thresh', action = 'store_true',  default=False, help='Subsample reads and plot score distribution to infer a good adapter match score threshold (Still experimental)')
		

		if inOpts is None :
			self.args = vars(self.parser.parse_args())
		else :
			self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# FASTQ Reader
########################################################################

class FQ(object):
	'''
	Handles FASTQ sequence format file input and output.
	'''

	def __init__(self, inFile=None):
		# Attributes
		self.inFile = inFile

		try:
			self.reader = open(self.inFile,'r')
		except:
			#File does not exist.
			print("ERROR: Cannot find file %s. Exiting!" % self.inFile, file=sys.stderr)
			sys.exit(1)
	
	def read(self):
		'''
		reads fastq file, reports tuple containing read info.
		'''

		with self.reader as fqLines:
			for header1 in fqLines:
				seq = next(fqLines).rstrip()
				header2 = next(fqLines).rstrip()
				qual = next(fqLines).rstrip()

				readID = (header1.split())[0][1:]
				data = (readID,seq,qual)
				yield(data)

	def subSample(self, num = 2000):
		'''
		Samples N number of fastq reads. Take in int as argument.
		Stores sub sampled reads as list in self.subSamples
		'''

		numSamples = num
		self.subSamples = list()

		for fqEntry in self.read():
			read, seq, qual = fqEntry
			self.subSamples.append(fqEntry)

			if len(self.subSamples)>= numSamples:
				break


########################################################################
# Fuzzy matching function
# 
########################################################################

def refAligner(seq):

	query = StripedSmithWaterman(seq,
		gap_extend_penalty=1)
	return query
	

def matchAdapters(query, ref):
	'''
	Preforms alignment using provided adapters.
	Adapters are list of tuples with (label, adapter seq)
	Returns various data regarding match of each adapter
	'''

	a = query(ref)
	#print(a)
	bb = StripedSmithWaterman(ref)
	b = bb("GTACTCTGCGTTGATACCACTGCTT")
	#print(regex.match('(%s){e<=10}' % ref, 'AATGTACTTCGTTCAGTTACGTATTGCT'), "r")
	#print(fuzz.partial_ratio('AATGTACTTCGTTCAGTTACGTATTGCT', ref), "F1")
	#print(fuzz.partial_ratio('AAGCAGTGGTATCAACGCAGAGTAC', ref), "F2")
	#print(fuzz.partial_ratio('GTACTCTGCGTTGATACCACTGCTT', ref), "F3")

	
	return (a.optimal_alignment_score, a.target_begin, a.target_end_optimal)
	
def pltHist(data, color, dataLabel):

	plt.hist(data, color=color, bins=30, range=(0,100), alpha=0.5, label=dataLabel)



def plotMatchScores(seqList, fivePrimeAdapters, threePrimeAdapters):
	'''
	Defines threshold for adapter match.
	'''

	for adapter in fivePrimeAdapters:
		name, adapterSeq = adapter
		query = adapterSeq

		matches = list()
		allSeqScores = list()
		threePrimeScores =  list()
		fivePrimeScores  = list()

		for seq in seqList:
			five_prime_seq = seq[:100]
			rest_of_seq = seq[100:]
			
			a = query(seq)
			allSeqScores.append(a.optimal_alignment_score)

			a = query(rest_of_seq)
			threePrimeScores.append(a.optimal_alignment_score)

			a = query(five_prime_seq)
			fivePrimeScores.append(a.optimal_alignment_score)


		# Plot histograms
		plt.clf()
		pltHist(threePrimeScores, "#1b9e77", "3' end")
		pltHist(allSeqScores, "#d95f02", "entire read")
		pltHist(fivePrimeScores, "#7570b3", "5' end")

		#Plot max null val
		plt.axvline(x=np.percentile(np.array(threePrimeScores),99), color="red")
		plt.text(np.percentile(np.array(threePrimeScores),99)-2, 400, "3' 99th percentile =  %s" % np.percentile(np.array(threePrimeScores),99), rotation = 90)
		

		plt.legend(title="Read segment queried", fontsize = 'small')
		plt.xlim(0,100)
		plt.xlabel("match score")
		plt.savefig("%s_%s_scoreDistributions.pdf" % (os.getcwd(),name), dpi=600)


	for adapter in threePrimeAdapters:
		name, adapterSeq = adapter
		query = adapterSeq

		matches = list()
		allSeqScores = list()
		threePrimeScores =  list()
		fivePrimeScores  = list()

		for seq in seqList:
			three_prime_seq = seq[-100:]
			rest_of_seq = seq[:-100]
			
			a = query(seq)
			allSeqScores.append(a.optimal_alignment_score)

			a = query(three_prime_seq)
			threePrimeScores.append(a.optimal_alignment_score)

			a = query(rest_of_seq)
			fivePrimeScores.append(a.optimal_alignment_score)


		# Plot histograms
		plt.clf()
		pltHist(allSeqScores, "#d95f02", "entire read")
		pltHist(fivePrimeScores, "#7570b3", "5' end")
		pltHist(threePrimeScores, "#1b9e77", "3' end")
		#Plot max null val
		plt.axvline(x=np.percentile(np.array(fivePrimeScores),99), color="red")
		plt.text(np.percentile(np.array(fivePrimeScores),99)-2, 400, "5' 99th percentile =  %s" % np.percentile(np.array(fivePrimeScores),99), rotation = 90)
		

		plt.legend(title="Read segment queried", fontsize = 'small')
		plt.xlim(0,100)
		plt.xlabel("match score")
		plt.savefig("%s_%s_scoreDistributions.pdf" % (os.getcwd(),name), dpi=600)



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
	
	seqFile = myCommandLine.args['input_file']
	isFasta = myCommandLine.args['isFasta']
	adapterFile = myCommandLine.args['adapter_fasta']
	plotThresh = myCommandLine.args['plot_thresh']
	scoreThresh = myCommandLine.args['score_threshold']
	doPolyA = myCommandLine.args['pA']


	if adapterFile == "ISPCR":
		# Use ISPCR primers         
		fivePrimeAdapters = [
							('ont',refAligner('AATGTACTTCGTTCAGTTACGTATTGCT')),
							('ispcrF',refAligner('AAGCAGTGGTATCAACGCAGAGTAC'))
							]

		threePrimeAdapters = [
							('ispcrR',refAligner('GTACTCTGCGTTGATACCACTGCTT')),
							]			#		  GTACTC   GTTGACG


	if isFasta:
		#Add option for FQ to read as fasta. 
		# Probably pass an argument as attribute to set
		# the reader function to either "read as" fastq or fastq
		pass

	else:
		fqObj = FQ(seqFile)


	# If plotThresh True, then run the sub sample, plot results and exit the program.
	if plotThresh:
		fqObj.subSample()
		plotMatchScores([x[1] for x in fqObj.subSamples], fivePrimeAdapters, threePrimeAdapters)
		sys.exit(1)

	# Run the adapter clipping.
	for num, fqEntry in enumerate(fqObj.read(),1):
		read, seq, qual = fqEntry
		#Match 5prime adapters
		fiveSeq = seq[:100]

		vals = [matchAdapters(x[1], fiveSeq) for x in fivePrimeAdapters]
		fivePrimeClip = max(np.asarray([vals[0][1], vals[1][1]], dtype=np.float32))
		
		#Match 3prime adapters
		threeSeq = seq[-50:]
		vals.extend([matchAdapters(x[1], threeSeq) for x in threePrimeAdapters])
		seqLen = len(seq)
		threePrimeClip = np.float32(vals[-1][0])
		#print(vals)
		print(read, "\t".join(",".join(str(x) for x in j) for j in vals), sep="\t")

		#break
		#fivePrimeClip = int(np.nan_to_num(fivePrimeClip))
		#threePrimeClip = int(np.nan_to_num(threePrimeClip))
		#print("@%s" % read,seq[fivePrimeClip:seqLen-(50-threePrimeClip)],"+" ,qual[fivePrimeClip:seqLen-(50-threePrimeClip)], sep="\n")
	# Additional analysis to look at polyA tails.
	# To be written...
	if doPolyA:
		pass


if __name__ == "__main__":
	main()
