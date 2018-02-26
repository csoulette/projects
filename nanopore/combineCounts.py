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
import networkx as nx
import matplotlib.pyplot as plt

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
		self.parser.add_argument('-i', "--input_bed", action = 'store', required=True, help='Input bed file.')
		self.parser.add_argument('-j', "--known_junctions", action = 'store', required=True, help='Input bed file.')
				
		if inOpts is None :
			self.args = vars(self.parser.parse_args())
		else :
			self.args = vars(self.parser.parse_args(inOpts))


########################################################################
# Gene
########################################################################

class Junction(object):
	'''
	Handles junction related data
	'''

	def __init__(self, name = None, strand = None):
		
		self.name = name
		self.strand = strand

		self.preceeding = set()
		self.succeeding = set()

		self.edges = set()

class Edge(object):
	def __init__(self, name = None, score = 0):
		
		self.name = name
		self.weight = score

########################################################################
# Functions
########################################################################

def reader(f):
	'''
	Bed file reader.
	'''

	with open(f,'r') as bedEntries:
		for entry in bedEntries:
			entry = entry.rstrip()
			columns = entry.split()
			columns[1] = int(columns[1])
			columns[2] = int(columns[2])
			columns[3] = int(columns[3])
			yield columns

def graphJunctions(cols, jDict, real):

	chrom, c1, c2, score = cols

	if chrom not in jDict:
		jDict[chrom] = dict()

	if c1 not in jDict[chrom]:
		jDict[chrom][c1] = Junction(c1)

	if c2 not in jDict[chrom]:
		jDict[chrom][c2] = Junction(c2)

	node1, node2 = jDict[chrom][c1], jDict[chrom][c2]

	node1.succeeding.add(node2)
	node2.preceeding.add(node1)

	edge = Edge((node1, node2))
	edge.weight = score

	if (c1, c2) in real:

		node1.edges.add(edge)
		#node2.edges.add(edge)

	return jDict

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
	
	inFile = myCommandLine.args['input_bed']
	knownSites = myCommandLine.args['known_junctions']

	real = set()
	with open(knownSites, 'r') as lines:
		for line in lines:
			chrom, c1, c2, strand = line.rstrip().split()
			real.add((int(c1),int(c2)))

	junctionDict = dict()
	for line in reader(inFile):
		junctionDict = graphJunctions(line, junctionDict, real)

	g = dict()
	for chrom, js in junctionDict.items():
		g[chrom] = nx.DiGraph()
		for pos, j in js.items():
				if len(j.edges)<1:
					continue
				outs = ",".join([str(x.name[1].name) for x in j.edges])

				weights = ",".join([str(x.weight) for x in j.edges])
				print(chrom, j.name, outs, weights, sep="\t")


if __name__ == "__main__":
	main()

