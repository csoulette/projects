#!/usr/bin/env python3


########################################################################
# File: permuteSampleGroups.py
#  executable: permuteSampleGroups.py
# Purpose: permute feature within and across sample sets.
#
#          
# Author: Cameron M. Soulette
# History:      cms 07/24/2017 Created
#
# This program was written in emacs.
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################
import random
import os, sys
import numpy as np
import matplotlib.pyplot as plt

from multiprocessing import Pool




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
        self.parser = argparse.ArgumentParser(description = '''permuteSampleGroups.py - a tool to permute within sample sets.''',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -i input_file -p N')
        # Add args
        self.parser.add_argument('-i', '--input_file', action = 'store', required=True, help='Input file. [Default : req*]')
        self.parser.add_argument('-p', '--num_threads', action = 'store', required=False, default=2, type=int, help='Num of threads [Default: 2]')
        
        
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))


########################################################################
# Functions
# 
########################################################################


def permute(x):
    
    
    gene, geneDict = x
    
    # Make list of sample IDs
    allSamples = geneDict['donors']
    donorList = list(allSamples.keys())

    # Get set of outlier mutants (observed)
    observedSet = geneDict['observed']
    
    # Calc fraction/absolute outlier number.
    observedFraction = float(len(observedSet)) / float(len(list(allSamples.keys())))
    observedNum = float(len(observedSet)) 

    
    
    nulls = list()
    iterations = 10000
    tests = 0

    for i in range(iterations):
            
        tests += 1

        # Sample random exon z-score.
        choices = [ random.choice(allSamples[donor]) for donor in donorList ] 
        
        #Some choices may be 'nan'. Adjust sample size.
        sampleSize = len([x  for x in choices if not np.isnan(x[1])])
        
        if sampleSize == 0:
            tests -= 1
            continue
        else:
            # Success must be z-score>3 and >=10 delta psi.
            success = float(sum([1 for x in choices 
                                 if abs(x[1])>=float(3) and
                                 x[0] == '>=10'])) / float(sampleSize)

            # add null observation to list
            nulls.append(success)

    # pvalue is equal to the number of times outlier z fraction was equal to or higher than experimental observed.
    pvalue = float(sum([1 for x in nulls if x>=observedFraction]))/float(tests)
    print(gene[0], gene[1], pvalue, observedNum, observedFraction, len(list(allSamples.keys())), len(geneDict['events']), sep="\t")




def readPermTable(inFile):

    '''
    Reads my special formatted table and returns data dict.
    '''
    geneDict = dict()
    
    with open(inFile,'r') as lines:
        for line in lines:
            cols = line.rstrip().split()
            ens, hugo, donorID, histo, psi, dPSI, zScore, event, numMutations, outlier = cols
            zScore = np.float32(zScore)

            geneKey = (ens, hugo)
            if geneKey not in geneDict:
                geneDict[geneKey] = {'donors':dict(), 'observed':set(), 'events':set()}
                
            if donorID not in geneDict[geneKey]['donors']:
                geneDict[geneKey]['donors'][donorID] = list()

            # Label ->   mut   ... non mut ... non mut
            # Value -> [zscor] ... [zscor] ... [zscore]
            
            # at this point I don't need to know the label, just how many mut labels per gene.
            # also, how many outlier (experimental observed)
            geneDict[geneKey]['donors'][donorID].append((dPSI,zScore,numMutations))

            geneDict[geneKey]['events'].add(event)

            if outlier == "outlier":
                geneDict[geneKey]['observed'].add(donorID)

        return geneDict


def main():
    '''
    TBD
    '''
    myCommandLine = CommandLine()
    inFile = myCommandLine.args['input_file']
    threads = myCommandLine.args['num_threads']
    
    p = Pool(threads)

    cmdList = list()

    geneDict = readPermTable(inFile)
    
    
    print("ensemblID","hugoID","raw_pvalue","observedNum","observedFraction","totalMutants","totalExons",sep="\t")


    # Filter for genes with more than 1 outlier AND
    # Filter for genes with more than 10000 unique combinations.

    outlierNumThreshold = 1
    iterationCutoff = 10000
    cmdList = [(gene,gDict) for gene,gDict in geneDict.items()
               if len(gDict['observed'])>outlierNumThreshold and
               (len(list(gDict['donors'].keys()))**len(gDict['events']))> iterationCutoff ]

    # Map it out.
    p.map(permute, cmdList)

    

########################################################################
# Main
# Here is the main program
# 
########################################################################

if __name__ == "__main__":
    main();
