#!/usr/bin/env python3


########################################################################
# File: geneExpressionZ_groups.py
#  executable: geneExpressionZ_groups.py
# Purpose: Convert gene expression table to TSV with z-score.
#
#          
# Author: Cameron M. Soulette
# History:      cms 06/01/2017 Created
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
        self.parser = argparse.ArgumentParser(description = 'geneExpressionZ_groups.py - a tool to convert gene expression matrix to z-score values.',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -o file_name -t expression_matrix -s sample_tsv -p N')
        # Add args
        self.parser.add_argument('-o', '--output_file', action = 'store', required=False, default=sys.stdout, help='Output file name. [Default : stdout]')
        self.parser.add_argument('-t', '--table_file', action = 'store', required=True, help='Input expression tableFastQ [Default: req*]')
        self.parser.add_argument('-s', '--sample_groups', action = 'store', required=True, help='Sample group file (TSV) [Default: req*]')
        self.parser.add_argument('-p', '--num_threads', action = 'store', required=False, default=2, type=int, help='Num of threads [Default: 2]')
        
        
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

        


########################################################################
# Functions
# 
########################################################################

def createGroups(gFile):
    '''
    Reads group TSV file and returns dictionary.
    '''

    groups = dict()
    with open(gFile,'r') as lines: 
        for line in lines:
            sample, group = line.rstrip().split()
            
            if group not in groups:
                groups[group] = set()
                
            
            groups[group].add(sample)

    return groups

def readMatrix(mFile, groups):
    '''
    Reads matrix file and returns matrix list.
    '''

    with open(mFile,'r') as lines:
        firstLine = next(lines).rstrip().split()
        headerDict = {head:col for col,head in enumerate(firstLine,0)}
        otherLines = list()
        
        for line in lines:
            otherLines.append((groups, headerDict, line.rstrip().split()))

    return otherLines


def zScore(dataTuple):

    sampleSizeCutoff, missingDataCutoff = 9, 0.7

    groups, matrixHeader, measuredValues = dataTuple
    
    allSamples = set(matrixHeader.keys())

    finalData = list()
    
    for g, sampleSet in groups.items():
        
        sampleSet = sampleSet.intersection(allSamples)

        if len(sampleSet)<sampleSizeCutoff:
            continue

        dataArray = np.asarray( [ measuredValues[matrixHeader[x]] for x in sampleSet ], dtype=np.float32)
        

        numNan = np.count_nonzero(~np.isnan(dataArray))
        percentNan = (len(dataArray)-numNan)/len(dataArray)

        if percentNan>=missingDataCutoff:
            continue
        
        mean = np.nanmean(dataArray)                                                                                                          
        std = np.nanstd(dataArray)   

        if std < 0.01:
            m = [float(0) if not np.isnan(x) else np.nan for x in dataArray]
        else:
            m = [(x - mean) / std if not np.isnan(x)  else np.nan for x in dataArray]  

        finalData.append((g,measuredValues[0],sampleSet, dataArray, m))
    return finalData


def main():
    '''
    Runs gene expression metric conversion from user specified input to user specified output.

    '''
    myCommandLine = CommandLine()
    
    # Get sample groups
    groups = createGroups(myCommandLine.args['sample_groups'])
    cmdList = readMatrix(myCommandLine.args['table_file'], groups)
    
    p = Pool(myCommandLine.args['num_threads'])

    zScores = p.map(zScore, cmdList)
    for i in zScores:
        for j in i:
            group, event, samples, rawValues, zValues = j
            [print(event, group, sample, rawValues[pos], zValues[pos], sep="\t") for pos, sample in enumerate(samples,0)]

########################################################################
# Main
# Here is the main program
# 
########################################################################

if __name__ == "__main__":
    main();
    

