#!/usr/bin/env python3


########################################################################
# File: gtf_parser.py
#  executable: gtf_parser.py
# Purpose: Preform various formatting functions on gtf file. Built around gencode gtf format.
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
        self.parser = argparse.ArgumentParser(description = 'gtf_parser.py - a tool for parsing gtf files. Built around gencode formatted GTF.',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -g gtf_file -f out_format ')
        # Add args
        self.parser.add_argument('-g', '--gtf_file', action = 'store', required=True, help='Input GTF file. [Default : stdin]')
        self.parser.add_argument('-f', '--format_type', action = 'store', required=True, help='Output format type.Input expression tableFastQ [Default: req*]')
        self.parser.add_argument('-0', '--0_based_coords', action = 'store_true', required=False, default=True, help='Output format coords are 0-based [Default: TRUE]')
                
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


    def __init__(self, gtfFile=None, analysis="exon ranges", zeroBased=True):
        

        self.zeroBased = zeroBased

        try:
            self.gtfFile = open(gtfFile, 'r')

        except:
            print("File not found.", file=sys.stderr)
            sys.exit(1)

        if analysis == "exon ranges":
            self.runAnalysis = self.eventRanges

    def eventRanges(self):

        transcriptDict = dict()
        for line in self.readGTF():

            cols = line.split("\t")
            chromo, caller, feature, start, end, phase, strand = cols[:7]
            dataCol = cols[-1]

            if self.zeroBased:
                start = str(int(start) - 1)

            if feature == 'transcript':

                # Set up data for transcript
                
                transcriptID = (re.search('(ENST[^"]+)', dataCol)).group(1) 
                geneID = (re.search('(ENSG[^"]+)', dataCol)).group(1) 
                hugoID =(re.search('gene_name \"([^"]+)', dataCol)).group(1) 
                
                transcriptDict[transcriptID] = {'transcriptData':tuple(), 'exonData':dict()}
                transcriptDict[transcriptID]['transcript data'] = (geneID, hugoID, chromo, start, end, strand)
                
                
            elif feature == 'exon':
                transcriptID = (re.search('(ENST[^"]+)', dataCol)).group(1) 
                exonID = (re.search('(ENSE[^"]+)', dataCol)).group(1) 
                exonNum = (re.search('exon_number ([0-9]+)', dataCol)).group(1) 
                
                transcriptDict[transcriptID]['exonData'][start] = (exonID, exonNum,  end)
                


        self.finalTranscripts = dict()
        for transcript, dataDict in transcriptDict.items():
            tData, eData = dataDict['transcript data'], dataDict['exonData']
            gID, hID, chromo, start, end, strand = tData
            
            self.finalTranscripts[transcript] = Transcript(transcript, gID, hID, chromo, start, end, strand, eData)
        
        return self.finalTranscripts

    def readGTF(self):

        with self.gtfFile as lines:
            for line in lines:
                
                if line[0] == "#":
                    continue

                yield line.rstrip()
                    

########################################################################
# Transcript
########################################################################

class Transcript(object):
    '''
    Holds important transcript data.
    Transcipt is called by the GTF class object.

    This class can be used to retrieve various transcript data. 
    '''

    def __init__(self, tensembl=None, gensembl=None, hugo=None, chromo=None, start=None, end=None, strand=None, exons=None):
        

        self.tensembl = tensembl
        self.gensbml = gensembl
        self.hugo = hugo
        
        self.chromosome = chromo
        self.strand = strand
        self.start = start
        self.end = end
        
        
        
        self.exons = self.defineExons(exons)

        
    def defineExons(self, exons):
        
        '''
        Takes in exon data and returns exon list.
        Exon data is in dict format. Key are exon
        start coordinates, values are:
        (exon num, exon name, end) tuple.
        '''

        exonDict = dict()
        exonItems = list(exons.keys())

        for num, start in enumerate(exons,0):

            exonNum, exonName, end = exons[start]
            start, end = int(start), int(end)

            exonDict[start] = Exon(exonName, exonNum, start, end)

            
        exonItems = list(sorted(exonDict.keys()))
        finalExons = dict()

        if len(exonItems) == 1:
            return None


        for num, start in enumerate(exonItems, 0):
            
            
            exonObj = exonDict[start]
            
            if num == 0:
                nextItem = exonItems[num+1]
                
                exonObj.addPrevExon(None)
                exonObj.addNextExon(exonDict[nextItem])
                
            elif num == len(exonItems)-1:
                prevItem = exonItems[num-1]

                exonObj.addPrevExon(exonDict[prevItem])
                exonObj.addNextExon(None)

            else:
                nextItem = exonItems[num+1]
                prevItem = exonItems[num-1]
                
                exonObj.addPrevExon(exonDict[prevItem])
                exonObj.addNextExon(exonDict[nextItem])
            
            finalExons[start] = exonObj
        return finalExons
########################################################################
# Exon
########################################################################


class Exon(object):
    '''
    Holds important exon data.
    Exon is called by the transcipt class object.
    '''

    def __init__(self, name=None, num=None, start=None, end=None):
        self.eensembl = name
        self.number = num

        self.start = start
        self.end = end

        self.previousExon = None
        self.nextExon = None

    def addNextExon(self, exonObj):

        self.nextExon = exonObj

    def addPrevExon(self, exonObj):
        
        self.previousExon = exonObj


########################################################################
# Functions
########################################################################


def getExonRanges(gtfObj):

    # Do whatever you want here.
    transcripts = gtfObj.runAnalysis()
    
    for transcript, tObj in transcripts.items():
        
        if tObj.exons == None:
            continue

        sortedKeys = sorted(list(tObj.exons.keys()))

        for exonStart in sortedKeys:
            eObj = tObj.exons[exonStart]
            

            if eObj.previousExon != None and eObj.nextExon != None:
                gene, chromo, strand, start, end = tObj.hugo, tObj.chromosome, tObj.strand, eObj.previousExon.start, eObj.nextExon.end
                e1c1, e1c2, e2c1, e2c2, e3c1, e3c = start, eObj.previousExon.end, eObj.start, eObj.end, eObj.nextExon.start, end
                exonCoordString = ":".join(str(x) for x in [e1c1, e1c2, e2c1, e2c2, e3c1, e3c])

                print(chromo, start, end, exonCoordString, gene, strand, sep="\t")
                      


def main():
    '''
    TDB
    '''
    myCommandLine = CommandLine()
    
    inFile = myCommandLine.args['gtf_file']
    analysis = myCommandLine.args['format_type']
    zeroBased = myCommandLine.args['0_based_coords']

    gtfObj = GTF(inFile, analysis, zeroBased)
    
    
    if analysis == "exon ranges":
        getExonRanges(gtfObj)
    

########################################################################
# Main
# Here is the main program
# 
########################################################################

if __name__ == "__main__":
    main();
    
