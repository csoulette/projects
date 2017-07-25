#!/usr/bin/env python3


########################################################################
# File: preProcessPermute.py
#  executable: preProcessPermute.py
# Purpose: Formate data for permutation test.
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
import os, sys
import numpy as np
from multiprocessing import Pool

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns

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
        self.parser = argparse.ArgumentParser(description = 'preProcessPermute.py - produces a formatted TSV for permuteSampleGroups.py.',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s ')
        # Add args
        self.parser.add_argument('-v', '--variant_file', action = 'store', required=False, default=sys.stdout, help='Output file name. [Default : stdout]')


        self.parser.add_argument('-lb', '--lower_bound', action = 'store', required=False, default=-50, type=int, help='Lower bound threshold. Defines the somatic mutation region of interest. (Inclusive) [Default: -50]')
        self.parser.add_argument('-ub', '--upper_bound', action = 'store', required=False, default=5, type=int, help='Upper bound threshold. Defines the somatic mutation region of interest. [Default: 5]')
        
        
        self.parser.add_argument('-p', '--num_threads', action = 'store', required=False, default=2, type=int, help='Num of threads [Default: 2]')
        
        
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

        
class Gene(object):
    '''
    Handles gene patient/mutation data.

    attributes:
    name --------------> ensembl gene identifier
    hugo --------------> hugo gene identifier

    dPSICutoff --------> delta psi outlier threshold
    zCutoff -----------> z-score outlier threshold

    mutantPatientData -> list of data tuples for each sample with mutation
    allMutants --------> sample ID set of mutants
    outlierMutants ----> sample ID set of outliers
    mutantHistos ------> histology/groups with mutants
    mutantEventDict ---> keeps track of mutated event

    methods:
    addMutant --------> adds sample somatic mutation to aggregate list.
    defOutlier -------> iterates through mutants and identifies outliers.
    addSpladderEvent -> adds gene spladder event to aggregate list.
    calcEventZ -------> goes through spladder events and caluclates z-scores
    zScoreIt ---------> special z-score function to handle np.nan data

    '''



    def __init__(self, geneID = None, hugo = None):
        self.name = geneID
        self.hugo = hugo

        # Cutoffs. 
        self.dPSIcutoff = 0.1
        self.zCutoff = float(3)


        self.allMutants = set()
        self.mutantHistos = set()
        self.outlierMutants = set()
        self.mutantEventDict = dict()
        self.mutantPatientData = list()
        


        self.spladderEvents = list()

        
        
    def addMutant(self, data = None):
        '''
        Appends sample mutation data to aggregate list.
        '''

        self.mutantPatientData.append(data)
       
        
    def defineOutliers(self):
        '''
        Separates outlier mutants and identifies affected histology types.
        '''

        
        for data in  self.mutantPatientData:
            donorID, histology, exonCoords, z, dPSI, snpPos, feature = data

            self.mutantHistos.add(histology)
            self.allMutants.add(donorID)
            
            # Is mutant outlier?
            if dPSI == ">=10" and abs(z)>=self.zCutoff:
                self.outlierMutants.add(donorID)
            
            # Keep track of events with nearby mutation.
            if donorID not in self.mutantEventDict:
                self.mutantEventDict[donorID] = list()

            self.mutantEventDict[donorID].append(exonCoords)
        
            

        
    def addSpladderEvent(self, event, key):
        '''
        Appends spladder event to aggregate list.
        '''

        self.spladderEvents.append((key, event))
                
        
    def calcEventZ(self, spladderHeader, histoTable):


        self.exons = dict()

        # Get the entire set of patients in the spladder table.
        spladderPatients = set(spladderHeader.keys())

        for histo in self.mutantHistos:
            
            patientHistoSet = histoTable[histo]
            spladderHistoSet = patientHistoSet.intersection(spladderPatients)
            mutantPatients = self.allMutants.intersection(spladderHistoSet)

            for eventTuple in self.spladderEvents:
                key, event = eventTuple 
                exonName = key

                if exonName not in self.exons:
                    self.exons[exonName] = dict()
                
                
                psiList = [event[spladderHeader[donorID]] for donorID in spladderHistoSet]
                psiArray = np.asarray(psiList, dtype=np.float32)
                                                                                                                                             
                if len(psiArray)<9:
                    continue
        
                percentNan = (len(psiArray)-np.isnan(psiArray).sum())/len(psiArray)                                                            
             

                if percentNan<0.7:
                    allZ = [np.nan for x in psiArray]
                    mean = np.nan 
                    
                else:        
                    allZ,mean = self.zScoreIt(psiArray)


                
                for mutPatient in mutantPatients:
    
                    pIndex = list(spladderHistoSet).index(mutPatient)
                    z = allZ[pIndex]
                    psi = psiArray[pIndex]
                    
                    if abs(psi-mean)>=self.dPSIcutoff:
                        dPSI = ">=10"
                    else:
                        dPSI = "<10"
                    
                    outlierEvents = self.mutantEventDict.get(mutPatient, [])
                    #eventString = ",".join(":".join(x) for x in outlierEvents))
                    numOfMutations = len(outlierEvents)


                    if mutPatient in self.outlierMutants:
                        dataTuple = (mutPatient, histo, psi, dPSI, z,":".join(exonName),numOfMutations, "outlier")
                    else:
                        dataTuple = (mutPatient, histo, psi, dPSI, z,":".join(exonName),numOfMutations, "non-outlier")
                    if mutPatient not in self.exons[exonName]:
                        self.exons[exonName][mutPatient] = dataTuple

                    elif abs(dataTuple[4])>abs(self.exons[exonName][mutPatient][4]):
                        self.exons[exonName][mutPatient] = dataTuple
 
    def zScoreIt(self, psiArray):

        '''
        Calculate a z-score.
        '''
        mean = np.nanmean(psiArray)                                                                                                    
        std = np.nanstd(psiArray)                                                                                                            

        if std < 0.01:
            m = [float(0) if not np.isnan(x) else np.nan for x in psiArray]
        else:
            m = [(x - mean) / std if not np.isnan(x)  else np.nan for x in psiArray]  
        return m, mean


########################################################################
# Functions
# 
########################################################################


def qcPlot(array, color, pltType, title, xlab, ylab, saveAs):
    '''
    Various QC plotting methods.
    '''

    if pltType == "count":
        plt.clf()
        plt.style.use('classic')
        f, ax = plt.subplots(figsize=(6, 3))
        sns.countplot(array, color=color)
        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title, size=8)
        plt.autoscale()
        plt.xticks(rotation='vertical', size=6)
        f.savefig(saveAs, transparent=True, bbox_inches='tight', pad_inches=0.25, dpi=600)
        plt.clf()
    
    if pltType == "hist":
        plt.clf()
        plt.style.use('classic')
        f, ax = plt.subplots(figsize=(6, 3))
        plt.hist(array, 50, color=color)
        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title, size=8)
        plt.autoscale()
        
        f.savefig(saveAs, transparent=True, bbox_inches='tight', pad_inches=0.25, dpi=600)
        plt.clf()

def getExonicCoords(refDict, cols, splicingType):
    '''
    Spladder Coords sometimes overlap the same exonic region.
    This function uses a reference table to link overlapping spladder events to a unique exonic coord.
    '''
    if splicingType == "es":
        spladderCoords = ":".join(cols[5:11])
            
    else:
        spladderCoords = ":".join(cols[5:13])
    chrom = cols[2]
    exonicCoords = refDict[(chrom,spladderCoords)]
    return exonicCoords


def filterMutations(mutationFile, exonReferenceCoords, lowerBound, upperBound):
    '''
    Filters pre computed somatic mutation data and returns list of mutations of interest.

    input file format details can be found here:
    /pod/pstore/groups/brookslab/csoulette/projects/splicing_mutations/somatic_mutation_file_format.txt
    '''
    geneObjects = dict()
    snpPositions = list()
    with open(mutationFile,'r') as lines:
        for line in lines:
            cols = line.rstrip().split()
            
            # Gene info.
            hugo, gene, chromosome = cols[0], cols[1], cols[2]
            
            # Patient/Mutation info.
            histology, donorID = cols[15], cols[14]
            snpPos, feature = int(cols[18]), cols[19]
            
            # Associated splicing event info.
            splicingType, psi, dPSI, z  = cols[4], cols[23], cols[22], cols[21]


            #### Filters ######
            
            # skip mutations outside of my window.
            if snpPos<lowerBound or snpPos>upperBound:
                continue
            
            if z == 'nan' or psi == 'nan':
                continue

            ### end filters ###


            # Get unique coords for overlapping spladder events
            exonCoords = getExonicCoords(exonReferenceCoords, cols, splicingType) 
            
            # Create gene object if not exists.
            if gene not in geneObjects:
                geneObjects[gene] = Gene(gene,hugo)
            
            # Add sample data as list of tuples.
            currentGene = geneObjects[gene]
            snpPositions.append(snpPos)
            dataTuple = (donorID, histology, exonCoords, float(z), dPSI, snpPos, feature)
            currentGene.addMutant(dataTuple)
        
        #QC plot
        title = "mutation relative position hist"
        xlab, ylab = "intron posiiton", "count"
        saveAs = 'snpPos_count.pdf'
        qcPlot(np.asarray(snpPositions,dtype=int), "gray", "count", title, xlab, ylab, saveAs)

        return geneObjects
    
def associateUniqueCoords(exonOverlapFile):
    '''
    Links unique exon coordinate to every spladder event, using spladder event coordinates.

    Substrate file was created by bedtools merging alt 5', 3' and skipped exon spladder coordinates.

    Skipped exon coordinate were filtered for annotated exons. 5'/3' exons were not.

    This function returns a dictionary where the keys are spladder coordinats, and the value is a unique exon coordinate.
    '''


    finalEvents = dict()
    eventCounts = list()
    with open(exonOverlapFile, 'r') as lines:
        for line in lines:
            line = line.rstrip()
            cols = line.split()
            
            chrom, uniqueC1, uniqueC2  = cols[0], cols[1], cols[2]
            val = (chrom, uniqueC1, uniqueC2)
        
            # This column contains a list of spladder events.
            events = cols[3].split(",")

            for event in events:
                finalEvents[(chrom,event)] = val
            
            # Tally the number of events per unique coord.
            eventCounts.append(len(events))
        
    # Plot the number of spladder events for each unique coord.
    # This is just to check any irregulaities in the data.
    eventCounts = np.asarray(eventCounts, dtype=np.int)
    title = "Total Exons:%s, Total Events:%s, Avg exons / event:%s" % (len(eventCounts), np.sum(eventCounts), np.mean(eventCounts))
    ylab, xlab = "number of exons", "number of spladder events"
    saveAs = 'spladderEvent_per_exon_hist.pdf'
    qcPlot(eventCounts, "gray", "hist", title, xlab, ylab, saveAs)

    return finalEvents


def defOutliers(geneObjects):

    '''
    Defines outliers using geneObjects function. 
    Returns updated gene objects and QC plots.
    '''

    [geneObj.defineOutliers() for gene, geneObj in geneObjects.items()]

    outlierNums = np.asarray([len(geneObj.outlierMutants) for gene, geneObj in geneObjects.items()],
                             dtype=np.int)
    

    title = "Total genes:%s" % (str(len(list(geneObjects.keys()))))
    ylab, xlab = "number of genes", "number of mutant outliers"
    saveAs = 'outlierMutants_per_gene_hist.pdf'
    qcPlot(outlierNums, "gray", "count", title, xlab, ylab, saveAs)
                                
    return geneObjects


def addSpladderEvents(sFile, exonReferenceCoords, geneObjects):
    '''
    Adds spladder events to gene objects.
    Returns updated objects.
    '''

    eventsPerGene = list()
    with open(sFile,'r') as lines:
        headerDict = dict()
        headerDict = {head:pos for pos,head in enumerate((next(lines).rstrip()).split(),0)}
    
        for line in lines:
            cols = line.rstrip().split()
            gene, chro, coordString = cols[0], cols[1],cols[3]
        
            if (chro, coordString) not in exonReferenceCoords:
                continue
                
            if gene not in geneObjects:
                continue
                        
            if len(geneObjects[gene].outlierMutants)<1:
                continue
                
            geneObjects[gene].addSpladderEvent(cols,exonReferenceCoords[(chro,coordString)])


    # QC plot
    eventsPerGene = np.asarray([len(geneObj.spladderEvents) for gene, geneObj in geneObjects.items()
                                if len(geneObj.outlierMutants)>0], dtype=np.int)
    title = "Avg. events per gene:%s, Median events per gene:%s" % (np.mean(eventsPerGene), np.median(eventsPerGene))
    ylab, xlab = "count", "spladder events per gene"
    saveAs = 'spladder_events_per_gene.pdf'
    qcPlot(eventsPerGene, "gray", "hist", title, xlab, ylab, saveAs)
    
    return geneObjects, headerDict


def getHist(hFile):

    histoDict = dict()
    
    with open(hFile,'r') as lines:
        for num,line in enumerate(lines,0):                                                                                            
            donorID,histo = line.rstrip().split()                                                                                                 
            
            if histo not in histoDict:                                                                                                       
                histoDict[histo] = set()
            
            histoDict[histo].add(donorID)
    return histoDict


########################################################################
# Main
# Here is the main program
# 
########################################################################

def main():
    '''
    TBD
    '''
    myCommandLine = CommandLine()
    exonOverlapFile = '/pod/pstore/groups/brookslab/csoulette/PCAWG/larvaAnalysis/allEvents_ranges.sorted.merge.bed'
    mutationFile = '/pod/pstore/groups/brookslab/csoulette/PCAWG/larvaAnalysis/allEvents_snv_outlier_hugo_bp_filtered.tsv'
    #mutationFile = 'temp'
    spladderFile = '/pod/pstore/groups/brookslab/csoulette/PCAWG/SplAdder/tables/allEvents_allOut.txt'
    histologyFile = '/pod/pstore/groups/brookslab/csoulette/PCAWG/metadata/donorID_histotype_ref.tsv'
    
    lowerBound = myCommandLine.args['lower_bound']
    upperBound = myCommandLine.args['upper_bound']


    # The general workflow for this process is as follows:
    ## 1. Identify somatic mutations of interest using some criteria
    ## 2. Compute z-scores for all exons in genes harboring somatic mutations.
    ## 3. Return data: gene, exon, sample, sample_group, mutation_criteria, z-score, psi, psi - groupMean_psi.

    
    #Step1. Link spladder events to unique exon coordinates.
    exonReferenceCoords = associateUniqueCoords(exonOverlapFile)
    
    #Step2. Get mutations of interest.
    geneObjects = filterMutations(mutationFile, exonReferenceCoords, lowerBound, upperBound)
    
    #Step3. Define outliers for each gene.
    geneObjects = defOutliers(geneObjects)


    #Step4a. Get all spladder events for genes of interest.
    geneObjects, spladderHeaderDict = addSpladderEvents(spladderFile, exonReferenceCoords, geneObjects)

    #Step4b. Create donor-histology reference dict.
    histoDict = getHist(histologyFile)
        
    #Step5. Calculate z-score and output to tsv.
    [ geneObj.calcEventZ(spladderHeaderDict, histoDict) for gene, geneObj in geneObjects.items() ]

    # Print out.
    for gene, geneObj in geneObjects.items():
        for exon, donors in geneObj.exons.items():
            for donor in donors:
                print(gene, geneObj.hugo, "\t".join(str(x) for x in geneObj.exons[exon][donor]), sep="\t")


if __name__ == "__main__":
    main();
    
