#!/usr/bin/env python3


########################################################################
# File: mutSS.py
#  executable: mutSS.py
# Purpose: Associate somatic variants with changes in splicing.
#
#          
# Author: Cameron M. Soulette
# History:      cms 09/04/2017 Created
#
# This program was written in emacs.
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################
from subprocess import Popen
from multiprocessing import Pool
from computeSSEventDist import calcDist
from computeSSEventDist import checkAnnotation
from computeSSEventDist import zScoreIt

import os, sys
import numpy as np
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
        self.parser = argparse.ArgumentParser(description = '''mutSS.py - a filtering tool to associate somatic variants with splicing quantification data.''',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s ')
        # Add args
        self.parser.add_argument('-vm', '--variant_manifest', action = 'store', required=True, help='Tab delimited file containing group label from which variants are derived from, followed by the direct path to variant file. [Default : req*]')

        self.parser.add_argument('-gtf', '--gtf_annotation', action = 'store', required=True,  help='Annotation GTF File. Used as to determine known/novel junctions, and used to define "exonic regions". [Default : req*]')

        self.parser.add_argument('-s', '--sample_groups', action = 'store', required=True,  help='Table delimited file containing sample ID followed by group. [Default : req*]')

        self.parser.add_argument('-p', '--num_threads', action = 'store', required=False, default=2, type=int,  help='Number of threads to run. [Default : 2]')

        self.group = self.parser.add_mutually_exclusive_group(required = True)

        self.group.add_argument('-spl', '--spladder_table', action = 'store',help='SplAdder formatted splicing quantification table. [Default : req*]')
        self.group.add_argument('-jb', '--juncbase_table', action = 'store', help='JuncBase formatted splicing quantification table. [Default : req*]')
        
        
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))


###

def getHisto(groupRef):
    histoDict = dict()
    with open(groupRef, 'r') as lines:
        for line in lines:
            line = line.rstrip()
            donorID,histo = line.split()
            if histo not in histoDict:
                histoDict[histo] = set()
            
                
            histoDict[histo].add(donorID)

    return histoDict


def runCMD(cmd):
    

    inp, out = cmd
    with open(out,'w') as out1:
        p = Popen(inp, shell=True, stdout=out1)
        p.wait()

    

def readOutputs(cmd):
    inp, out = cmd
    tempDict = dict()
    group = (out.split("_"))[0]
    
    with open(out,'r') as lines:
        for line in lines:
            cols = line.rstrip().split()
            chromo, ensembl = cols[0], cols[6]
            if chromo not in tempDict:
                tempDict[chromo] = dict()
            if ensembl not in tempDict[chromo]:
                tempDict[chromo][ensembl] = list()
                
            tempDict[chromo][ensembl] = line.rstrip()

    return (group, tempDict)
    

def getGenicSNVs(genesByChr, snvManifest, threads):
    
    
    allGenesBed = 'temp_Genes.bed'
    with open(allGenesBed, 'w') as out1:
        for chro, genes in genesByChr.items():
            for gene, data in genes.items():
                start = min(data['starts'])
                stop = max(data['ends'])
                hugoID, chromosome, strand = data['info']
                print(chromosome, start, stop, gene, hugoID, strand, sep="\t", file=out1)

    sortedAllGenesBed = 'temp_sorted_Genes.bed'
    with open(sortedAllGenesBed,'w') as out2:
        p = Popen('bedtools sort -i %s' % allGenesBed, shell=True, stdout=out2)
        p.wait()

    
    cmdList = list()
    intersectCommand = 'bedtools intersect -wa -wb -b %s -a ' % sortedAllGenesBed 
    
    

    with open(snvManifest,'r') as lines:
        for line in lines:
            group, snvFile = line.rstrip().split()
            groupOut = "groupTemp_%s_genesIntersect.bed" % group

            cmdList.append((intersectCommand + snvFile, groupOut)) 
    
    p = Pool(threads)
    
    for i, _ in enumerate(p.imap_unordered(runCMD, cmdList), 1):
        sys.stderr.write('Intersecting SNVs with gene bins...\rdone {0:%}'.format(i/len(cmdList)))
    sys.stderr.write('\n')

        
    groupSNVs = p.map(readOutputs, cmdList)
    p.close()
    return cmdList #groupSNVs

def readGTF(gtf):
    
    '''
    Takes GTF as argument. Returns GTF data structured as dictionary.
    Also returns file from gtf_parser with annotated exon events.
 
    Dictionary includes GTF parsed information about exons for each gene.
    
    genesByChr contains gene exon info from GTF file. 
    1. First key value pair is chromosome, geneID. 
    2. Second key value pair is geneID, data type. (Data type can be end or start position information)
    3. Third key valye pair is data type, list/set of data. (Data can be end or or start position information)

    
    '''

    genesByChr = dict()
    with open(gtf, 'r') as lines:
        for line in lines:
            if line[0] == '#':
                continue

            cols = line.rstrip().split("\t")

            if cols[2] != 'exon':
                continue
                
            chromosome, start, end, strand, dataCol = cols[0], int(cols[3])-1, int(cols[4]), cols[6], cols[-1]
            
            transcriptType = (re.search('transcript_type \"([^"]+)', dataCol)).group(1) 
            
            

            
            geneID = (re.search('(ENSG[^"]+)', dataCol)).group(1) 
            hugoID =(re.search('gene_name \"([^"]+)', dataCol)).group(1) 

            if chromosome not in genesByChr:
                genesByChr[chromosome] = dict()

            if geneID not in genesByChr[chromosome]:
                genesByChr[chromosome][geneID] = {"starts":list(),
                                                  "ends":list(),
                                                  "info": (hugoID, chromosome, strand),
                                                  "allExons":set(),
                                                  "nonIRExons":set()}
                
            genesByChr[chromosome][geneID]["starts"].append(start)
            genesByChr[chromosome][geneID]["ends"].append(end)
            genesByChr[chromosome][geneID]["allExons"].add( (start, end) )
            
            if transcriptType != 'retained_intron':
                genesByChr[chromosome][geneID]["nonIRExons"].add( (start, end) )
            

    annotatedEvents = 'temp_annotatedExons.bed'
    with open(annotatedEvents, 'w') as out1:

        p = Popen('pyton3 gtf_parser.py -i %s -f "exon events"' % (gtf), shell=True, stdout=out1)
        p.wait()
        

            
    return genesByChr, annotatedEvents
            

def snvSplicingEventIntersect(snvGeneIntersectData, splicingEventFile, threads):

    
    dataList = list()
    cmdList = list()
    cmdList_nonIntersecting = list()
    for geneIntersect in snvGeneIntersectData:
        geneIntersectCMD, intersectOut = geneIntersect
        
        group = (intersectOut.split("_"))[1]
        

        cmd = 'bedtools intersect -wa -wb -a %s -b %s' % (intersectOut, splicingEventFile)
        outFile = "groupTemp_%s_genesIntersect_splicingEventIntersect.bed" % group
        cmdList.append( (cmd, outFile ) )

        cmd_nonIntersection = 'bedtools intersect -v -a %s -b %s' % (intersectOut, splicingEventFile)
        outFile_nonIntersection = "groupTemp_%s_genesIntersect_NOsplicingEvent.bed" % group
        cmdList.append( (cmd_nonIntersection, outFile_nonIntersection ) )
        

        dataList.append( (group, intersectOut, splicingEventFile, outFile, outFile_nonIntersection) )



    p = Pool(threads)
    
    # Code for progress bar retrieved from here - goo.gl/wvV6q5
    for i, _ in enumerate(p.imap_unordered(runCMD, cmdList), 1):
        sys.stderr.write('Intersecting SNVs with splicing events...\rdone {0:%}'.format(i/len(cmdList)))
    sys.stderr.write('\n')
    p.close()

    return dataList

def exonSNVAssociate(bedTuple):
    
    group, bedFile, genes = bedTuple

    

    tempSNVDict = dict()
    temp_out = "groupTemp_%s_snvOut.txt" % group
    out1 = open(temp_out, 'w') 
    with open(bedFile, 'r') as lines:
        for line in lines:
            cols = line.rstrip().split()
            
            # SNP Data
            chromosome, snpC1, snpC2, snpData = cols[:4]
            donorID, histo, null, snpType = snpData.split(":")
            refAllele, refAllele, mutAllele = snpType.split(",")

            # Gene Data
            ensembl, hugo, strand = cols[7:10]

            # Splicing Data
            eventCoords, eventType = cols[13:15]

            # annotated Exons
            jcnTupleSet = set(genes[chromosome][ensembl]['starts'])
            jcnTupleSet = jcnTupleSet.union(set(genes[chromosome][ensembl]['ends']))

            ####
            if ensembl not in tempSNVDict:
                tempSNVDict[ensembl] = dict()
            
            if eventCoords not in tempSNVDict[ensembl]:
                tempSNVDict[ensembl][eventCoords] = dict()

            if donorID not in tempSNVDict[ensembl][eventCoords]:
                tempSNVDict[ensembl][eventCoords][donorID] = list()
            ###


            relativeDistance, closestFeature, affectedExon = calcDist(eventType, eventCoords, snpC1, strand)
            checkedExons = checkAnnotation(jcnTupleSet, eventCoords, eventType, strand, ensembl)

            dataTuple = (chromosome, snpC1, snpC2, ".", 
                         strand, donorID, histo,
                         refAllele, mutAllele, ensembl, 
                         hugo, eventCoords, ":".join(str(x) for x in affectedExon), "\t".join(checkedExons),
                         eventType, closestFeature, relativeDistance)
            
            print("\t".join(str(x) for x in dataTuple), file=out1) 
            tempSNVDict[ensembl][eventCoords][donorID].append(dataTuple)
    out1.close()
    

def appendZScores(cmd):
    
    group, groupSampleSet, snvOutFile, splicingQuantFile = cmd

    tmpEventDict = dict()
    
    with open(snvOutFile,'r') as lines:
        for line in lines:
            cols = line.rstrip().split()
            
            ensembl, eventCoords = cols[9], cols[11]
            donorID = cols[5]

            if ensembl not in tmpEventDict:
                tmpEventDict[ensembl] = dict()
            if eventCoords not in tmpEventDict[ensembl]:
                tmpEventDict[ensembl][eventCoords] = dict()

            tmpEventDict[ensembl][eventCoords][donorID] = line.rstrip()


    

    with open(splicingQuantFile, 'r') as lines:
        headerList = next(lines).rstrip().split()
        headerDict = {value:pos for pos,value in enumerate(headerList,0)}
        allGroupSet = set(headerList)
        groupSubSet = groupSampleSet.intersection(allGroupSet)

        if len(groupSubSet)<10:
            return

        out1 = open("groupTemp_%s_snvOut_zscore.tsv" % group, 'w')
        
        for line in lines:
            cols = line.rstrip().split()
            ens = cols[0]

            if ens not in tmpEventDict:
                continue

            eventCoords = cols[4]

            if eventCoords not in tmpEventDict[ens]:
                continue

            psiArray = np.asarray([cols[headerDict[sample]] for sample in groupSubSet], 
                                  dtype=np.float32)
            
            percentNan = float(np.count_nonzero(~np.isnan(psiArray))) / len(psiArray)
            
            if percentNan < 0.7:
                zscores = [np.nan for i in psiArray]
                mean = np.nan
                stdev = np.nan
            else:
                zscores, mean, stdev = zScoreIt(psiArray)

            
            mutantPatients = set(tmpEventDict[ens][eventCoords].keys())
            mutantsWithData = mutantPatients.intersection(groupSubSet)
            
            for i in list(mutantsWithData):
                indexNum = list(groupSubSet).index(i)
                psi = psiArray[indexNum]
                zScore = zscores[indexNum]
                print(tmpEventDict[ens][eventCoords][i], psi, mean, stdev, zScore,  sep="\t", file=out1 )
            
    out1.close()


def intersectNonEventSNVs(cmd):
    group, nonIntersectingSNVFile, exonBed = cmd

    temp_out = "groupTemp_%s_closestExon_nonEventSNV.bed" % group
    with open(temp_out, 'w') as out1:
        p = Popen('bedtools intersect -wa -wb -a %s -b %s' % (nonIntersectingSNVFile, exonBed), shell=True, stdout=out1)
        p.wait()
    
    uniqueDict = dict()
    with open(temp_out, 'r') as lines:
        for line in lines:
            cols = line.rstrip().split()

            # SNP Data
            chromosome, snpC1, snpC2, snpData = cols[:4]
            donorID, histo, null, snpType = snpData.split(":")
            refAllele, refAllele2, mutAllele = snpType.split(",")
            
            # Gene Data
            ensembl, hugo, strand = cols[7:10]

            # Splicing Data
            eventCoords, eventType = cols[13], 'es'


            relativeDistance, closestFeature, affectedExon = calcDist(eventType, eventCoords, snpC1, strand)
            
            dataTuple = (chromosome, snpC1, snpC2, ".", 
                         strand, donorID, histo,
                         refAllele, mutAllele, ensembl, 
                         hugo, eventCoords, ":".join(str(x) for x in affectedExon), 'K', 'K',
                         eventType, closestFeature, relativeDistance)


            key = (donorID, chromosome, snpC1)
            
            if key not in uniqueDict:
                uniqueDict[key] = dataTuple
                
            elif abs(uniqueDict[key][-1])>abs(relativeDistance):
                uniqueDict[key] = dataTuple



    with open('groupTemp_%s_nonIntersectingSNVs_exonIntersect.bed' % group, 'w') as out1:
        for i,v in uniqueDict.items():
            print("\t".join(str(x) for x in v), file=out1)
            

            

def associateSNVs(dataListOfTuples, genes, splicingQuantFile, histoDict, annotatedEvents, threads):


    cmdList = list()
    for groupData in dataListOfTuples:
        group, geneIntersectFile, splicingEventFile, splicingIntersectFile, nonIntersectingSNVFile = groupData
        #exonSNVAssociate(splicingIntersectFile)
        cmdList.append( (group, splicingIntersectFile, genes ) )
        

    p = Pool(threads)
    # Code for progress bar retrieved from here - goo.gl/wvV6q5
    for i, _ in enumerate(p.imap_unordered(exonSNVAssociate, cmdList), 1):
        sys.stderr.write('Computing SNV distance from closest feature...\rdone {0:%}'.format(i/len(cmdList)))
    sys.stderr.write('\n')
    p.close()

    newCMDList = list()
    for cmd in cmdList:
        group, splicingIntersectFile, genes = cmd
        snvOutFile = "groupTemp_%s_snvOut.txt" % group
        groupSampleSet = histoDict[group]

        newCMDList.append( (group, groupSampleSet, snvOutFile, splicingQuantFile) )


    p = Pool(threads)
    # Code for progress bar retrieved from here - goo.gl/wvV6q5
    for i, _ in enumerate(p.imap_unordered(appendZScores, newCMDList), 1):
        sys.stderr.write('Appending z-scores to splicing events...\rdone {0:%}'.format(i/len(newCMDList)))
    sys.stderr.write('\n')
    p.close()
        
        
    cmdList = list()
    for groupData in dataListOfTuples:
        group, geneIntersectFile, splicingEventFile, splicingIntersectFile, nonIntersectingSNVFile = groupData
        cmdList.append( (group, nonIntersectingSNVFile, annotatedEvents) )

    
    p = Pool(threads)
    # Code for progress bar retrieved from here - goo.gl/wvV6q5
    for i, _ in enumerate(p.imap_unordered(intersectNonEventSNVs, cmdList), 1):
        sys.stderr.write('Intersecting splicing non-event associated SNVs...\rdone {0:%}'.format(i/len(cmdList)))
    sys.stderr.write('\n')
    p.close()
    


####

def main():


    myCommandLine = CommandLine()
    
    snvGroupManifest  = myCommandLine.args['variant_manifest']
    splicingEventFile = myCommandLine.args['spladder_table']
    gtfFile           = myCommandLine.args['gtf_annotation']
    histoReference    = myCommandLine.args['sample_groups']
    pthreads          = myCommandLine.args['num_threads']
    

    # Create gene bins.
#    gtfFile = '/pod/pstore/groups/brookslab/csoulette/annotations/gencode.v19.annotation.gtf'
    genes, annotatedEventFile = readGTF(gtfFile)


    # Intersect SNVs with gene Bins.
#    snvGroupManifest = '/pod/pstore/groups/brookslab/PCAWG/Oct2016_Freeze/snv_by_histo_new/histo_snv_manifest.txt'
    snvsGeneIntersections = getGenicSNVs(genes, snvGroupManifest, pthreads)

    # Intersect SNVs with splicing quantifications.
#    splicingEventFile = '/pod/pstore/groups/brookslab/csoulette/PCAWG/SplAdder/tables/sorted_allEvents_allOut_170809.txt'
    intersectionDataList = snvSplicingEventIntersect(snvsGeneIntersections, splicingEventFile, pthreads)
    
    # Get patients per histology
#    histoReference = "/pod/pstore/groups/brookslab/csoulette/PCAWG/metadata/donorID_histotype_ref.tsv"
    histoDict = getHisto(histoReference)


    # Associate SNVs 
#    splicingQuantFile = '/pod/pstore/groups/brookslab/csoulette/PCAWG/SplAdder/tables/allEvents_allOut_170809.tsv'
    associateSNVs(intersectionDataList, genes, splicingEventFile, histoDict, annotatedEventFile, pthreads)


if __name__ == "__main__":
    main()


