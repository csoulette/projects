from subprocess import Popen
from multiprocessing import Pool
from computeSSEventDist import calcDist
from computeSSEventDist import checkAnnotation
from computeSSEventDist import zScoreIt

import os, sys
import numpy as np
import re



###

def getHisto():
    histoDict = dict()
    with open("/pod/pstore/groups/brookslab/csoulette/PCAWG/metadata/donorID_histotype_ref.tsv",'r') as lines:
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
    

def getGenicSNVs(genesByChr, snvManifest):
    
    
    allGenesBed = 'tempGenes.bed'
    with open(allGenesBed, 'w') as out1:
        for chro, genes in genesByChr.items():
            for gene, data in genes.items():
                start = min(data['starts'])
                stop = max(data['ends'])
                hugoID, chromosome, strand = data['info']
                print(chromosome, start, stop, gene, hugoID, strand, sep="\t", file=out1)

    sortedAllGenesBed = 'sortedTempGenes.bed'
    with open(sortedAllGenesBed,'w') as out2:
        p = Popen('bedtools sort -i %s' % allGenesBed, shell=True, stdout=out2)
        p.wait()

    
    cmdList = list()
    intersectCommand = 'bedtools intersect -wa -wb -b %s -a ' % sortedAllGenesBed 
    
    

    with open(snvManifest,'r') as lines:
        for line in lines:
            group, snvFile = line.rstrip().split()
            groupOut = "%s_genesIntersect.bed" % group

            cmdList.append((intersectCommand + snvFile, groupOut)) 
    
    p = Pool(30)
    
    for i, _ in enumerate(p.imap_unordered(runCMD, cmdList), 1):
        sys.stderr.write('Intersecting SNVs with gene bins...\rdone {0:%}'.format(i/len(cmdList)))
    sys.stderr.write('\n')

        
    groupSNVs = p.map(readOutputs, cmdList)
    p.close()
    return cmdList #groupSNVs

def readGTF(gtf):
    

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
            
            
    return genesByChr
            

def snvSplicingEventIntersect(snvGeneIntersectData, splicingEventFile):

    
    dataList = list()
    cmdList = list()
    cmdList_nonIntersecting = list()
    for geneIntersect in snvGeneIntersectData:
        geneIntersectCMD, intersectOut = geneIntersect
        
        group = (intersectOut.split("_"))[0]
        

        cmd = 'bedtools intersect -wa -wb -a %s -b %s' % (intersectOut, splicingEventFile)
        outFile = "%s_genesIntersect_splicingEventIntersect.bed" % group
        cmdList.append( (cmd, outFile ) )

        cmd_nonIntersection = 'bedtools intersect -v -a %s -b %s' % (intersectOut, splicingEventFile)
        outFile_nonIntersection = "%s_genesIntersect_NOsplicingEvent.bed" % group
        cmdList.append( (cmd_nonIntersection, outFile_nonIntersection ) )
        

        dataList.append( (group, intersectOut, splicingEventFile, outFile, outFile_nonIntersection) )



    p = Pool(30)
    
    # Code for progress bar retrieved from here - goo.gl/wvV6q5
    for i, _ in enumerate(p.imap_unordered(runCMD, cmdList), 1):
        sys.stderr.write('Intersecting SNVs with splicing events...\rdone {0:%}'.format(i/len(cmdList)))
    sys.stderr.write('\n')
    p.close()

    return dataList

def exonSNVAssociate(bedTuple):
    
    group, bedFile, genes = bedTuple

    

    tempSNVDict = dict()
    temp_out = "%s_snvOut.txt" % group
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

        out1 = open("%s_snvOut_zscore.tsv" % group, 'w')
        
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

    temp_out = "%s_closestExon_nonEventSNV.bed" % group
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



    with open('%s_nonIntersectingSNVs_exonIntersect' % group, 'w') as out1:
        for i,v in uniqueDict.items():
            print("\t".join(str(x) for x in v), file=out1)
            

            

def associateSNVs(dataListOfTuples, genes, splicingQuantFile, histoDict, annotatedEvents):


    cmdList = list()
    for groupData in dataListOfTuples:
        group, geneIntersectFile, splicingEventFile, splicingIntersectFile, nonIntersectingSNVFile = groupData
        #exonSNVAssociate(splicingIntersectFile)
        cmdList.append( (group, splicingIntersectFile, genes ) )
        

    p = Pool(30)
    # Code for progress bar retrieved from here - goo.gl/wvV6q5
    for i, _ in enumerate(p.imap_unordered(exonSNVAssociate, cmdList), 1):
        sys.stderr.write('Computing SNV distance from closest feature...\rdone {0:%}'.format(i/len(cmdList)))
    sys.stderr.write('\n')
    p.close()

    newCMDList = list()
    for cmd in cmdList:
        group, splicingIntersectFile, genes = cmd
        snvOutFile = "%s_snvOut.txt" % group
        groupSampleSet = histoDict[group]

        newCMDList.append( (group, groupSampleSet, snvOutFile, splicingQuantFile) )


    p = Pool(30)
    # Code for progress bar retrieved from here - goo.gl/wvV6q5
    for i, _ in enumerate(p.imap_unordered(appendZScores, newCMDList), 1):
        sys.stderr.write('Appending z-scores to splicing events...\rdone {0:%}'.format(i/len(newCMDList)))
    sys.stderr.write('\n')
    p.close()
        
        
    cmdList = list()
    for groupData in dataListOfTuples:
        group, geneIntersectFile, splicingEventFile, splicingIntersectFile, nonIntersectingSNVFile = groupData
        cmdList.append( (group, nonIntersectingSNVFile, annotatedEvents) )

    
    p = Pool(30)
    # Code for progress bar retrieved from here - goo.gl/wvV6q5
    for i, _ in enumerate(p.imap_unordered(intersectNonEventSNVs, cmdList), 1):
        sys.stderr.write('Intersecting splicing non-event associated SNVs...\rdone {0:%}'.format(i/len(cmdList)))
    sys.stderr.write('\n')
    p.close()
    


####

def main():


    # Create gene bins.
    gtfFile = '/pod/pstore/groups/brookslab/csoulette/annotations/gencode.v19.annotation.gtf'
    genes = readGTF(gtfFile)


    # Intersect SNVs with gene Bins.
    snvGroupManifest = '/pod/pstore/groups/brookslab/PCAWG/Oct2016_Freeze/snv_by_histo_new/histo_snv_manifest.txt'
    snvsGeneIntersections = getGenicSNVs(genes, snvGroupManifest)

    # Intersect SNVs with splicing quantifications.
    splicingEventFile = '/pod/pstore/groups/brookslab/csoulette/PCAWG/SplAdder/tables/sorted_allEvents_allOut_170809.txt'
    intersectionDataList = snvSplicingEventIntersect(snvsGeneIntersections, splicingEventFile)
    
    # Get patients per histology
    histoDict = getHisto()


    # Associate SNVs 
    annotatedEvents = 'sorted.annotated_events.bed'
    splicingQuantFile = '/pod/pstore/groups/brookslab/csoulette/PCAWG/SplAdder/tables/allEvents_allOut_170809.tsv'
    associateSNVs(intersectionDataList, genes, splicingQuantFile, histoDict, annotatedEvents)


if __name__ == "__main__":
    main()


