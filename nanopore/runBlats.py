#!/usr/bin/env python3.4
import subprocess
from multiprocess import pool
import sys
import random

#manifestFile = sys.argv[1]
#splitNumber = int(sys.argv[2])
reference = "/pod/pstore/groups/brookslab/reference_indices/blat_indices/GRCh38.primary_assembly.chr-only.fa.2bit"
oocFile = "/pod/pstore/groups/brookslab/reference_indices/blat_indices/hg38_11.occ"

class CommandLine(object):

 def __init__(self, inOpts=None) :
        '''                                                                                                                                                                                                                                                                                                    
        CommandLine constructor.                                                                                                                                                                                                                                                                               
        Implements a parser to interpret the command line argv string using argparse.                                                                                                                                                                                                                          
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'runBlats.pt - a tool for splitting fasta files and running multiple instances of blat.',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null',
                                             add_help = True, #default is True                                                                                                                                                                                                                                 
                                             prefix_chars = '-',
                                             usage = '%(prog)s -p N --manifest fastaManifest.tsv'
                                              )
        # Set arguments.                                                                                                                                                                                                                                                                              
        self.parser.add_argument('--manifest', action='store', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Fasta manifest file. Tab delimited file containing condition information in the first column, and direct path to file in second column (Default: stdin)')
        self.parser.add_argument('-p', '--threads', action='store', default=2, type=int, help='Number of blats to run simultaneously  (Default: 2)')
        
        # inOpts can be used to set default arguments and is useful for development/debugging.                                                                                                                                                                                                                 
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

            
        self.manifestFile = self.args['manifest']
        self.threads = self.args['threads']





def readManifest(tsvFile):
    '''
    Takes in tsv file with condition/path delimeted with tab. 
    Each line is a separate fasta file.
    '''
    fileDict = dict()
    with tsvFile as lines:
        for line in lines:
            condition,filePath = line.rstrip().split("\t")
            fileDict[condition] = filePath
            print(filePath)
    return fileDict


def splitAndBlat(fileDict,splitBy):
    '''
    For each file, split it into smaller batches, run blat, then concatenate the results.
    '''

    for condition,directPath in fileDict.items():
        process = subprocess.Popen('grep -c ">" %s ' % directPath, shell=True, stdout=subprocess.PIPE)
        process.wait()
        totalLines = int(process.communicate()[0].decode("utf-8").rstrip())
        fileList = splitFiles(condition,directPath,splitBy)
        
        # Check lines of new files
        totalLinesFromSplit = 0
        for files in fileList:
            process = subprocess.Popen('grep -c ">" %s ' % files, shell=True, stdout=subprocess.PIPE)
            process.wait()
            totalLinesFromSplit += int(process.communicate()[0].decode("utf-8").rstrip())
        
        if totalLines != totalLinesFromSplit:
            sys.exit(1)

        else:
            runBlats(fileList,condition)
            # Cleanup
            for files in fileList:
                process = subprocess.Popen('rm %s ' % files, shell=True)
                process.wait()

def runBlats(fileList,condition):
    
    cmd = "/pod/pstore/groups/brookslab/bin/blat -noHead -ooc=%s -q=rna -fine -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 " % oocFile
    cmd = cmd + "%s " % reference
    waitList = list()
    outBlats = list()
    for file in fileList:
        outFile = file + "blatOut.psl"
        waitList.append(subprocess.Popen("%s %s %s" % (cmd,file,outFile), shell=True, stderr=subprocess.DEVNULL))
        outBlats.append(outFile)
                        
    for process in waitList:
        process.wait()
                        

    subprocess.Popen("cat %s > %s_blat_out.psl" % (" ".join([file for file in outBlats]),condition), shell = True)
                        
    

def splitFiles(filePrefix,filePath,numFiles):
    
    random.seed(100)
    
    fileList = [filePrefix] * numFiles
    fileList = [prefix + "_" + str(random.randint(0,9999999)) for prefix in fileList]
    writeList = [open(outFile,'w') for outFile in fileList]
 
    with open(filePath,'r') as lines:
        faEntry = str()
        
        for line in lines:
            if line[0] == ">":
                if len(faEntry):
                    print(faEntry, file= random.choice(writeList))
                    faEntry = line.rstrip() + "\n"
                else:
                    faEntry = line.rstrip() + "\n"
            else:
                faEntry += line.rstrip()
                
        if len(faEntry):
            print(faEntry, file= random.choice(writeList))
            faEntry = line.rstrip() + "\n"
    [outFile.close for outFile in writeList]

    return fileList


def main():
    '''
    run main program.
    '''
    myCommandLine = CommandLine()
    manifestFile,splitNumber = myCommandLine.manifestFile,myCommandLine.threads
    fileDict = readManifest(manifestFile)
    splitAndBlat(fileDict,splitNumber)



main()
