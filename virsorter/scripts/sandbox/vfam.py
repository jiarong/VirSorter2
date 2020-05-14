#!/usr/bin/env python

"""
SCRIPT:  profileHMMsFromFASTA.py
AUTHOR:  Peter Skewes-Cox
UPDATED:  February 2014

DESCRIPTION:
This script takes as input a FASTA file containing protein
sequences from which profile-HMMs are ultimately built.
Redundant sequences are collapsed, and are BLASTed in
an 'all-by-all' fashion. The BLAST results are used as
input to MCL to generate clusters of related sequences. These
sequence clusters can be filtered further, and multiple
sequence alignments are generated for each cluster. The MSAs
are input into hmmbuild, and the resulting profile-HMMs are
concatenated into a single HMMER3-compatible profile-HMM.

USAGE:    
  profileHMMsFromFASTA.py -f <FASTA> [-m <#>] [-p <#>] [-c <#>] [-C] [-I <#>] [-n <#>] [-a <#>] [-o <string>] [-h]

OPTIONS:
  -f      FASTA file sequences from which to generate clusters
  -m      minimum sequence length [default = 1]
  -p      minumum fraction identity for initial sequence collapsing [default = 1.0]
  -c      minimum fraction coverage for initial sequence collapsing [default = 0.0]
  -C      impose fraction coverage heuristics for inclusion of sequences in MSAs
  -I      inflation number for cluster expansion in mcl [default = 2.0]
  -n      minimum number of sequences allowed in a MSA [default = 2]
  -a      number of cores on which to run all processes [default = 8]
  -o      output prefix for cluster names (default based on filename)
  -h      print help message

"""

import sys,getopt,os
import errno

def main(args=None):

    fastaFile = None
    minSequenceLength = 1
    fractionID = 1.0
    fractionCov = None
    coverageHeuristics = False
    inflationNum = None
    minSequences = 2
    numCores = 8
    prefix = None
    shortOptions = "f:m:p:c:CI:n:a:o:h"
    longOptions = []

    try:
        opts, args = getopt.getopt(args, shortOptions, longOptions)
    except getopt.error as msg:
        print(msg)
        print(__doc__)
        return(-2)

    try:
        for option,optionArg in opts:
            if option=='-f':
                fastaFile = optionArg
            elif option=='-m':
                minSequenceLength = int(optionArg)
            elif option=='-p':
                fractionID = float(optionArg)
            elif option=='-c':
                fractionCov = float(optionArg)
            elif option=='-C':
                coverageHeuristics = True
            elif option=='-I':
                inflationNum = float(optionArg)
            elif option=='-n':
                minSequences = int(optionArg)
            elif option=='-a':
                numCores = int(optionArg)
            elif option=='-o':
                prefix = optionArg
            elif option=='-h':
                print(__doc__)
                return(0)
            else:
                print("%s option not implemented" % option)

    except Exception as msg:
        print(msg)
        print(__doc__)
        return(-2)
    
    argProblem = False

    if len(args) > 0:
        print("What are these extra arguments?: %s" % (' '.join(args)))
        argProblem = True
    if fastaFile == None:
        print("You must provide a FASTA file!")
        argProblem = True
    elif not os.path.isfile(fastaFile):
        print("%s does not exist!" % fastaFile)
        argProblem = True
    if minSequenceLength < 1:
        print("Minimum sequence length must be at least 1!")
        argProblem = True
    if fractionID < 0.0 or fractionID > 1.0:
        print("Fraction identity must be between 0 and 1 inclusive!")
        argProblem = True
    elif fractionID < 0.5:
        print("WARNING: CD-HIT recommends not dropping below 50% identity for sequence collapsing!")
    if fractionCov != None and (fractionCov < 0.0 or fractionCov > 1.0):
        print("Fraction coverage must be between 0 and 1 inclusive!")
        argProblem = True
    if inflationNum != None:
        if inflationNum < 0.0:
            print("The inflation number must be positive!")
            argProblem = True
        elif inflationNum <= 1.0 or inflationNum > 6.0:
            print("WARNING: The recommended inflation values are > 1.0 and <= 6.0!")
    if minSequences < 1:
        print("The minimum number of sequences must be at least 1 to build a multiple sequence alignment!")
        argProblem = True
    elif minSequences == 1:
        print("Though you can technically build a multiple sequence alignment from a single sequence, you may want to reconsider!")
    if numCores < 1:
        print("The number of cores must be at least 1!")
        argProblem = True
    if prefix != None and '.' in prefix:
        print("Do not use a period ('.') in your prefix!")
        argProblem = True
    if argProblem:
        print(__doc__)
        return(-1)

    # first make sure we have all the external dependencies satisfied
    print("Looking for necessary executables...")
    missingExecutables = testExecutables(["cd-hit","makeblastdb","blastp",["mcxload","-h"],"mcl","hmmbuild"])
    if len(missingExecutables) != 0:
        for missing in missingExecutables:
            print("Unable to call %s!  This will not work without this dependency resolved!" % missing)
        sys.exit(-2)
    # then let the pipeline ride
    print("Removing exact sequence duplicates from FASTA file...")
    filteredFASTA = filterFASTA(fastaFile,minSequenceLength)
    print("Collapsing sequences at specified redundancy levels...")
    collapsedFASTA = collapseSequence(filteredFASTA,fractionID,fractionCov=fractionCov)
    print("Performing all-by-all BLAST of collapsed FASTA file...")
    blastResults = allByAllBLAST(collapsedFASTA,numCores)
    print("Filtering out polyprotein sequences...")
    blackList = filterPolyproteins(blastResults)
    # blastToMcl has its own verbosity to report status messages
    mclResults = blastToMcl(blastResults,blackList,inflationNum=inflationNum)
    print("Splitting FASTA file based on MCL results...")
    fastaFiles = mclToFASTA(mclResults,collapsedFASTA,prefix=prefix)
    print("Filtering cluster FASTA files...")
    filteredFiles = filterOnCoverage(fastaFiles,coverageHeuristics=coverageHeuristics)
    unalignedFiles = filterOnNumber(filteredFiles,minSequences=minSequences)
    print("Running MUSCLE on all cluster FASTA files...")
    alignedFiles = batchMuscleCall(unalignedFiles,prefix=prefix)
    print("Running hmmbuild on all alignment files...")
    hmmFiles = batchHmmbuildCall(alignedFiles,numCores=numCores,prefix=prefix)
    # get ready to concatenate the final output HMM
    fastaParts = os.path.splitext(fastaFile)
    if prefix == None:
        outputFile = fastaParts[0]+".hmm"
    else:
        outputFile = prefix+".hmm"
    print("Building master HMM file...")
    concatentateHMMs(hmmFiles,outputFile)
    print("%s contains profile-HMMs from all input clusters." % outputFile)
    print("Creating directories and moving files...")

    # do some clean-up of the intermediate files
    allFastaFiles = set()
    allFastaFiles.add(filteredFASTA)
    allFastaFiles.add(collapsedFASTA)
    for ff in fastaFiles: allFastaFiles.add(ff)
    for ff in filteredFiles: allFastaFiles.add(ff)
    for ff in unalignedFiles: allFastaFiles.add(ff)
    try:
        os.mkdir("msaFiles")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for m in alignedFiles:
        os.rename(m,"msaFiles/"+m)

    try:
        os.mkdir("fastaFiles")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for f in allFastaFiles:
        os.rename(f,"fastaFiles/"+f)

    try:
        os.mkdir("hmmFiles")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for h in hmmFiles:
        os.rename(h,"hmmFiles/"+h)

    try:
        os.mkdir("logFiles")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for l in [x for x in os.listdir(os.getcwd()) if (x.startswith(fastaParts[0]) or (prefix != None and x.startswith(prefix))) and x.endswith(".log")]:
        os.rename(l,"logFiles/"+l)

    try:
        os.mkdir("clusteringFiles")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for c in [x for x in os.listdir(os.getcwd()) if x.startswith(fastaParts[0]) and \
                    (("allByAll" in x) or (x.endswith(".clstr")) or \
                      (x[-4:] in ('.pin','.phr','.psq')))]:
        os.rename(c,"clusteringFiles/"+c)
    
    print("Done!")
    
    return(0)

def testExecutables(executables):
    """
    This function takes in a list of executables whose presence
    is to be detected before ostensibly making system calls 
    to them at a later time.  It returns a list of executables
    that were NOT identified as callable.
    """
    from subprocess import call,STDOUT
    from os import devnull
    badExecutables = []
    for executable in executables:
        try:
            returnVal = call(executable,stdout=open(devnull,'w'),stderr=STDOUT)
        except OSError:
            badExecutables.append(executable)
    return badExecutables

def filterFASTA(fastaFile,minSequenceLength):
    """
    This function removes exact sequence duplicates from a FASTA file,
    as well as those sequence that do not meet the minimum length.
    """
    fastaFih = FastaIterator(open(fastaFile))
    fastaFnp = os.path.splitext(fastaFile)
    fastaFon = fastaFnp[0]+"_noDupes_minL"+str(minSequenceLength)+fastaFnp[1]
    fastaFoh = open(fastaFon,'w')
    # storing the sequence hashes is far easier than storing the
    # the sequences themselves
    seqHashes = {}
    for record in fastaFih:
        if hash(record.sequence) not in seqHashes:
            if len(record.sequence) >= minSequenceLength:
                fastaFoh.write(str(record)+'\n')
            seqHashes[hash(record.sequence)] = None
    fastaFoh.close()
    return fastaFon

def collapseSequence(fastaFile,fractionID,fractionCov=None):
    """
    This function takes in a FASTA file, a minimum fraction coverage,
    and a minimum fraction idenitity and calls CD-HIT 
    <http://www.bioinformatics.org/cd-hit/> to collapse the input
    sequences into non-redundant representatives at the specified levels.
    It returns a collapsed FASTA file name.
    """
    from subprocess import call
    if fractionCov == None:
        coverageArgs = ''
        coverageTitle = ''
    else:
        coverageArgs = " -s %s " % fractionCov
        coverageTitle = "_s%s" % (str(int(100*fractionCov)))
    outputParts = os.path.splitext(fastaFile)
    identityTitle = "_c%s" % (str(int(100*fractionID)))
    outputFile = outputParts[0] + identityTitle + coverageTitle + outputParts[1]
    cdhitCmd = "cd-hit -i %s -o %s %s -c %s" % (fastaFile,outputFile,coverageArgs,fractionID)
    returnVal = call(cdhitCmd.split())
    return outputFile

def clust(fastaFile,fractionID,fractionCov=None):
    """
    This function takes in a FASTA file, a minimum fraction coverage,
    and a minimum fraction idenitity and calls CD-HIT 
    <http://www.bioinformatics.org/cd-hit/> to collapse the input
    sequences into non-redundant representatives at the specified levels.
    It returns a collapsed FASTA file name.

    word size 5 is for thresholds 0.7 ~ 1.0
    word size 4 is for thresholds 0.6 ~ 0.7
    word size 3 is for thresholds 0.5 ~ 0.6
    word size 2 is for thresholds 0.4 ~ 0.5 (also see psi-cd-hit)
    below 0.4 psi-cd-hit
    """
    from subprocess import call
    if fractionCov == None:
        coverageArgs = ''
        coverageTitle = ''
    else:
        coverageArgs = " -s %s " % fractionCov
        coverageTitle = "_s%s" % (str(int(100*fractionCov)))
    outputParts = os.path.splitext(fastaFile)
    identityTitle = "_c%s" % (str(int(100*fractionID)))
    outputFile = outputParts[0] + identityTitle + coverageTitle + outputParts[1]

    cdhitCmd = "cd-hit -i %s -o %s %s -c %s" % (fastaFile,outputFile,coverageArgs,fractionID)
    returnVal = call(cdhitCmd.split())
    return outputFile

def allByAllBLAST(fastaFile,numCores):
    """
    This function takes in a FASTA file and a number of cores, formats
    the file using formatdb to be a BLAST protein database, and runs
    BLAST on the file to itself as the BLAST-formatted database.  It
    returns an tab-delimited formatted BLAST results file name.
    """
    import os
    from subprocess import call
    formatdbCmd = "makeblastdb -in %s -dbtype prot" % fastaFile
    returnVal = call(formatdbCmd.split())
    blastResults = os.path.splitext(fastaFile)[0] + "_allByAll_blastp.br"
    blastCmd = "blastp -query %s -out %s -db %s -outfmt 6 -num_threads %s" % (fastaFile,blastResults,fastaFile,numCores)
    returnVal = call(blastCmd.split())
    return blastResults    

def filterPolyproteins(blastFile):
    """
    This function takes a BLAST results file as input and uses a number
    of heuristics to make a call on sequences likely to be polyprotein
    sequences. The function returns a 'black list' of proteins to be 
    removed at the next step.
    """
    blastFh = open(blastFile)
    # store all hits in alignInfo, with structure
    # {query:(subject,qStart,qEnd)}
    alignInfo = {}
    # queryInfo infers sequence length from lines where
    # query and subject are one and the same
    queryInfo = {}
    for line in blastFh:
        data = line.split('\t')
        # check for self vs. self blast
        if data[0] == data[1]:
            queryInfo[data[0]] = int(data[3])
        else:
            if data[0] not in alignInfo:
                alignInfo[data[0]] = [(data[1],int(data[6]),int(data[7]))]
            else:
                alignInfo[data[0]].append((data[1],int(data[6]),int(data[7])))
    blackList = []
    for query,queryLength in list(queryInfo.items()):
        # we don't want to remove shorter sequences
        if queryLength > 400 and query in alignInfo:
            queryAlignments = alignInfo[query]
            if len(queryAlignments) > 1:
                # remove alignments where the subject isn't >70% covered
                filteredAlignments = []
                for alignment in queryAlignments:
                    subject = alignment[0]
                    if subject in queryInfo:
                        subjectLength = queryInfo[subject]
                        # don't want to consider alignments where
                        # query and subject have similar lengths
                        if subjectLength < 0.7*queryLength:
                            subjectCoverage = float(abs(alignment[1]-\
                                             alignment[2]))/float(subjectLength)
                            if subjectCoverage >= 0.7:
                                filteredAlignments.append((alignment[1],\
                                                           alignment[2]))
                # blacklist if the query is > 90% covered
                queryCoverage = {}
                for qRange in filteredAlignments:
                    for nt in range(qRange[0],qRange[1]):
                        queryCoverage[nt] = None
                if len(queryCoverage) > 0.8*queryLength:
                    blackList.append(query)
                del queryCoverage

    return blackList

def blastToMcl(blastResults,blackList,inflationNum=None,prefix=None):
    """
    This function takes in an -m 8 formatted BLAST results file, a
    blast list of sequences to be ignored, and an inflation number.
    It converts the good sequences in the BLAST results file to an 
    .abc file, calls mcxload on the .abc file to generate a .mci and
    .tab file, and calls mcl on the latter two files to generate 
    newline-separated clusters containing tab-delimited sequence titles.
    """
    import os
    from subprocess import call
    if prefix == None:
        prefix = os.path.splitext(blastResults)[0]
    abcFilename = prefix + ".abc"
    mciFilename = prefix + ".mci"
    tabFilename = prefix + ".tab"
    blastFile = open(blastResults)
    abcFile = open(abcFilename,'w')
    # convert BLAST results to abc format
    print("Converting BLAST results to .abc format...")
    for line in blastFile:
        data = line.rstrip().split()
        if data[0] not in blackList and data[1] not in blackList:
            abcLine = '\t'.join([data[0],data[1],data[10]])+'\n'
            abcFile.write(abcLine)
    abcFile.close()
    # run mcxload to generate matrix
    print("Running mcxload on .abc file...")
    mcxloadCmd = "mcxload -abc %s --stream-mirror --stream-neg-log10 -stream-tf "'ceil(200)'" -o %s -write-tab %s" % (abcFilename,mciFilename,tabFilename)
    returnVal= call(mcxloadCmd.split())
    mclFilename = prefix + ".mcl"
    if inflationNum == None:
        inflationArgs = ''
    else:
        inflationArgs = '-I %s ' % inflationNum
    print("Running mcl on generate matrix files...")
    mclCmd = "mcl %s -use-tab %s %s -o %s" % (mciFilename,tabFilename,inflationArgs,mclFilename)
    returnVal = call(mclCmd.split())
    return mclFilename

def mclToFASTA(mclFile,fastaFile,prefix=None):
    """
    This function takes in the MCL clusters and a FASTA file
    and outputs a bunch of numbered FASTA files and returns
    a list of the file names which will be used later. The files
    are named based on the input FASTA file or the provided
    prefix.
    """
    fastaParts = os.path.splitext(fastaFile)
    lineNum = 0
    mclDict = {}

    mclFi = open(mclFile)
    for line in mclFi:
        lineNum += 1
        if prefix == None:
            fastaFon = fastaParts[0]+'_cluster'+str(lineNum)+fastaParts[1]
        else:
            fastaFon = prefix+"_"+str(lineNum)+fastaParts[1]
        for title in line.rstrip().split('\t'):
            mclDict[title] = fastaFon

    fastaFih = FastaIterator(open(fastaFile))
    for record in fastaFih:
        if record.title.split()[0] in mclDict:
            fastaFon = mclDict[record.title.split()[0]]
            fastaFoh = open(fastaFon,'a')
            if prefix == None:
                cluster = fastaFon.split("_cluster")[-1].split(".")[0]
                clusterText = "cluster"+cluster
            else:
                clusterText = os.path.splitext(fastaFon)[0]
                recordItems = record.title.split('|')
            recordItems = record.title.split('|')
            recordItems.insert(4,clusterText)
            record.title = '|'.join(recordItems)    
            fastaFoh.write(str(record)+'\n')
            fastaFoh.close()

    return list(set(mclDict.values()))

def filterSequencesOnCoverage(fastaFile,covHeurDict):
    """
    This function takes in a FASTA file and a dictionary containing
    coverage heuristics dictionary of the form {k1:v1,k2:v2,...,kn:vn}
    where kx represents 100 residues, and vx represents the minimum
    fraction relative length between each sequence and the median sequence
    length of the cluster.  It returns None if there exist no sequences
    meeting the minimum/maximum coverage requirements (an admittedly
    extremely rare, yet completely possible condition), and returns the
    name of the filtered FASTA file.
    """
    recordLengths = {}
    fastaRecords = FastaIterator(open(fastaFile))
    for record in fastaRecords:
        recordLengths[record.title] = len(record)
    lengths = list(recordLengths.values())
    lengths.sort()
    lenLengths = len(lengths)
    if lenLengths % 2 != 0:
        median = lengths[(lenLengths-1)/2]
    else:
        upper = lengths[lenLengths/2]
        lower = lengths[(lenLengths/2)-1]
        median = float(upper+lower)/2
    coverageKey = int(median/100)
    # We want bi-directional coverage, so a 60% coverage cut-off
    # means that the length of the smallest sequence is no less
    # than 60% of the length of the largest sequence.  This will
    # be centered around the median.
    if coverageKey in covHeurDict:
        coverageThreshold = (1.0+covHeurDict[coverageKey])/2
    else:
        coverageThreshold = (1.0+max(covHeurDict.values()))/2
    for title,length in list(recordLengths.items()):
            if length < coverageThreshold * median or length * coverageThreshold > median:
                del recordLengths[title]
    if len(recordLengths) == 0:
        return None
    elif len(recordLengths) == lenLengths:
        return fastaFile
    else:
        fastaRecords = FastaIterator(open(fastaFile))
        fastaParts = os.path.splitext(fastaFile)
        newFastaFilename = fastaParts[0]+"_filtered"+fastaParts[1]
        newFastaFile = open(newFastaFilename,'w')
        for record in fastaRecords:
            if record.title in recordLengths:
                newFastaFile.write(str(record)+'\n')
        newFastaFile.close()
        return newFastaFilename

def filterOnCoverage(fastaFiles,coverageHeuristics=False):
    """
    This function takes in a list of FASTA files, and a boolean
    to determine whether or not to apply coverage heuristics. 
    If the coverage heuristics are applied, the sequences that
    don't meet the requirement are filtered out. This function
    returns a list of the filtered FASTA files.
    """    
    # make new FASTA files when necessary for those clusters
    # containing sequences that don't meet coverage requirements
    if coverageHeuristics:
        # these heuristics derive from the FlowerPower publication:
        # http://www.ncbi.nlm.nih.gov/pubmed/17288570
        covHeurDict = {0:0.6,1:0.65,2:0.7,3:0.75,4:0.8,5:0.85}
        newFastaFiles = []
        for fastaFile in fastaFiles:
            newFastaFile = filterSequencesOnCoverage(fastaFile,covHeurDict)
            if newFastaFile != None:
                newFastaFiles.append(newFastaFile)
        fastaFiles = newFastaFiles
    return fastaFiles

def filterOnNumber(fastaFiles,minSequences=2):
    """
    This function takes in a list of FASTA files and a threashold
    for the minumum number of sequences required to keep a cluster.
    Clusters that don't meet the number of sequences requirement
    are removed. This function returns a filtered list of FASTA files.
    """
    if minSequences != 1:
        newFastaFiles = []
        for fastaFile in fastaFiles:
            count = 0
            fastaFh = open(fastaFile)
            fIn = FastaIterator(fastaFh,raw=True)
            for record in fIn:
                count += 1
                if count >= minSequences:
                    newFastaFiles.append(fastaFile)
                    fastaFh.close()
                    break
        fastaFiles = newFastaFiles
    return fastaFiles

def batchMuscleCall(fastaFiles,prefix=None):
    """
    This function takes in a list of FASTA files makes multiple
    sequence alignments using MUSCLE. It returns a list of the
    alignment files generated in FASTA format.
    """
    from subprocess import call

    # do the alignments
    alignedFiles = []
    if prefix == None:
        masterFn = fastaFiles[0].split("_cluster")[0]+"_inProfiles.fasta"
    else:
        masterFn = prefix+"_inProfiles.fasta"
    masterFh = open(masterFn,'w')
    for fastaFile in fastaFiles:
        fastaFh = open(fastaFile)
        for line in fastaFh:
            masterFh.write(line)
        fastaFh.close()
        fastaParts = os.path.splitext(fastaFile)
        alignmentFile = fastaParts[0]+".msa"+fastaParts[1]
        alignedFiles.append(alignmentFile)
        logFile = alignmentFile+".log"
        muscleCmd ="muscle -in %s -out %s -log %s -quiet" % (fastaFile,alignmentFile,logFile) 
        returnVal = call(muscleCmd.split())
    masterFh.close()
    return alignedFiles

def batchHmmbuildCall(alignmentFiles,numCores=16,prefix=None):
    """
    This function takes in a list of FASTA MSA files, and 
    builds HMMs for each of them using HMMer3's hmmbuild 
    function.
    """
    from subprocess import call
    hmmFiles = []
    for alignmentFile in alignmentFiles:
        if prefix == None:
            hmmName = alignmentFile.rsplit('.msa', 1)[0]
            #hmmName = os.path.splitext(alignmentFile)[0]
        else:
            hmmName = alignmentFile.rsplit('.msa', 1)[0]
            #hmmName = '_'.join(alignmentFile.split('_')[0:2]).split('.')[0]
        hmmFile = hmmName+".hmm"
        hmmFiles.append(hmmFile)
        logFile = hmmFile + ".log"
        hmmbuildCmd = "hmmbuild --informat afa -n %s -o %s --cpu %s %s %s" % (hmmName,logFile,numCores,hmmFile,alignmentFile)
        returnVal = call(hmmbuildCmd.split())
    return hmmFiles

def concatentateHMMs(hmmFiles,outputFile):
    """
    This function takes in a list of HMM files and an 
    output file name and creats a single HMM file with 
    the provided name.
    """    
    outputFh = open(outputFile,'w')
    for hmmFile in hmmFiles:
        hmmFh = open(hmmFile)
        for line in hmmFh:
            outputFh.write(line)
        hmmFh.close()
    outputFh.close()

# here is the class and associated functions necessary for parsing
# FASTA sequence files

class FastaRecord:

    def __init__ (self,title='',sequence='',colwidth=60):
        self.colwidth=colwidth
        self.title=title
        self.sequence=sequence

    def __str__(self):
        s = []
        s.append('>%s' % self.title)
        i = 0
        while i < len(self):
            s.append(self[i:i+self.colwidth])
            i = i + self.colwidth
        return os.linesep.join(s)
    
    def __len__(self):
        return len(self.sequence)

    def __getitem__(self,item):
        return self.sequence[item]
   
def FastaIterator(fh,raw=False):

    def readToTitle(fh):
        preLines = []
        while True:
            l = fh.readline()
            if l.startswith('>'):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)

    preLines,nextTitleLine = readToTitle(fh)

    while nextTitleLine != None:
        title = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine = readToTitle(fh)
        if raw:
            yield (title,''.join(preLines))
        else:
            rec = FastaRecord()
            rec.title = title
            rec.sequence = ''.join([x.rstrip() for x in preLines])
            yield rec

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
