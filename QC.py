#!/usr/bin/env python

__author__ = "Sivasubramani S"
__email__ = "s.sivasubramani@cgiar.org"
__version__ = "0.0.1"

import argparse
import sys
from collections import defaultdict

# CONSTANTS
NEWLINE = '\n'
CSV = ','
TAB = "\t"
SPACE = " "
IUPAC = ("A", "T", "G", "C")
tasks = ['GETSTATS', 'LGC2HMP', 'LGC2MTX', 'LGC2FLAPJACK', 'F1CHECK', 'CONSENSUS', 'RENAMESAMPLE']
HMPHEAD = "rs#" + TAB \
          + "alleles" + TAB \
          + "chrom" + TAB \
          + "pos" + TAB \
          + "strand" + TAB \
          + "assembly#" + TAB \
          + "center" + TAB \
          + "protLSID" + TAB \
          + "assayLSID" + TAB \
          + "panelLSID" + TAB \
          + "QCcode"
FJGHEAD = "# fjFile = GENOTYPE"
FJMHEAD = "# fjFile = MAP"
F1SUMMARYHEAD = "F1" + CSV \
                + "ParentA" + CSV \
                + "ParentB" + CSV \
                + "CountPolymorphic" + CSV \
                + "CountHet" + CSV \
                + "CountHomoPA" + CSV \
                + "CountHomoPB" + CSV \
                + "CountN" + CSV \
                + "PercentageHet" + CSV \
                + "Comment" + CSV \
                + "Marker [CallA/CallB/CallF1]"

CONSENSUSSUMHEAD = "GenotypeID" + CSV \
                   + "MarkerID" + CSV \
                   + "ConsensusBase" + CSV \
                   + "ConsensusScore" + CSV \
                   + "A" + CSV \
                   + "C" + CSV \
                   + "G" + CSV \
                   + "K" + CSV \
                   + "M" + CSV \
                   + "N" + CSV \
                   + "S" + CSV \
                   + "R" + CSV \
                   + "T" + CSV \
                   + "W" + CSV \
                   + "Y"

parser = argparse.ArgumentParser(description="QC pipeline: Processes the LGC file for tasks involved in purity check")
parser.add_argument("--task", dest="task", type=str,
                    metavar="<STRING>",
                    help="Task to be performed. Tasks are: 'LGC2HMP', 'LGC2FLAPJACK', 'LGC2MTX', 'F1CHECK', 'CONSENSUS', 'RENAMESAMPLE'")
parser.add_argument("--snp-info", dest="snpInfoFile",
                    metavar="<FILENAME>", help="SNP information file.")
parser.add_argument("--grid-file", dest="gridFile",
                    metavar="<FILENAME>", help="LGC Grid Matrix File")
parser.add_argument("--parent-info", dest="parentInfoFile",
                    metavar="<FILENAME>", help="Parent information file")
parser.add_argument("--sample-info", dest="sampleInfoFile",
                    metavar="<FILENAME>", help="Sample information file")
parser.add_argument("--out", dest="outPrefix",
                    metavar="<STRING>", help="Output filename prefix")
parser.add_argument("--cutOff", dest="cutOff",
                    metavar="<Integer>", help="Percentage expected heterozygosity for F1 verification")

# Parse commandline arguments
options = parser.parse_args()


def tree():
    """
    Initialize nested dictionary
    :return:
    """
    return defaultdict(tree)


def getConsensusBase(baseCount):
    i = 0
    count = []
    allele = []
    baseDict = init_bases()
    for base in sorted(baseCount, key=baseCount.get, reverse=True):
        if base != 'N':
            count.append(baseCount[base])
            allele.append(base)
            i += 1
        # if i == 1 and base != 'N':
        #     return base
    if count[0] == 0:
        return 'N', 0
    if count[0] == count[1]:
        tempBase = allele[0] + allele[1]
        if tempBase in baseDict:
            base = baseDict[tempBase]
        else:
            if (allele[0] not in IUPAC) and (allele[1] not in IUPAC):
                tempBase = 'N'
            elif allele[0] not in IUPAC:
                tempBase = allele[1]
            elif allele[1] not in IUPAC:
                tempBase = allele[0]
            base = baseDict[tempBase]
        purity = 50
    else:
        base = allele[0]
        purity = int((count[0] * 100) / sumOfList(count, len(count)))

    return base, purity


def sumOfList(listName, listSize):
    if (listSize == 0):
        return 0
    else:
        return listName[listSize - 1] + sumOfList(listName, listSize - 1)


def initBaseCount():
    baseCount = tree()
    baseCount['A'] = 0
    baseCount['T'] = 0
    baseCount['G'] = 0
    baseCount['C'] = 0
    baseCount['N'] = 0
    baseCount['R'] = 0
    baseCount['S'] = 0
    baseCount['Y'] = 0
    baseCount['K'] = 0
    baseCount['W'] = 0
    baseCount['K'] = 0
    baseCount['M'] = 0
    return baseCount


def init_bases():
    """
    Initialize a dictionary of 2 letter -> 1 letter IUPAC codes.
    :return:
    """
    bases = defaultdict(dict)
    bases['AA'] = 'A'
    bases['TT'] = 'T'
    bases['GG'] = 'G'
    bases['CC'] = 'C'
    bases['NN'] = 'N'
    bases['AT'] = 'W'
    bases['TA'] = 'W'
    bases['AG'] = 'R'
    bases['GA'] = 'R'
    bases['AC'] = 'M'
    bases['CA'] = 'M'
    bases['TG'] = 'K'
    bases['GT'] = 'K'
    bases['TC'] = 'Y'
    bases['CT'] = 'Y'
    bases['GC'] = 'S'
    bases['CG'] = 'S'
    bases['Uncallable'] = 'N'
    bases['Unused'] = 'N'
    bases['?'] = 'N'
    bases['NTC'] = 'N'
    bases['A'] = 'A'
    bases['T'] = 'T'
    bases['G'] = 'G'
    bases['C'] = 'C'
    bases['N'] = 'N'
    bases['S'] = 'S'
    bases['R'] = 'R'
    bases['M'] = 'M'
    bases['Y'] = 'Y'
    bases['K'] = 'K'
    bases['W'] = 'W'

    return bases


def init_OneLetterDict():
    bases = tree()
    bases['A'] = 'AA'
    bases['T'] = 'TT'
    bases['G'] = 'GG'
    bases['C'] = 'CC'
    bases['N'] = 'NN'
    bases['S'] = 'GC'
    bases['R'] = 'AG'
    bases['M'] = 'AC'
    bases['Y'] = 'CT'
    bases['K'] = 'GT'
    bases['W'] = 'AT'

    return bases


def die(errorMessage):
    """
    Print error message and exit from process. DIRTY METHOD
    :param errorMessage:
    :return:
    """
    print(errorMessage)
    sys.exit()


def writeFJGT(outFile, gtDataGM, markerIDs):
    """
    Method to write Flapjack Genotype file
    :param outFile:
    :param gtDataGM:
    :param markerInfo:
    :return:
    """
    # Opening output file handle
    outFileHandle = open(outFile, 'w')
    bases = init_OneLetterDict()
    # Build and write the first/header line for Flapjack Genotype
    for gtID in sorted(gtDataGM):
        outFileHandle.write(FJGHEAD + NEWLINE)
        for i in range(1, len(markerIDs)):
            outFileHandle.write(TAB + markerIDs[i])
        outFileHandle.write(NEWLINE)
        break
    # Build and write marker lines for Flapjack Genotype
    for gtID in gtDataGM:
        for gtIDX in sorted(gtDataGM[gtID]):
            outFileHandle.write(gtID)
            for i in range(1, len(markerIDs)):
                outBase = str(gtDataGM[gtID][gtIDX][markerIDs[i]])
                if outBase in bases:
                    outBaseCall = bases[outBase]
                else:
                    outBaseCall = '--'
                #
                outFileHandle.write(TAB + str(outBaseCall))
            outFileHandle.write(NEWLINE)


def writeGridFile(outFile, gtDataGM, markerIDs):
    """
    Method to write Grid Genotype file
    :param outFile:
    :param gtDataGM:
    :param markerInfo:
    :return:
    """
    # Opening output file handle
    outFileHandle = open(outFile, 'w')
    bases = init_OneLetterDict()
    # Build and write the first/header line for Flapjack Genotype
    for gtID in sorted(gtDataGM):
        outFileHandle.write("DNA \\ Assay")
        for i in range(1, len(markerIDs)):
            outFileHandle.write(CSV + markerIDs[i])
        outFileHandle.write(NEWLINE)
        break
    # Build and write marker lines for Flapjack Genotype
    for gtID in gtDataGM:
        for gtIDX in sorted(gtDataGM[gtID]):
            outFileHandle.write(gtID)
            for i in range(1, len(markerIDs)):
                outBase = str(gtDataGM[gtID][gtIDX][markerIDs[i]])
                if outBase in bases:
                    outBaseCall = bases[outBase]
                    outBaseCall = outBaseCall[0] + ":" + outBaseCall[1]
                else:
                    outBaseCall = '-:-'
                #
                outFileHandle.write(CSV + str(outBaseCall))
            outFileHandle.write(NEWLINE)


def writeFJKMP(outFile, markerInfo, markerIDs):
    """
    Method to write Flapjack Map file
    :param outFile:
    :param markerInfo:
    :return:
    """
    # Opening output file handle
    outFileHandle = open(outFile, 'w')
    # Build and write the first/header line for Flapjack Map
    outFileHandle.write(FJMHEAD + NEWLINE)
    # Build and write marker lines for Flapjack Map
    for i in range(1, len(markerIDs)):
        outFileHandle.write(markerIDs[i])
        outFileHandle.write(TAB + markerInfo[markerIDs[i]]['chr'])
        outFileHandle.write(TAB + toMB(markerInfo[markerIDs[i]]['pos']))
        outFileHandle.write(NEWLINE)


def toMB(position):
    """
    Method to convert physical position to MBs
    :param position:
    :return:
    """
    return str(float(position) / 1000000)


def writeHMP(outFile, gtDataMG, markerInfo, markerIDs):
    """
    From LGCdata dictionary and Marker Info dictionary create Hapmap file
    :param outFile: Output Hapmap file name
    :param gtDataMG: LGC data dictionary. marker -> Genotype -> AlleleCall
    :param markerInfo: Marker Info Dictionary. marker -> chr/pos
    :return:
    """

    # Opening output file handle
    outFileHandle = open(outFile, 'w')
    # Build and write the first/header line for Hapmap
    for i in range(1, len(markerIDs)):
        outFileHandle.write(HMPHEAD)
        for gtID in sorted(gtDataMG[markerIDs[i]]):
            for gtIDX in sorted(gtDataMG[markerIDs[i]][gtID]):
                outFileHandle.write(TAB + gtID)
        outFileHandle.write(NEWLINE)
        break

    # Build and write marker lines for Hapmap
    for i in range(1, len(markerIDs)):

        if markerIDs[i] in markerInfo:
            chr = markerInfo[markerIDs[i]]['chr']
            pos = markerInfo[markerIDs[i]]['pos']
        else:
            chr = "NA"
            pos = "NA"

        outFileHandle.write(markerIDs[i] + TAB \
                            + "N/N" + TAB \
                            + chr + TAB \
                            + pos + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" \
                            )
        for gtID in sorted(gtDataMG[markerIDs[i]]):
            for gtIDX in sorted(gtDataMG[markerIDs[i]][gtID]):
                outFileHandle.write(TAB + gtDataMG[markerIDs[i]][gtID][gtIDX])
        outFileHandle.write(NEWLINE)


def getString(gtDataGM, parent):
    """
    Get string from dictionary
    :param gtDataGM: dictionary[genotype][marker] = baseCall
    :param parent: GenotypeName of Parent
    :return: String of base calls
    """
    parentString = ""
    for marker in sorted(gtDataGM[parent]):
        parentString = parentString + gtDataGM[parent][marker]
    return parentString


def checkParents(outF1resultsFileHandle, gtDataMG, parentA, parentB, fOne, cutOff):
    """
    Get list of polymorphic markers between the parents based on LGC data

    :param gtDataMG: dictionary[marker][genotype] = baseCall
    :param parentA: GenotypeName of Parent A
    :param parentB: GenotypeName of Parent B
    :return: List of polymorphic markers
    """
    polyMorphicList = []
    countB = 0
    countA = 0
    countHet = 0
    countPM = 0
    countN = 0
    for marker in sorted(gtDataMG):
        bases = init_bases()
        callA = gtDataMG[marker][parentA][1]
        callB = gtDataMG[marker][parentB][1]
        callFone = gtDataMG[marker][fOne][1]
        if IUPAC.__contains__(callA) and IUPAC.__contains__(callB):
            if callA != callB:
                polyMorphicList.append(marker + "[" + callA + "/" + callB + "/" + callFone + "]")
                countPM += 1
                if bases[callA + callB] == callFone:
                    countHet += 1
                if callA == callFone:
                    countA += 1
                if callB == callFone:
                    countB += 1
                if callFone == "N":
                    countN += 1
    # return polyMorphicList, str(countPM), str(countHet), str(countA), str(countB)
    if (countPM - countN) > 0:
        expHet = int((countHet * 100) / (countPM - countN))
    else:
        expHet = 0
    if expHet == 100:
        trueF1 = "SuccessfulF1"
    elif expHet > int(cutOff):
        trueF1 = "Inconclusive"
    else:
        trueF1 = "UnsuccessfulF1"
    outF1resultsFileHandle.write(
        fOne + CSV + parentA + CSV + parentB + CSV + str(countPM) + CSV + str(countHet) + CSV + str(countA) + CSV + str(
            countB) + CSV + str(countN) + CSV + str(expHet) + CSV + trueF1 + CSV + SPACE.join(
            polyMorphicList) + NEWLINE)


def getMarkerInfoDict(snpInfoFile):
    """

    :param snpInfoFile:
    :return: disctionaty of [MarkerID][chr/pos]
    """
    with open(snpInfoFile) as snpInfoHandle:
        markerInfo = defaultdict(dict)
        for line in snpInfoHandle:
            line = line.strip()
            lineEntries = line.split(TAB)
            if lineEntries[2].isalpha():
                continue
            markerInfo[lineEntries[0]]['chr'] = lineEntries[1]
            markerInfo[lineEntries[0]]['pos'] = lineEntries[2]
    return markerInfo


def readLGCgridFile(gridFile):
    """

    :param gridFile:
    :return: dictionary of [Marker][Genotype], [Genotype][Marker] and List of [MarkerIDs]
    """
    # Initialize 2 letter nucleotide dictionary
    bases = init_bases()
    # Reading LGC Grid file. Delimiter = Comma
    with open(gridFile) as gridHandle:
        lineNo = 0
        dataFlag = 0
        gtDataMG = tree()
        gtDataGM = tree()
        gtCount = tree()
        for line in gridHandle:
            line = line.rstrip()
            lineEntries = line.split(CSV)
            lineNo += 1
            if lineEntries[0] == "DNA \\ Assay":
                dataFlag = 1
                markerIDs = lineEntries
                continue
            if line == "":
                continue
            if dataFlag == 0:
                continue
            if lineEntries[0] == "NTC":
                continue
            if lineEntries[0] in gtCount:
                gtCount[lineEntries[0]] = str(int(gtCount[lineEntries[0]]) + 1)
            else:
                gtCount[lineEntries[0]] = "1"
            gtIdx = gtCount[lineEntries[0]]
            for i in range(1, len(lineEntries)):
                if not lineEntries[i] == "":
                    gtDataMG[markerIDs[i]][lineEntries[0]][gtIdx] = bases[lineEntries[i].replace(":", "")]
                    gtDataGM[lineEntries[0]][gtIdx][markerIDs[i]] = bases[lineEntries[i].replace(":", "")]
                else:
                    gtDataMG[markerIDs[i]][lineEntries[0]][gtIdx] = bases['NN']
                    gtDataGM[lineEntries[0]][gtIdx][markerIDs[i]] = bases['NN']
    return gtDataMG, gtDataGM, markerIDs


def readFjGTfile(flapjackGtFile):
    """

        :param flapjackGtFile:
        :return: dictionary of [Marker][Genotype], [Genotype][Marker] and List of [MarkerIDs]
        """
    # Initialize 2 letter nucleotide dictionary
    bases = init_bases()
    # Reading Flapjack Genotype File file. Delimiter = Comma
    with open(flapjackGtFile) as gridHandle:
        lineNo = 0
        dataFlag = 0
        gtDataMG = tree()
        gtDataGM = tree()
        gtCount = tree()
        for line in gridHandle:
            line = line.rstrip()
            lineEntries = line.split(CSV)
            lineNo += 1
            if line == FJGHEAD:
                dataFlag = 1
                continue
            if line == "":
                continue
            if dataFlag == 0:
                continue
            if dataFlag == 1:
                markerIDs = lineEntries
                dataFlag = 2
                continue
            if lineEntries[0] == "NTC":
                continue
            if lineEntries[0] in gtCount:
                gtCount[lineEntries[0]] = str(int(gtCount[lineEntries[0]]) + 1)
            else:
                gtCount[lineEntries[0]] = "1"
            gtIdx = gtCount[lineEntries[0]]
            for i in range(1, len(lineEntries)):
                if not lineEntries[i] == "":
                    gtDataMG[markerIDs[i]][lineEntries[0]][gtIdx] = bases[lineEntries[i].replace(":", "")]
                    gtDataGM[lineEntries[0]][gtIdx][markerIDs[i]] = bases[lineEntries[i].replace(":", "")]
                else:
                    gtDataMG[markerIDs[i]][lineEntries[0]][gtIdx] = bases['NN']
                    gtDataGM[lineEntries[0]][gtIdx][markerIDs[i]] = bases['NN']
    return gtDataMG, gtDataGM, markerIDs


def processReplicates(gtDataMG, outConsensusSummaryFile):
    # global consensusMG, consensusGM, bases, mID, gID, baseCount, idx, base, baseList
    consensusMG = tree()
    consensusGM = tree()
    bases = init_OneLetterDict()
    basesTwoToOne = init_bases()
    outConsensusSummaryFileHandle = open(outConsensusSummaryFile, 'w')
    outConsensusSummaryFileHandle.write(CONSENSUSSUMHEAD + NEWLINE)
    for mID in gtDataMG:
        for gID in gtDataMG[mID]:
            baseCount = initBaseCount()
            for idx in sorted(gtDataMG[mID][gID]):
                base = str(gtDataMG[mID][gID][idx])
                # baseList = bases[base]
                # baseCount[baseList[0]] += 1
                # baseCount[baseList[1]] += 1
                baseCount[base] += 1
            (base, purity) = getConsensusBase(baseCount)
            if len(gtDataMG[mID][gID]) > 1:
                # print(mID + TAB + gID + str(base))
                outConsensusSummaryFileHandle.write(gID + CSV + mID + CSV + str(base) + CSV + str(purity) + CSV
                                                    + CSV.join(map(str, list(baseCount.values()))) + NEWLINE)

            consensusGM[gID][1][mID] = base
            consensusMG[mID][gID][1] = base
    outConsensusSummaryFileHandle.close()
    return consensusGM, consensusMG


# Here we go
"""__main__"""

# Task to be performed
task = options.task

# Return help message if the task is not provided or given wrong
if task not in tasks:
    print("Invalid/No Task provided.")
    parser.print_help()


if task == "GETSTATS":
    required = "gridFile"
    if options.__dict__[required] is None:
        print("Option/Value Missing for", required)
        parser.print_help()
        sys.exit(1)
    
    gridFile = options.gridFile
    gtDataMG, gtDataGM, markerIDs = readLGCgridFile(gridFile)
    print("Number of Markers: ", len(markerIDs))
    print("Number of samples:", len(gtDataGM.keys()))

# To convert LGC Grid file to Hapmap file
if task == "LGC2HMP":
    # Look for the mandatory options are parsed/provided.
    required = "gridFile snpInfoFile outPrefix".split(" ")
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)

    gridFile = options.gridFile
    snpInfoFile = options.snpInfoFile
    parentInfoFile = options.parentInfoFile
    outHmpFile = options.outPrefix + ".hmp.txt"

    gtDataMG, gtDataGM, markerIDs = readLGCgridFile(gridFile)
    markerInfo = getMarkerInfoDict(snpInfoFile)

    # Calling the method to write output Hapmap file
    writeHMP(outHmpFile, gtDataMG, markerInfo, markerIDs)

# To convert LGC Grid file to Flapjack file(s)
if task == "LGC2FLAPJACK":
    # Look for the mandatory options are parsed/provided.
    required = "gridFile outPrefix".split(" ")
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)

    gridFile = options.gridFile
    snpInfoFile = options.snpInfoFile
    outFjkGTFile = options.outPrefix + ".genotype"
    outFjkMPFile = options.outPrefix + ".map"

    gtDataMG, gtDataGM, markerIDs = readLGCgridFile(gridFile)
    writeFJGT(outFjkGTFile, gtDataGM, markerIDs)

    if snpInfoFile is not None:
        # Reading Snp information file. Delimiter = TAB. rsID<TAB>chr<TAB>pos
        markerInfo = getMarkerInfoDict(snpInfoFile)
        writeFJKMP(outFjkMPFile, markerInfo, markerIDs)

if task == "CONSENSUS":
    # Look for the mandatory options are parsed/provided.
    required = "gridFile outPrefix".split(" ")
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)
    gridFile = options.gridFile
    outFjkGTFile = options.outPrefix + ".genotype"
    outConsensusSummaryFile = options.outPrefix + "_ConsensusSummary.csv"
    outGridFile = options.outPrefix + "_Grid.csv"
    gtDataMG, gtDataGM, markerIDs = readLGCgridFile(gridFile)
    consensusGM, consensusMG = processReplicates(gtDataMG, outConsensusSummaryFile)
    writeFJGT(outFjkGTFile, consensusGM, markerIDs)
    writeGridFile(outGridFile, consensusGM, markerIDs)
    if options.snpInfoFile is not None:
        snpInfoFile = options.snpInfoFile
        outHmpFile = options.outPrefix + ".hmp.txt"
        markerInfo = getMarkerInfoDict(snpInfoFile)
        writeHMP(outHmpFile, consensusMG, markerInfo, markerIDs)

if task == "F1CHECK":
    # Look for the mandatory options are parsed/provided.
    required = "gridFile parentInfoFile outPrefix".split(" ")
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)

    gridFile = options.gridFile
    parentInfoFile = options.parentInfoFile
    outFjkGTFile = options.outPrefix + ".genotype"
    outF1resultsFile = options.outPrefix + "_F1summary.csv"
    outConsensusSummaryFile = options.outPrefix + "_ConsensusSummary.csv"
    cutOff = options.cutOff
    if cutOff is None:
        cutOff = 100
    gtDataMG, gtDataGM, markerIDs = readLGCgridFile(gridFile)
    consensusGM, consensusMG = processReplicates(gtDataMG, outConsensusSummaryFile)
    writeFJGT(outFjkGTFile, consensusGM, markerIDs)

    # Opening Parent Info file and making a dictionary of dict[F1] = listOf(pA, pB)
    with open(parentInfoFile) as parentInfoHandle:
        fOneInfo = defaultdict(dict)
        for line in parentInfoHandle:
            line = line.rstrip()
            lineEntries = line.split(TAB)
            if len(lineEntries) < 3:
                continue
            if lineEntries[1] == "--" or lineEntries[1] is None:
                continue
            fOneInfo[lineEntries[0]] = [lineEntries[1], lineEntries[2]]

    outFileHandle = open(outF1resultsFile, 'w')
    outFileHandle.write(F1SUMMARYHEAD + NEWLINE)
    for fOne in fOneInfo:
        parentA = fOneInfo[fOne][0]
        parentB = fOneInfo[fOne][1]
        # print("Processing F1: " + fOne)
        if not consensusGM[fOne]:
            print(fOne + "=" + parentA + "x" + parentB + TAB + "F1: " + fOne + " Data Missing")
            continue
        if not consensusGM[parentA]:
            print(fOne + "=" + parentA + "x" + parentB + TAB + "Parent: " + parentA + " Data Missing")
            continue
        if not consensusGM[parentB]:
            print(fOne + "=" + parentA + "x" + parentB + TAB + "Parent: " + parentB + " Data Missing")
            continue
        if consensusGM[parentA] and consensusGM[parentB]:
            # Get list of polymorphic markers between the parents based on LGC data
            checkParents(outFileHandle, consensusMG, parentA, parentB, fOne, cutOff)
        else:
            print("Skipping F1:" + fOne)
    outFileHandle.close()

if task == "RENAMESAMPLE":
    # Look for the mandatory options are parsed/provided.
    required = "gridFile sampleInfoFile outPrefix".split(" ")
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)

    gridFile = options.gridFile
    sampleInfoFile = options.sampleInfoFile
    outGridFile = options.outPrefix + "_Grid.csv"

    flag = 0

    with open(sampleInfoFile) as sampleInfoFileHandle:
        sampleInfo = tree()
        for line in sampleInfoFileHandle:
            line = line.rstrip()
            lineEntries = line.split(TAB)
            if lineEntries[0] == "" or lineEntries[1] == "":
                continue
            sampleInfo[lineEntries[0]] = lineEntries[1]
    outGridFileHandle = open(outGridFile, 'w')

    with open(gridFile) as gridFileHandle:
        for line in gridFileHandle:
            line = line.rstrip()
            lineEntries = line.split(CSV)
            if flag == 0:
                outGridFileHandle.write(line + NEWLINE)
                if lineEntries[0] == "DNA \ Assay":
                    flag = 1
                continue
            if lineEntries[0] in sampleInfo:
                lineEntries[0] = sampleInfo[lineEntries[0]]
            outGridFileHandle.write(CSV.join(lineEntries) + NEWLINE)
