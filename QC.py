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
IUPAC = ("A", "T", "G", "C")
tasks = ['LGC2HMP', 'LGC2MTX', 'F1CHECK']
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

parser = argparse.ArgumentParser(description="QC pipeline: Processes the LGC file for tasks involved in purity check")
parser.add_argument("--task", dest="task", type=str,
                  metavar="<string>", help="Task to be performed. Eg: LGC2HMP")
parser.add_argument("--snp-info", dest="snpInfoFile",
                  metavar="<FILENAME>", help="SNP information file.")
parser.add_argument("--grid-file", dest="gridFile",
                  metavar="<FILENAME>", help="LGC Grid Matrix File")
parser.add_argument("--parent-info", dest="parentInfoFile",
                  metavar="<FILENAME>", help="Parent information file")
parser.add_argument("--out-HMP", dest="outHmpFile",
                  metavar="<FILENAME>", help="Output Hapmap File")

# Parse commandline arguments
options = parser.parse_args()

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
    bases['?'] = 'N'
    bases['NTC'] = 'N'
    return bases

def die(errorMessage):
    """
    Print error message and exit from process. DIRTY METHOD
    :param errorMessage:
    :return:
    """
    print(errorMessage)
    sys.exit()

def writeHMP(outFile, lgcDataMG, markerInfo):
    """
    From LGCdata dictionary and Marker Info dictionary create Hapmap file
    :param outFile: Output Hapmap file name
    :param lgcDataMG: LGC data dictionary. marker -> Genotype -> AlleleCall
    :param markerInfo: Marker Info Dictionary. marker -> chr/pos
    :return:
    """

    # Opening output file handle
    outFileHandle = open(outFile, 'w')
    # Build and write the first/header line for Hapmap
    for markerID in sorted(lgcDataMG):
        outFileHandle.write(HMPHEAD)
        for gtID in sorted(lgcDataMG[markerID]):
            outFileHandle.write(TAB + gtID)
        outFileHandle.write(NEWLINE)
        break

    # Build and write marker lines for Hapmap
    for markerID in sorted(lgcDataMG):
        outFileHandle.write(markerID + TAB \
                            + "N/N" + TAB \
                            + markerInfo[markerID]['chr'] + TAB \
                            + markerInfo[markerID]['pos'] + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" + TAB \
                            + "NA" \
                            )
        for gtID in sorted(lgcDataMG[markerID]):
            outFileHandle.write("\t" + lgcDataMG[markerID][gtID])

        outFileHandle.write(NEWLINE)

def getString(lgcDataGM, parent):
    """
    Get string from dictionary
    :param lgcDataGM: dictionary[genotype][marker] = baseCall
    :param parent: GenotypeName of Parent
    :return: String of base calls
    """
    parentString = ""
    for marker in sorted(lgcDataGM[parent]):
        parentString = parentString + lgcDataGM[parent][marker]
    return parentString

def checkParents(lgcDataMG, parentA, parentB):
    """
    Get list of polymorphic markers between the parents based on LGC data

    :param lgcDataMG: dictionary[marker][genotype] = baseCall
    :param parentA: GenotypeName of Parent A
    :param parentB: GenotypeName of Parent B
    :return: List of polymorphic markers
    """
    polyMorphicList = []
    for marker in sorted(lgcDataMG):
        callA = lgcDataMG[marker][parentA]
        callB = lgcDataMG[marker][parentB]
        if IUPAC.__contains__(callA) and IUPAC.__contains__(callB):
            if  callA != callB :
                polyMorphicList.append(marker + "[" + callA + "/" + callB + "]")
    return polyMorphicList


# Here we go
"""__main__"""

# Task to be performed
task = options.task

# Return help message if the task is not provided or given wrong
if not task in tasks:
    print("Invalid/No Task provided.")
    parser.print_help()

# To convert LGC Grid file to Hapmap file

if task == "LGC2HMP":
    # Look for the mandatory options are parsed/provided.
    required = "gridFile snpInfoFile parentInfoFile outHmpFile".split(" ")
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)

    gridFile = options.gridFile
    snpInfoFile = options.snpInfoFile
    parentInfoFile = options.parentInfoFile
    outHmpFile = options.outHmpFile

    # Initialize 2 letter nucleotide dictionary
    bases = init_bases()
    # Reading LGC Grid file. Delimiter = Comma
    with open(gridFile) as gridHandle:
        lineNo = 0
        dataFlag = 0
        lgcDataMG = defaultdict(dict)
        lgcDataGM = defaultdict(dict)
        for line in gridHandle:
            line = line.rstrip()
            lineEntries = line.split(CSV)
            lineNo += 1
            if lineEntries[0] == "DNA \\ Assay":
                dataFlag = 1
                markerIDs = lineEntries
                continue
            if dataFlag == 0:
                continue
            if lineEntries[0] == "NTC":
                continue
            for i in range(1, len(lineEntries)):
                if not lineEntries[i] == "" :
                    lgcDataMG[markerIDs[i]][lineEntries[0]] = bases[lineEntries[i].replace(":","")]
                    lgcDataGM[lineEntries[0]][markerIDs[i]] = bases[lineEntries[i].replace(":", "")]
                else :
                    lgcDataMG[markerIDs[i]][lineEntries[0]] = bases['NN']
                    lgcDataGM[lineEntries[0]][markerIDs[i]] = bases['NN']
    # Reading Snp information file. Delimiter = TAB. rsID<TAB>chr<TAB>pos
    with open(snpInfoFile) as snpInfoHandle:
        markerInfo = defaultdict(dict)
        for line in snpInfoHandle:
            line = line.strip()
            lineEntries = line.split(TAB)
            markerInfo[lineEntries[0]]['chr'] = lineEntries[1]
            markerInfo[lineEntries[0]]['pos'] = lineEntries[2]

    # Calling the method to write output Hapmap file
    writeHMP(outHmpFile, lgcDataMG, markerInfo)

    # Opening Parent Info file and making a dictionary of dict[F1] = listOf(pA, pB)
    with open(parentInfoFile) as parentInfoHandle:
        fOneInfo = defaultdict(dict)
        for line in parentInfoHandle:
            line = line.rstrip()
            lineEntries = line.split(TAB)
            fOneInfo[lineEntries[0]] = [lineEntries[1],lineEntries[2]]

    polyMorphic = defaultdict(dict)
    for fOne in fOneInfo:
        parentA = fOneInfo[fOne][0]
        parentB = fOneInfo[fOne][1]
        if not lgcDataGM[parentA]:
            print("Parent: " + parentA + " Data Missing")
        if not lgcDataGM[parentB]:
            print("Parent: " + parentB + " Data Missing")
        if lgcDataGM[parentA] and lgcDataGM[parentB]:
            # Get list of polymorphic markers between the parents based on LGC data
            polyMorphic[parentA][parentB] = checkParents(lgcDataMG, parentA, parentB)
            print(parentA + TAB + parentB + TAB + str(len(polyMorphic[parentA][parentB])) + TAB + CSV.join(polyMorphic[parentA][parentB]))
        else:
            print("Skipping F1:" + fOne)





