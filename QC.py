#!/usr/bin/env python
import string
import sys

__author__ = "Sivasubramani S"
__email__ = "s.sivasubramani@cgiar.org"
__version__ = "0.0.1"

import argparse
from collections import defaultdict

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


options = parser.parse_args()

def init_bases():
    """

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
    print(errorMessage)
    sys.exit()

def writeHMP(outFile, lgcDataMG, markerInfo):
    outFileHandle = open(outFile, 'w')
    for markerID in sorted(lgcDataMG):
        outFileHandle.write(HMPHEAD)
        for gtID in sorted(lgcDataMG[markerID]):
            outFileHandle.write(TAB + gtID)
        outFileHandle.write(NEWLINE)
        break

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
                            + "NA" + TAB \
                            )
        for gtID in sorted(lgcDataMG[markerID]):
            outFileHandle.write("\t" + lgcDataMG[markerID][gtID])
        outFileHandle.write(NEWLINE)

"""__main__"""

task = options.task

if not task in tasks:
    print("Invalid/No Task provided.")
    parser.print_help()

if task == "LGC2HMP":
    required = "gridFile snpInfoFile outHmpFile".split(" ")
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)

    gridFile = options.gridFile
    snpInfoFile = options.snpInfoFile
    outHmpFile = options.outHmpFile

    bases = init_bases()
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
                lgcDataMG[markerIDs[i]][lineEntries[0]] = bases[lineEntries[i].replace(":","")]
                lgcDataGM[lineEntries[0]][markerIDs[i]] = bases[lineEntries[i].replace(":", "")]


    with open(snpInfoFile) as snpInfoHandle:
        markerInfo = defaultdict(dict)
        for line in snpInfoHandle:
            line = line.strip()
            lineEntries = line.split(TAB)
            markerInfo[lineEntries[0]]['chr'] = lineEntries[1]
            markerInfo[lineEntries[0]]['pos'] = lineEntries[2]

    writeHMP(outHmpFile, lgcDataMG, markerInfo)





