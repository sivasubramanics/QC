#!/usr/bin/env python3

__author__ = "Sivasubramani S"
__email__ = "s.sivasubramani@cgiar.org"
__version__ = "0.2"


# Add marker quality
# Convert inconclusive to high-possible F1

# from QC import TAB
import os, sys, argparse
from collections import defaultdict
# from LGCutils import *

tasks = ['LGC2GRID', 'GETSTATS', 'LGC2FJK', 'F1CHECK', 'CONSENSUS', 'RENAME', 'MERGEDATA', 'FWDBREED']
MAPSEP = "\t"
TAB = "\t"
NEWLINE = "\n"
SPACE = " "
CSV = ','
GRIDSEP = ","
FJGHEAD = "# fjFile = GENOTYPE"
F1SUMMARYHEAD = "F1" + GRIDSEP \
                + "ParentA" + GRIDSEP \
                + "ParentB" + GRIDSEP \
                + "TotalMarkers" + GRIDSEP \
                + "CountPolymorphic" + GRIDSEP \
                + "CountHet" + GRIDSEP \
                + "CountHomoPA" + GRIDSEP \
                + "CountHomoPB" + GRIDSEP \
                + "CountN" + GRIDSEP \
                + "PercentageHet" + GRIDSEP \
                + "Comment" + GRIDSEP \
                + "ParentAvaialable" + GRIDSEP \
                + "Marker [CallA/CallB/CallF1]"
IUPAC = ("A", "T", "G", "C")
MISSING_CALLS = ('N:N', '?:?', "?", "Uncallable", "Unused", "missing", "Empty", "NTC")
sampleMap = defaultdict(list)
parentMap = defaultdict(list)
gtData = defaultdict(list)
consensusData = defaultdict(list)
consensusOutData = defaultdict(list)
consensusCutOff = 0.5
selectedMarkers = []

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
    bases['Empty'] = 'N'
    bases['?'] = 'N'
    bases['NTC'] = 'N'
    bases['A'] = 'AA'
    bases['T'] = 'TT'
    bases['G'] = 'GG'
    bases['C'] = 'CC'
    bases['N'] = 'NN'
    bases['S'] = 'CG'
    bases['R'] = 'AG'
    bases['M'] = 'AC'
    bases['Y'] = 'CT'
    bases['K'] = 'GT'
    bases['W'] = 'AT'
    return bases

def init_baseCount():
    baseCount = defaultdict()
    baseCount['A'] = 0
    baseCount['T'] = 0
    baseCount['G'] = 0
    baseCount['C'] = 0
    baseCount['N'] = 0
    baseCount['R'] = 0
    baseCount['S'] = 0
    baseCount['M'] = 0
    baseCount['K'] = 0
    baseCount['Y'] = 0
    baseCount['W'] = 0
    baseCount['A'] = 0
    return baseCount

def addColon(seq):
    tmpSeq = list(seq)
    return ":".join(tmpSeq)

def addSlash(seq):
    # tmpSeq = list(seq)
    return seq.replace(":", "/")
    
def writeGridFile(outFile, gtData, markerIDs):
    """
    Method to write Grid Genotype file
    :param outFile:
    :param gtDataGM:
    :param markerInfo:
    :return:
    """
    # Opening output file handle
    outFileHandle = open(outFile, 'w')
    bases = init_bases()
    # Build and write the first/header line for Flapjack Genotype
    outFileHandle.write("DNA \\ Assay" + GRIDSEP + GRIDSEP.join(markerIDs))
    outFileHandle.write(NEWLINE)
    for sampleID in gtData:
        outLine = []
        outLine.append(sampleID)
        for baseCall in gtData[sampleID]:
            outLine.append(baseCall)
        outFileHandle.write(GRIDSEP.join(outLine))
        outFileHandle.write(NEWLINE)

def writeFjGTFile(outFile, gtData, markerIDs):
    # Opening output file handle
    outFileHandle = open(outFile, 'w')
    bases = init_bases()
    # Build and write the first/header line for Flapjack Genotype
    outFileHandle.write(FJGHEAD + NEWLINE)
    outFileHandle.write(MAPSEP + MAPSEP.join(markerIDs) + NEWLINE)
    for sampleID in gtData:
        outLine = []
        outLine.append(sampleID)
        for baseCall in gtData[sampleID]:
            outLine.append(addSlash(baseCall))
        outFileHandle.write(MAPSEP.join(outLine))
        outFileHandle.write(NEWLINE)

def writePurityTable(outFile, parentMap, consensusOutData):
    outFileHandle = open(outFile, 'w')
    outFileHandle.write("SAMPLE_NAME" + GRIDSEP + "TotalMarkers" + GRIDSEP + "noMatch" + GRIDSEP + "noMismatch" + GRIDSEP + "noBothMissing" + GRIDSEP + "noConsensusMissing" + GRIDSEP + "noSampleMissing" + GRIDSEP + "PercentagePurity")
    outFileHandle.write(NEWLINE)
    for parentID in parentMap:
        for sampleID in parentMap[parentID]:
            if sampleID not in consensusOutData:
                continue
            # total, match, mismatch, bothN, conN, sampleN
            count = [0, 0, 0, 0, 0, 0, 0]
            for i in range(0, len(markerIDs)):
                count[0] += 1
                if consensusOutData[parentID][i] == consensusOutData[sampleID][i]:
                    if consensusOutData[parentID][i] == 'N:N':
                        count[3] += 1
                    else:
                        count[1] += 1
                else:
                    if consensusOutData[parentID][i] == 'N:N':
                        count[4] += 1
                    else:
                        if consensusOutData[sampleID][i] == 'N:N':
                            count[5] += 1
                        else:
                            count[2] += 1
            count[6] = round(count[1]/(count[0] - count[3]) * 100, 2)
            outFileHandle.write(sampleID + GRIDSEP + GRIDSEP.join(map(str, count)))
            outFileHandle.write(NEWLINE)

def readGridFile(inFile, sampleMap):
    bases = init_bases()
    with open(inFile) as gridHandle:
        for line in gridHandle:
            line = line.strip()
            lineEntries = line.split(GRIDSEP)
            if lineEntries[0] == "DNA \\ Assay":
                dataFlag = 1
                markerIDs = lineEntries[1:]
                continue
            if line == "":
                continue
            if dataFlag == 0:
                continue
            if lineEntries[0] == "NTC":
                continue
            if sampleMap[lineEntries[0]]:
                sampleId = sampleMap[lineEntries[0]][0]
            else:
                print(lineEntries[0],"does not present in sampleMap. Ignoring the sample")
                continue
            for i in range(1, len(lineEntries)):
                if not lineEntries[i] == "":
                    lineEntries[i] = bases[lineEntries[i].replace(":", "")]
                else:
                    lineEntries[i] = bases['NN']
            gtData[sampleId]=lineEntries[1:]
    return markerIDs, gtData

def readMapFile(inFile):
    mapDict = defaultdict(list)
    mapList = []
    lineNo = 0
    with open(inFile) as sm:
        for line in sm:
            lineNo += 1
            if lineNo == 1:
                continue
            line = line.strip()
            lineEntries = line.split(MAPSEP)
            if len(lineEntries) == 1:
                continue
            if lineEntries[1] == "--" or lineEntries[1] is None:
                continue
            mapDict[lineEntries[0]].append(lineEntries[1])
            if lineEntries[0] not in mapList:
                    mapList.append(lineEntries[0])
    print("LOG: Processed " + str(lineNo - 1) + " entries from file " + inFile + " .!")
    return mapDict, mapList

def next_key(tmpList, current_key):
    # temp = list(test_dict)
    try:
        res = tmpList[tmpList.index(current_key) + 1]
        if res == 'N':
            next_key(tmpList, res)
    except (ValueError, IndexError):
        res = None
    return res

def sort_dict(unsorted_dict):
    sorted_values = sorted(unsorted_dict.values(), reverse=True) # Sort the values
    sorted_dict = {}
    for i in sorted_values:
        for k in unsorted_dict.keys():
            if unsorted_dict[k] == i:
                sorted_dict[k] = unsorted_dict[k]
                break
    return sorted_dict

def getConsensusCalls(parentList, parentMap, markerIDs, gtData, outFile, consensusCutOff, task):
    outFileHandle = open(outFile, 'w')
    outFileHandle.write("parent" + GRIDSEP + "noReps" + GRIDSEP + "noMarkers" + GRIDSEP + "allNs" + GRIDSEP + "noMajorCalls" + GRIDSEP + "noMinorCalls")
    outFileHandle.write(NEWLINE)
    for parentID in parentMap:
        if task == "F1CHECK" or task == "FWDBREED":
            if parentID not in parentList:
                continue
        consensusBase = []
        count = [0,0,0]
        for markerNum in range(0, len(markerIDs)):
            # baseCount = init_baseCount()
            baseCount = defaultdict(list)
            countNs = 0
            for sampleID in parentMap[parentID]:
                if sampleID not in gtData:
                    continue
                baseCall = gtData[sampleID][markerNum]
                if baseCall == 'N:N':
                    countNs += 1
                else:
                    if baseCall in baseCount:
                        baseCount[baseCall] += 1
                    else:
                        baseCount[baseCall] = 1
            if countNs == len(parentMap[parentID]):
                consensusBase.append('N:N')
                count[0] += 1
            else:    
                # baseCount = sort_dict(baseCount)
                # print(baseCount)
                sortedBases = sorted(baseCount, key=baseCount.get, reverse=True)
                # print(sortedBases)
                # for base in sorted(baseCount, key=baseCount.get, reverse=True):
                for base in sortedBases:
                    if base == 'N:N':
                        continue
                    if baseCount[base]/(len(parentMap[parentID]) - countNs) > consensusCutOff:
                        consensusBase.append(base)
                        count[1] += 1
                    else:
                        if baseCount[base]/(len(parentMap[parentID]) - countNs) == consensusCutOff:
                            nextbase = next_key(sortedBases, base)
                            if baseCount[base] == baseCount[nextbase]:
                                consensusBase.append('N:N')
                                count[2] += 1
                            else:
                                consensusBase.append(base)
                                count[1] += 1
                        else:
                            consensusBase.append('N:N')
                            count[2] += 1
                    break
        consensusData[parentID] = consensusBase
        outFileHandle.write(parentID + GRIDSEP + str(len(parentMap[parentID])) + GRIDSEP + str(len(markerIDs)) + GRIDSEP + GRIDSEP.join(map(str, count)))
        outFileHandle.write(NEWLINE)
    return consensusData

def readParentInfo(inFile):
    with open(inFile) as parentInfoHandle:
        parentInfo = defaultdict(dict)
        lineNo = 0
        for line in parentInfoHandle:
            line = line.rstrip()
            lineNo += 1
            if lineNo == 1:
                continue
            lineEntries = line.split(MAPSEP)
            if len(lineEntries) < 3:
                continue
            if lineEntries[1] == "--" or lineEntries[1] is None:
                continue
            parentInfo[lineEntries[0]] = [lineEntries[1], lineEntries[2]]
        print("LOG: Processed " + str(lineNo - 1) + " entries from file " + inFile + " .!")
    return parentInfo

def getPolymorphicMarkers(outFileHandle, processedData, parentA, parentB, fOneName, markerIDs, cutOff):
    availability = "BOTH"
    if parentA in processedData and parentB in processedData:
        listA = processedData[parentA]
        listB = processedData[parentB]
        listX = processedData[fOneName]
        polyMorphicList = []
        bases = init_bases()
        # total(0), countPM(1), countHet(2), countA(3), countB(4), countN(5)
        count = [0, 0, 0, 0, 0, 0]
        for i in range(0, len(listX)):
            count[0] += 1
            if (listA[i] != listB[i]) and (listA[i] not in MISSING_CALLS) and (listB[i] not in MISSING_CALLS):
                parA = listA[i].split(":")
                parB = listB[i].split(":")
                fone = listX[i].split(":")
                if parA[0] == parA[1] and parB[0] == parB[1]:
                    polyMorphicList.append(markerIDs[i] + "[" + listA[i] + "/" + listB[i] + "/" + listX[i] + "]")
                    count[1] += 1
                    if listX[i] in MISSING_CALLS:
                        count[5] += 1
                    if fone[0] == fone[1]:
                        if fone[0] == parA[0]:
                            count[3] += 1
                        if fone[0] == parB[0]:
                            count[4] += 1
                    if fone[0] == parA[0] and fone[1] == parB[0]:
                        count[2] += 1
                    if fone[0] == parB[0] and fone[1] == parA[0]:
                        count[2] += 1
        if (count[1] - count[5]) > 0:
            expHet = int((count[2] * 100) / (count[1] - count[5]))
        else:
            expHet = 0
        if expHet == 100:
            trueF1 = "SuccessfulF1"
        elif expHet > int(cutOff):
            trueF1 = "Inconclusive"
        else:
            trueF1 = "UnsuccessfulF1"
    if parentA in processedData and parentB not in processedData:
        availability = "PARENTA"
        listA = processedData[parentA]
        # listB = processedData[parentB]
        listX = processedData[fOneName]
        polyMorphicList = []
        bases = init_bases()
        # total(0), countPM(1), countHet(2), countA(3), countB(4), countN(5)
        count = [0, 0, 0, 0, 0, 0]
        for i in range(0, len(listX)):
            count[0] += 1
            if listA[i] not in MISSING_CALLS and listX[i] not in MISSING_CALLS:
                polyMorphicList.append(markerIDs[i] + "[" + listA[i] + "/" + "-:-" + "/" + listX[i] + "]")
                parA = listA[i].split(":")
                fone = listX[i].split(":")
                if parA[0] == parA[1]:
                    count[1] += 1
                    if listX[i] in MISSING_CALLS:
                        count[5] += 1
                    elif listX[i] == listA[i]:
                        count[3] += 1
                    elif fone[0] != fone[1]:
                        count[2] += 1
                    else:
                        count[4] += 1
        if (count[1] - count[5]) > 0:
            expHet = int((count[2] * 100) / (count[1] - count[5]))
        else:
            expHet = 0
        if expHet == 100:
            trueF1 = "SuccessfulF1"
        elif expHet > int(cutOff):
            trueF1 = "Inconclusive"
        else:
            trueF1 = "UnsuccessfulF1"
    if parentA not in processedData and parentB in processedData:
        availability = "PARENTB"
        # listA = processedData[parentA]
        listB = processedData[parentB]
        listX = processedData[fOneName]
        polyMorphicList = []
        bases = init_bases()
        # total(0), countPM(1), countHet(2), countA(3), countB(4), countN(5)
        count = [0, 0, 0, 0, 0, 0]
        for i in range(0, len(listX)):
            count[0] += 1
            if listB[i] not in MISSING_CALLS and listX[i] not in MISSING_CALLS:
                polyMorphicList.append(markerIDs[i] + "[" + "-:-" + "/" + listB[i] + "/" + listX[i] + "]")
                parB = listB[i].split(":")
                fone = listX[i].split(":")
                if parB[0] == parB[1]:
                    count[1] += 1
                    if listX[i] in MISSING_CALLS:
                        count[5] += 1
                    elif listX[i] == listB[i]:
                        count[4] += 1
                    elif fone[0] != fone[1]:
                        count[2] += 1
                    else:
                        count[3] += 1
        if (count[1] - count[5]) > 0:
            expHet = int((count[2] * 100) / (count[1] - count[5]))
        else:
            expHet = 0
        if expHet == 100:
            trueF1 = "SuccessfulF1"
        elif expHet > int(cutOff):
            trueF1 = "Inconclusive"
        else:
            trueF1 = "UnsuccessfulF1"
    outFileHandle.write(
            fOneName + GRIDSEP + parentA + GRIDSEP + parentB + GRIDSEP + GRIDSEP.join(map(str, count)) + GRIDSEP + str(expHet) + GRIDSEP + trueF1 + GRIDSEP + availability + GRIDSEP + SPACE.join(
                polyMorphicList))
    outFileHandle.write(NEWLINE)

def pedigreeVerification(outFile, parentInfo, processedData, markerIDs, cutOff):
    outFileHandle = open(outFile, 'w')
    outFileHandle.write(F1SUMMARYHEAD)
    outFileHandle.write(NEWLINE)
    for fOneName in parentInfo:
        parentA = parentInfo[fOneName][0]
        parentB = parentInfo[fOneName][1]
        # print("Processing F1: " + fOne)
        if fOneName not in processedData:
            print(fOneName + MAPSEP + parentA + MAPSEP + parentB + MAPSEP + "F1 Data Missing")
            continue
        elif parentA not in processedData and parentB not in processedData:
            print(fOneName + MAPSEP + parentA + MAPSEP + parentB + MAPSEP + "Both Data Missing")
        else:
            if parentA not in processedData:
                print(fOneName + MAPSEP + parentA + MAPSEP + parentB + MAPSEP + "Parent A Data Missing")
            if parentB not in processedData:
                print(fOneName + MAPSEP + parentA + MAPSEP + parentB + MAPSEP + "Parent B Data Missing")
        # if processedData[parentA] and processedData[parentB]:
            # Get list of polymorphic markers between the parents based on LGC data
            getPolymorphicMarkers(outFileHandle, processedData, parentA, parentB, fOneName, markerIDs, cutOff)
        # else:
        #     print("Skipping F1:" + fOneName)

def readMetaData(inFile):
    with open(inFile) as md:
        lineNo = 0
        sampleMap = defaultdict(list)
        parentMap = defaultdict(list)
        parentInfo = defaultdict(dict)
        parentList = []
        for line in md:
            line = line.strip()
            lineNo += 1
            if lineNo == 1:
                continue
            lineEntries = line.split(MAPSEP)
            if len(lineEntries) == 1:
                continue
            if lineEntries[1] == "--" or lineEntries[1] is None:
                continue
            sampleMap[lineEntries[0]].append(lineEntries[2])
            parentMap[lineEntries[1]].append(lineEntries[2])
            if len(lineEntries) > 3:
                if lineEntries[3] == "--" or lineEntries[3] is None:
                    continue
                parentInfo[lineEntries[2]] = [lineEntries[3], lineEntries[4]]
                if lineEntries[3] not in parentList:
                    parentList.append(lineEntries[3])
                if lineEntries[4] not in parentList:
                    parentList.append(lineEntries[4])
    print("LOG: Processed " + str(lineNo - 1) + " entries from file " + inFile + " .!")
    return sampleMap, parentMap, parentInfo, parentList

def readGridCSV(inFile, selectedMarkers):
    gtData = defaultdict(list)
    dataFlag = 0
    markerIDs = []
    with open(inFile) as fh:
        lineCount = 0
        markerIDs = []
        indexList = []
        bases = init_bases()
        for line in fh:
            line = line.strip()
            lineCount += 1
            lineEntries = line.split(GRIDSEP)
            if lineEntries[0] == "DNA \\ Assay":
                dataFlag = 1
                if selectedMarkers:
                    for i in range(1, len(lineEntries)):
                        if lineEntries[i] in selectedMarkers:
                            markerIDs.append(lineEntries[i])
                            indexList.append(i)
                else:
                    markerIDs = lineEntries[1:]
                continue
            if line == "":
                continue
            if dataFlag == 0:
                continue
            if lineEntries[0] == "NTC":
                continue
            if selectedMarkers:
                for j in range(0, len(indexList)):
                    i = indexList[j]
                    if not lineEntries[i] == "":
                        if lineEntries[i] in MISSING_CALLS:
                            lineEntries[i] = "N:N"
                        else:
                            lineEntries[i] = str(lineEntries[i])
                    else:
                        lineEntries[i] = "N:N"
                    gtData[lineEntries[0]].append(lineEntries[i])
            else:
                for i in range(1, len(lineEntries)):
                    if not lineEntries[i] == "":
                        if lineEntries[i] in MISSING_CALLS:
                            lineEntries[i] = "N:N"
                        else:
                            lineEntries[i] = str(lineEntries[i])
                    else:
                        lineEntries[i] = "N:N"
                gtData[lineEntries[0]] = lineEntries[1:]
    return gtData, markerIDs

def readMetaDataFile(gtData, inFile):
    parentList = []
    sampleMap = defaultdict(list)
    parentInfo = defaultdict(list)
    parentMap = defaultdict(list)
    with open(inFile) as fh:
        lineCount = 0
        for line in fh:
            line = line.strip()
            lineCount += 1
            if lineCount == 1:
                continue
            lineEntries = line.split(MAPSEP)
            # sampleID[0], designation[1], sampleName[2], parentA[3], parentB[4]
            if lineEntries[0] not in gtData:
                continue
            if len(lineEntries) < 3:
                print("No entry for " + lineEntries[0] + " in metaData. Skipping..!")
                continue
            # print(lineEntries)
            sampleMap[lineEntries[0]] = lineEntries[2]
            parentMap[lineEntries[1]].append(lineEntries[2])
            if len(lineEntries) == 5:
                parentInfo[lineEntries[2]] = [lineEntries[3], lineEntries[4]]
                parentList.append(lineEntries[3])
                parentList.append(lineEntries[4])
        parentList = list(set(parentList))
        for parentName in parentList:
            if parentName not in parentMap:
                print("Do not have data for the Parent: ", parentName)
    return sampleMap, parentMap, parentInfo, parentList

def updateGTdata(gtData, sampleMap):
    updatedGTdata = defaultdict(list)
    for sampleID in sampleMap:
        if sampleID in gtData:
            updatedGTdata[sampleMap[sampleID]] = gtData[sampleID]
        else:
            print("No entry for " + sampleID + " in metaData. Skipping..!")
    return updatedGTdata

def readSnpFile(inFile):
    markerList = []
    with open(inFile) as fh:
        for line in fh:
            line = line.strip()
            if line not in markerList:
                markerList.append(line)
    return markerList

def Average(lst):
    return sum(lst) / len(lst)

def markerAssistedSelection(outFile, processedData, markerIDs, qtlData, favAlleleData, parentMap):
    outFileHandle = open(outFile, 'w')
    # outFileHandle.write("SAMPLE_NAME" + GRIDSEP 
    #     + "QTL_Name" + GRIDSEP 
    #     + "noMarkers" + GRIDSEP 
    #     + "noMissing" + GRIDSEP 
    #     + "PartialMatch" + GRIDSEP 
    #     + "CompleteMatch" + GRIDSEP 
    #     + "FAV|GT" + NEWLINE)
    outLine = []
    outLine.append("SAMPLE_NAME")
    outLine.append("PARENT_A")
    outLine.append("PARENT_B")
    outLine.append("TotalMarkers")
    outLine.append("TotalFavourable")
    outLine.append("TotalMissing")
    outLine.append("AvgPartialMatch")
    outLine.append("AvgCompleteMatch")
    for qtl in qtlData:
        outLine.append(qtl + "_noMarker")
        outLine.append(qtl + "_noMissing")
        outLine.append(qtl + "_partialMatch")
        outLine.append(qtl + "_completeMatch")
        outLine.append(qtl + "_markerString")
    # outLine.append("AveragePartialMatch")
    # outLine.append("AverageCompleteMatch")
    outFileHandle.write(CSV.join(outLine) + NEWLINE)
    for sample in processedData:
        outLine = []
        alleleCalls = processedData[sample]
        totMarker = 0
        totFav = 0
        totMissing = 0
        # totHet = 0
        outLine.append(sample)
        if sample in parentMap:
            outLine.append(parentMap[sample][0])
            outLine.append(parentMap[sample][1])
        else :
            outLine.append("na")
            outLine.append("na")
        for qtl in qtlData:
            markerScore = []
            noMissing = 0
            noMarker = 0
            partialMatch = 'NaN'
            completeMatch = 'NaN'
            markerString = ""
            for marker in qtlData[qtl]:
                if marker in markerIDs:
                    noMarker += 1
                    totMarker += 1
                    index = markerIDs.index(marker)
                    if markerString:
                        markerString = markerString + "; " + favAlleleData[marker].split(":")[0] + "|" + alleleCalls[index]
                    else:
                        markerString = favAlleleData[marker].split(":")[0] + "|" + alleleCalls[index]
                    if alleleCalls[index] == favAlleleData[marker]:
                        totFav += 1
                        markerScore.append(1)
                    elif alleleCalls[index].split(":")[0] == favAlleleData[marker].split(":")[0]:
                        markerScore.append(0.5)
                    elif alleleCalls[index].split(":")[1] == favAlleleData[marker].split(":")[1]:
                        markerScore.append(0.5)
                    elif alleleCalls[index] == 'N:N':
                        # markerScore.append(0)
                        totMissing += 1
                        noMissing += 1
                    else:
                        markerScore.append(0)
            if markerScore:
                partialMatch = Average(markerScore)
                completeMatch = min(markerScore)
            outLine.append(noMarker)
            outLine.append(noMissing)
            outLine.append(partialMatch)
            outLine.append(completeMatch)
            outLine.append(markerString)
            # outFileHandle.write(sample + GRIDSEP 
            #     + qtl + GRIDSEP 
            #     + str(noMarker) + GRIDSEP 
            #     + str(noMissing) + GRIDSEP 
            #     + str(partialMatch) + GRIDSEP 
            #     + str(completeMatch) + GRIDSEP
            #     + markerString + NEWLINE)
        parMatchAvg = 0
        completeMatchAvg = 0
        x = 0
        for i in range(5,len(outLine)):
            if i % 5 == 0 and i > 0:
                if outLine[i] == 'NaN':
                    continue
                j = i + 1
                x += 1
                parMatchAvg += outLine[i]
                completeMatchAvg += outLine[j]
        outLine.insert(3, totMarker)
        outLine.insert(4, totFav)
        outLine.insert(5, totMissing)
        outLine.insert(6, round(parMatchAvg/x, 2))
        outLine.insert(7, round(completeMatchAvg/x, 2))
        # outLine.append(round(parMatchAvg/x, 2))
        # outLine.append(round(completeMatchAvg/x, 2))
        outFileHandle.write(CSV.join(map(str, outLine)) + NEWLINE)
    outFileHandle.close()

def readQtlFile(inFile):
    qtlData = defaultdict(list)
    favAlleleData = defaultdict(list)
    lineNo = 0
    with open(inFile) as fh:
        for line in fh:
            lineEntries = []
            line = line.strip()
            if line.startswith('#'):
                continue
            lineNo += 1
            if lineNo == 1:
                continue
            lineEntries = line.split(TAB)
            qtlData[lineEntries[0]].append(lineEntries[1])
            favAlleleData[lineEntries[1]] = str(lineEntries[4])+":"+str(lineEntries[4])
    return qtlData,favAlleleData

def initializeList(listName, listLength):
    for i in range(1, listLength+1):
        # print(i, listLength)
        listName.append('N:N')
    return listName

def readLGCcsv(lgcFile):
    gtData = defaultdict(list)
    markerIDs = []
    snpFlag = 0
    dataFlag = 0
    noMarkers = 0
    noData = 0
    with open(lgcFile) as fh:
        for line in fh:
            line = line.strip()
            if line == 'SNPs':
                snpFlag = 1
                continue
            if line == "":
                snpFlag = 0
                dataFlag = 0
                continue
            if line == 'Data':
                dataFlag = 1
                continue
            if snpFlag == 1:
                noMarkers += 1
                if noMarkers == 1:
                    continue
                lineEntries = line.split(CSV)
                markerIDs.append(lineEntries[0])
            if dataFlag == 1:
                noData += 1
                if noData == 1:
                    continue
                lineEntries = line.split(CSV)
                if lineEntries[7] in MISSING_CALLS:
                    continue
                if lineEntries[7] not in gtData:
                    gtData[lineEntries[7]] = initializeList(gtData[lineEntries[7]], len(markerIDs))
                if lineEntries[6] not in markerIDs:
                    print(lineEntries[6], "Marker infor does not present in data. Please check the data.")
                idx = markerIDs.index(lineEntries[6])
                gtData[lineEntries[7]][idx] = lineEntries[3]
    return gtData, markerIDs

parser = argparse.ArgumentParser(description="QC pipeline: Processes the LGC file for tasks involved in purity check")
parser.add_argument("--task", dest="task", type=str, metavar="<STRING>", help="Task to be performed. Tasks are: "+SPACE.join(tasks))
parser.add_argument("--lgc-file", dest="lgcFile", metavar="<FILENAME>", help="LGC raw data File")
parser.add_argument("--grid-file", dest="gridFile", metavar="<FILENAME>", help="LGC Grid Matrix File")
parser.add_argument("--grid-files", dest="gridFiles", metavar="<FILENAME>", help="Comma seperated LGC Grid Matrix Files")
parser.add_argument("--meta-data", dest="metaDataFile", metavar="<FILENAME>", help="Parent information file")
parser.add_argument("--marker-list", dest="snpFile", metavar="<FILENAME>", help="File with list of snps")
parser.add_argument("--qtl-file", dest="qtlFile", metavar="<FILENAME>", help="QTL file in GOBii format")
parser.add_argument("--out", dest="outPrefix", metavar="<STRING>", help="Output filename prefix")
parser.add_argument("--f1het-cutOff", dest="cutOff", default=60, metavar="<Integer>", help="Percentage expected heterozygosity for F1 verification")
parser.add_argument("--consensus-cutOff", dest="consensusCutOffInput", default=50, metavar="<Integer>", help="Percentage propotion to be considered to call consensus")

# Parse commandline arguments
options = parser.parse_args()

# Task to be performed
task = options.task

# Return help message if the task is not provided or given wrong
if task not in tasks:
    print("Invalid/No Task provided.")
    parser.print_help()

if task == "F1CHECK":
    required = "gridFile metaDataFile outPrefix".split(SPACE)
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)
            
    gridFile = options.gridFile
    metaDataFile = options.metaDataFile
    outPrefix = options.outPrefix
    cutOff = options.cutOff
    if cutOff is None:
            cutOff = 100

    snpFile = options.snpFile
    if snpFile is not None:
        selectedMarkers = readSnpFile(snpFile)
    
    consensusCutOffInput = options.consensusCutOffInput
    if consensusCutOffInput is not None:
        consensusCutOff = int(consensusCutOffInput) / 100

    gtData, markerIDs = readGridCSV(gridFile, selectedMarkers)
    sampleMap, parentMap, parentInfo, parentList = readMetaDataFile(gtData, metaDataFile)
    # sampleMap, sampleList = readMapFile(sampleMapFile)
    # parentMap, parentList = readMapFile(parentMapFile)
    # markerIDs, gtData = readGridFile(gridFile, sampleMap)
    gtData = updateGTdata(gtData, sampleMap)
    consensusData = getConsensusCalls(parentList, parentMap, markerIDs, gtData, outPrefix + "_" + "ConsensusSummary.csv", consensusCutOff, task)
    processedData = {**gtData, **consensusData}

    consensusOutData = defaultdict(list)
    for parent in consensusData:
        consensusOutData[parent] = consensusData[parent]
        for sampleID in parentMap[parent]:
            consensusOutData[sampleID] = gtData[sampleID]
            if sampleID in processedData:
                processedData.pop(sampleID)
            
    writeGridFile(outPrefix + "_" + "Grid_consensus.csv", consensusOutData, markerIDs)
    writeFjGTFile(outPrefix + "_" + "Grid_processed_FJ.data", processedData, markerIDs)
    writeGridFile(outPrefix + "_" + "Grid_processed.csv", processedData, markerIDs)
    writePurityTable(outPrefix + "_" + "Puritytable.csv", parentMap, consensusOutData)
    pedigreeVerification(outPrefix + "_" + "F1summary.csv", parentInfo, processedData, markerIDs, cutOff)

if task == "CONSENSUS":
    required = "gridFile metaDataFile outPrefix".split(SPACE)
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)
    
    gridFile = options.gridFile
    outPrefix = options.outPrefix
    metaDataFile = options.metaDataFile
    # sampleMapFile = options.sampleMap
    # parentMapFile = options.parentMap

    snpFile = options.snpFile
    if snpFile is not None:
        selectedMarkers = readSnpFile(snpFile)
    
    consensusCutOffInput = options.consensusCutOffInput
    if consensusCutOffInput is not None:
        consensusCutOff = int(consensusCutOffInput) / 100

    gtData, markerIDs = readGridCSV(gridFile, selectedMarkers)
    sampleMap, parentMap, parentInfo, parentList = readMetaDataFile(gtData, metaDataFile)
    # sampleMap, sampleList = readMapFile(sampleMapFile)
    # parentMap, parentList = readMapFile(parentMapFile)
    # markerIDs, gtData = readGridFile(gridFile, sampleMap)
    gtData = updateGTdata(gtData, sampleMap)
    consensusData = getConsensusCalls(parentList, parentMap, markerIDs, gtData, outPrefix + "_" + "ConsensusSummary.csv", consensusCutOff, task)
    processedData = {**gtData, **consensusData}

    consensusOutData = defaultdict(list)
    for parent in consensusData:
        consensusOutData[parent] = consensusData[parent]
        for sampleID in parentMap[parent]:
            consensusOutData[sampleID] = gtData[sampleID]
            if sampleID in processedData:
                processedData.pop(sampleID)
    
    writeGridFile(outPrefix + "_" + "Grid_consensus.csv", consensusOutData, markerIDs)
    writeFjGTFile(outPrefix + "_" + "Grid_processed_FJ.data", processedData, markerIDs)
    writeGridFile(outPrefix + "_" + "Grid_processed.csv", processedData, markerIDs)
    writePurityTable(outPrefix + "_" + "Puritytable.csv", parentMap, consensusOutData)

if task == "GETSTATS":
    required = "gridFile"
    if options.__dict__[required] is None:
        print("Option/Value Missing for", required)
        parser.print_help()
        sys.exit(1)
    
    gridFile = options.gridFile
    snpFile = options.snpFile
    if snpFile is not None:
        selectedMarkers = readSnpFile(snpFile)

    gtData, markerIDs = readGridCSV(gridFile, selectedMarkers)
    print("Number of Markers: ", len(markerIDs))
    print("Number of samples:", len(gtData.keys()))

if task == "LGC2FJK":
    required = "gridFile outPrefix".split(SPACE)
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)
    
    gridFile = options.gridFile
    outPrefix = options.outPrefix
    snpFile = options.snpFile
    if snpFile is not None:
        selectedMarkers = readSnpFile(snpFile)

    gtData, markerIDs = readGridCSV(gridFile, selectedMarkers)
    writeFjGTFile(outPrefix + "_FJ.data", gtData, markerIDs)

if task == "RENAME":
    required = "gridFile metaDataFile outPrefix".split(SPACE)
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)
    
    snpFile = options.snpFile
    if snpFile is not None:
        selectedMarkers = readSnpFile(snpFile)

    gridFile = options.gridFile
    outPrefix = options.outPrefix
    metaDataFile = options.metaDataFile    
    gtData, markerIDs = readGridCSV(gridFile, selectedMarkers)
    sampleMap, parentMap, parentInfo, parentList = readMetaDataFile(gtData, metaDataFile)
    gtData = updateGTdata(gtData, sampleMap)
    writeFjGTFile(outPrefix + "_FJ.data", gtData, markerIDs)
    writeGridFile(outPrefix + "_Grid.csv", gtData, markerIDs)

if task == "MERGEDATA":
    required = "gridFiles outPrefix".split(SPACE)
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)
    gridFiles = options.gridFiles
    outPrefix = options.outPrefix
    outGridFile = outPrefix + "_Grid.csv"
    outFJKFile = outPrefix + "_FJK.data"
    snpFile = options.snpFile
    selectedMarkers = []
    if snpFile is not None:
        selectedMarkers = readSnpFile(snpFile)

    gtDataList = []
    markerIdsList = []

    # gtData, markerIDs = readGridCSV(gridFile, selectedMarkers)
    litFiles = gridFiles.split(",")
    for fileName in litFiles:
        # print(fileName)
        if selectedMarkers:
            gtData, markerIDs = readGridCSV(fileName, selectedMarkers)
        else:
            gtData, markerIDs = readGridCSV(fileName, None)
            print(fileName, len(markerIDs))
        gtDataList.append(gtData)
        markerIdsList.append(markerIDs)

    firstList = []
    for  listMarkers in markerIdsList:
        secondList = listMarkers

        set_1 = set(firstList)
        set_2 = set(secondList)

        list_2_items_not_in_list_1 = list(set_2 - set_1)
        firstList = firstList + list_2_items_not_in_list_1
    combinedMarkerList = firstList
    combinedMarkerList.sort()
    combinedGtData = defaultdict(list)
    for fileName in litFiles:
        gtData = defaultdict(list)
        markerIDs = []
        gtData, markerIDs = readGridCSV(fileName, None)
        for sample in gtData:
            genotypeData = gtData[sample]
            if sample not in combinedGtData:
                combinedGtData[sample] = []
                combinedGtData[sample] = initializeList(combinedGtData[sample], len(combinedMarkerList))
            for marker in markerIDs:
                if marker in combinedMarkerList:
                    q_idx = markerIDs.index(marker)
                    s_idx = combinedMarkerList.index(marker)
                    combinedGtData[sample][s_idx] = genotypeData[q_idx]
    writeGridFile(outGridFile, combinedGtData, combinedMarkerList)
    writeFjGTFile(outFJKFile, combinedGtData, combinedMarkerList)

if task == "FWDBREED":
    required = "gridFile metaDataFile qtlFile outPrefix".split(SPACE)
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)
            
    gridFile = options.gridFile
    metaDataFile = options.metaDataFile
    qtlFile = options.qtlFile
    outPrefix = options.outPrefix
    cutOff = options.cutOff
    if cutOff is None:
            cutOff = 100

    snpFile = options.snpFile
    if snpFile is not None:
        selectedMarkers = readSnpFile(snpFile)

    consensusCutOffInput = options.consensusCutOffInput
    if consensusCutOffInput is not None:
        consensusCutOff = int(consensusCutOffInput) / 100


    qtlData, favAlleleData = readQtlFile(qtlFile)
    gtData, markerIDs = readGridCSV(gridFile, selectedMarkers)
    sampleMap, parentMap, parentInfo, parentList = readMetaDataFile(gtData, metaDataFile)
    # sampleMap, sampleList = readMapFile(sampleMapFile)
    # parentMap, parentList = readMapFile(parentMapFile)
    # markerIDs, gtData = readGridFile(gridFile, sampleMap)
    gtData = updateGTdata(gtData, sampleMap)
    consensusData = getConsensusCalls(parentList, parentMap, markerIDs, gtData, outPrefix + "_" + "ConsensusSummary.csv", consensusCutOff, task)
    processedData = {**gtData, **consensusData}

    consensusOutData = defaultdict(list)
    for parent in consensusData:
        consensusOutData[parent] = consensusData[parent]
        for sampleID in parentMap[parent]:
            consensusOutData[sampleID] = gtData[sampleID]
            if sampleID in processedData:
                processedData.pop(sampleID)
            
    writeGridFile(outPrefix + "_" + "Grid_consensus.csv", consensusOutData, markerIDs)
    writeFjGTFile(outPrefix + "_" + "Grid_processed_FJ.data", processedData, markerIDs)
    writeGridFile(outPrefix + "_" + "Grid_processed.csv", processedData, markerIDs)
    writePurityTable(outPrefix + "_" + "Puritytable.csv", parentMap, consensusOutData)
    markerAssistedSelection(outPrefix + "_" + "MAS.csv", processedData, markerIDs, qtlData, favAlleleData, parentInfo)

if task == "LGC2GRID":
    required = "lgcFile outPrefix".split(SPACE)
    for req in required:
        if options.__dict__[req] is None:
            print("Option/Value Missing for", req)
            parser.print_help()
            sys.exit(1)
    lgcFile = options.lgcFile
    outPrefix = options.outPrefix
    outGridFile = outPrefix + "_Grid.csv"

    gtData, markerIDs = readLGCcsv(lgcFile) 
    writeGridFile(outGridFile, gtData, markerIDs)

