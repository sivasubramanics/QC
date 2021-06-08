from collections import defaultdict

MAPSEP = "\t"
NEWLINE = "\n"
SPACE = " "
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
                + "Marker [CallA/CallB/CallF1]"

IUPAC = ("A", "T", "G", "C")


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
    """
    
    """
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
    tmpSeq = list(seq)
    return "/".join(tmpSeq)

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
            outLine.append(addColon(bases[baseCall]))
        outFileHandle.write(GRIDSEP.join(outLine))
        outFileHandle.write(NEWLINE)

def writeFjGTFile(outFile, gtData, markerIDs):
    # Opening output file handle
    outFileHandle = open(outFile, 'w')
    bases = init_bases()
    # Build and write the first/header line for Flapjack Genotype
    outFileHandle.write(FJGHEAD + NEWLINE)
    outFileHandle.write(MAPSEP + MAPSEP.join(markerIDs))
    for sampleID in gtData:
        outLine = []
        outLine.append(sampleID)
        for baseCall in gtData[sampleID]:
            outLine.append(addSlash(bases[baseCall]))
        outFileHandle.write(MAPSEP.join(outLine))
        outFileHandle.write(NEWLINE)

def writePurityTable(outFile, parentMap, consensusOutData, markerIDs):
    outFileHandle = open(outFile, 'w')
    outFileHandle.write("SAMPLE_NAME" + GRIDSEP + "total" + GRIDSEP + "match" + GRIDSEP + "mismatch" + GRIDSEP + "bothN" + GRIDSEP + "conN" + GRIDSEP + "sampleN")
    outFileHandle.write(NEWLINE)
    for parentID in parentMap:
        for sampleID in parentMap[parentID]:
            if sampleID not in consensusOutData:
                continue
            # total, match, mismatch, bothN, conN, sampleN
            count = [0, 0, 0, 0, 0, 0]
            for i in range(0, len(markerIDs)):
                count[0] += 1
                if consensusOutData[parentID][i] == consensusOutData[sampleID][i]:
                    if consensusOutData[parentID][i] == 'N':
                        count[3] += 1
                    else:
                        count[1] += 1
                else:
                    if consensusOutData[parentID][i] == 'N':
                        count[4] += 1
                    else:
                        if consensusOutData[sampleID][i] == 'N':
                            count[5] += 1
                        else:
                            count[2] += 1
            outFileHandle.write(sampleID + GRIDSEP + GRIDSEP.join(map(str, count)))
            outFileHandle.write(NEWLINE)

def readGridFile(inFile, sampleMap):
    bases = init_bases()
    gtData = defaultdict(list)
    markerIDs = []
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

def getConsensusCalls(parentList, parentMap, markerIDs, gtData, outFile, task):
    consensusData = defaultdict(list)
    outFileHandle = open(outFile, 'w')
    outFileHandle.write("parent" + GRIDSEP + "noReps" + GRIDSEP + "noMarkers" + GRIDSEP + "allNs" + GRIDSEP + "noMajorCalls" + GRIDSEP + "noMinorCalls")
    outFileHandle.write(NEWLINE)
    for parentID in parentMap:
        if task == "F1CHECK":
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
                if baseCall == 'N':
                    countNs += 1
                else:
                    if baseCall in baseCount:
                        baseCount[baseCall] += 1
                    else:
                        baseCount[baseCall] = 1
            if countNs == len(parentMap[parentID]):
                consensusBase.append('N')
                count[0] += 1
            else:    
                # baseCount = sort_dict(baseCount)
                # print(baseCount)
                sortedBases = sorted(baseCount, key=baseCount.get, reverse=True)
                # print(sortedBases)
                # for base in sorted(baseCount, key=baseCount.get, reverse=True):
                for base in sortedBases:
                    if base == 'N':
                        continue
                    if baseCount[base]/(len(parentMap[parentID]) - countNs) > 0.5:
                        consensusBase.append(base)
                        count[1] += 1
                    else:
                        if baseCount[base]/(len(parentMap[parentID]) - countNs) == 0.5:
                            nextbase = next_key(sortedBases, base)
                            if baseCount[base] == baseCount[nextbase]:
                                consensusBase.append('N')
                                count[2] += 1
                            else:
                                consensusBase.append(base)
                                count[1] += 1
                        else:
                            consensusBase.append('N')
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
    listA = processedData[parentA]
    listB = processedData[parentB]
    listX = processedData[fOneName]
    polyMorphicList = []
    bases = init_bases()
    # total(0), countPM(1), countHet(2), countA(3), countB(4), countN(5)
    count = [0, 0, 0, 0, 0, 0]
    for i in range(0, len(listX)):
        count[0] += 1
        if IUPAC.__contains__(listA[i]) and IUPAC.__contains__(listB[i]):
            if listA[i] != listB[i]:
                polyMorphicList.append(markerIDs[i] + "[" + listA[i] + "/" + listB[i] + "/" + listX[i] + "]")
                count[1] += 1
                if bases[listA[i] + listB[i]] == listX[i]:
                    count[2] += 1
                if listA[i] == listX[i]:
                    count[3] += 1
                if listB[i] == listX[i]:
                    count[4] += 1
                if listX[i] == "N":
                    count[5] += 1
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
            fOneName + GRIDSEP + parentA + GRIDSEP + parentB + GRIDSEP + GRIDSEP.join(map(str, count)) + GRIDSEP + str(expHet) + GRIDSEP + trueF1 + GRIDSEP + SPACE.join(
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
        if parentA not in processedData:
            print(fOneName + MAPSEP + parentA + MAPSEP + parentB + MAPSEP + "Parent A Data Missing")
            continue
        if parentB not in processedData:
            print(fOneName + MAPSEP + parentA + MAPSEP + parentB + MAPSEP + "Parent B Data Missing")
            continue
        if processedData[parentA] and processedData[parentB]:
            # Get list of polymorphic markers between the parents based on LGC data
            getPolymorphicMarkers(outFileHandle, processedData, parentA, parentB, fOneName, markerIDs, cutOff)
        else:
            print("Skipping F1:" + fOneName)

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

def readGridCSV(inFile):
    gtData = defaultdict(list)
    markerIDs = []
    with open(inFile) as fh:
        lineCount = 0
        markerIDs = []
        bases = init_bases()
        for line in fh:
            line = line.strip()
            lineCount += 1
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
            for i in range(1, len(lineEntries)):
                if not lineEntries[i] == "":
                    lineEntries[i] = bases[lineEntries[i].replace(":", "")]
                else:
                    lineEntries[i] = bases['NN']
            gtData[lineEntries[0]]=lineEntries[1:]
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
            sampleMap[lineEntries[0]] = lineEntries[2]
            parentMap[lineEntries[1]].append(lineEntries[2])
            if len(lineEntries) == 5:
                parentInfo[lineEntries[2]] = [lineEntries[3], lineEntries[4]]
                parentList.append(lineEntries[3])
                parentList.append(lineEntries[4])
    return sampleMap, parentMap, parentInfo, parentList

def updateGTdata(gtData, sampleMap):
    updatedGTdata = defaultdict(list)
    for sampleID in sampleMap:
        if sampleID in gtData:
            updatedGTdata[sampleMap[sampleID]] = gtData[sampleID]
        else:
            print("No entry for " + sampleID + " in metaData. Skipping..!")
    return updatedGTdata
