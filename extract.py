import sys
from collections import defaultdict

inMapFile = sys.argv[1]
inGridFile = sys.argv[2]
outGridFile = sys.argv[3]
TAB = "\t"

missingData = ('NTC', 'EMPTY')
gtData = defaultdict(list)
markerIds = []
dataFlag = 0
with open(inGridFile) as fh:
	for line in fh:
		line = line.strip()
		lineEntries = line.split(",")

		if lineEntries[0] == "DNA \\ Assay":
			dataFlag = 1
			markerIds = lineEntries[1:]
			continue

		if dataFlag == 0:
			continue

		if lineEntries[0] in missingData:
			continue

		gtData[lineEntries[0]] = lineEntries[1:]

with open(inMapFile) as fh:
	with open(outGridFile, 'w') as outFileHandle:
		outFileHandle.write("DNA \\ Assay" + "," + ",".join(markerIds) + "\n")
		for line in fh:
			line = line.strip()
			lineEntries = line.split(TAB)
			if lineEntries[0] in gtData:
				outFileHandle.write(lineEntries[0] + "," + ",".join(gtData[lineEntries[0]]) + "\n")
			else:
				print(lineEntries[0] + " Does not exist.")
	outFileHandle.close()