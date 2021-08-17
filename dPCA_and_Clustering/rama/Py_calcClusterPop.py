#Call python Py_calcClusterPop.py seqName s1/s2
#
#For example, python calcClusterPop.py AGAGAG s1

import os
import sys

dirname = sys.argv[2]+"c"+sys.argv[1]+"_phipsi"
files = os.listdir(dirname)

totNumFrames = 0
clustFrames = {}
for name in files:
	if name[-3:] == "txt" and name[0:3] != "log":
		clusterFile = open(dirname + "/" + name, "r")
		clustFrames[int(name.split(".")[0][7:])] = 0
		for line in clusterFile:
			clustFrames[int(name.split(".")[0][7:])] += 1
			totNumFrames += 1

clustDict = {}
dictList = []
sortList = []
for index, clust in enumerate(clustFrames.keys()):
	clustDict[clust] = float(clustFrames[clust]) / float(totNumFrames)
	dictList.append(index)

numClusters = len(dictList)
for i in range(numClusters):
	currMaxVal = -1
	currMaxClust = -1
	for j in range(len(dictList)):
		if clustDict[dictList[j]] > currMaxVal:
			currMaxClust = dictList[j]
			currMaxVal = clustDict[dictList[j]]
	sortList.append(currMaxClust)
	dictList.remove(currMaxClust)


clustCount = 1
for index, num in enumerate(sortList):
	if num != 0:
		os.rename(sys.argv[2] + "c" + sys.argv[1] + "_cluster" + str(num) + ".png", sys.argv[2] + sys.argv[1] + "_cluster" + str(clustCount) + ".png")
		print "cluster " + str(clustCount) + ": " + str(clustDict[num])
		clustCount += 1
	else:
		os.rename(sys.argv[2] + "c" + sys.argv[1] + "_cluster" + str(num) + ".png", sys.argv[2] + sys.argv[1] + "_cluster" + str(0) + ".png")
		print "cluster " + str(0) + ": " + str(clustDict[num])
