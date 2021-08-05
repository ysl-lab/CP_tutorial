#Aidan Fike
#June 27, 2019
#
#Program to find the file where the first cis peptide bond(s) were formed

cisFile = open("cisOut.txt", "r")

firstTime = [10000000]*6
firstTraj = [-1]*6

lastTime = [0]*6
lastTraj = [-1]*6

for line in cisFile:
	words=line.split()
	if len(words) != 4:
		if float(words[5]) < firstTime[int(words[6])]:
			firstTime[int(words[6])] = float(words[5])
			firstTraj[int(words[6])] = words[3] + "/" + words[4]
			
		if float(words[5]) > lastTime[int(words[6])]:
			lastTime[int(words[6])] = float(words[5])
			lastTraj[int(words[6])] = words[3] + "/" + words[4]
	else:
		print "No cis angles found"

for index, name in enumerate(firstTraj):
	if name != -1:
		print "Cis dihedral found at position ", index
		print "firstTraj: ", firstTraj[index], firstTime[index]
		print "lastTraj: ", lastTraj[index], lastTime[index], "\n"
