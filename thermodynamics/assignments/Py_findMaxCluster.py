with open("assignments.txt", "r") as asFile:
	maxClus = 0
	for line in asFile:
		words = line.split()
		if int(words[4]) > maxClus:
			maxClus = int(words[4])
	
	print "maxClus", maxClus 
