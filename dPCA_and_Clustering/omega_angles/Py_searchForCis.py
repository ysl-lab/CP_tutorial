#Aidan Fike
#June 28, 2019

#Program to go through measured omega angle trajectories and find all cis peptide bonds 

import os

cisFound = False
for direct in ["s1", "s2"]:
	for filename in os.listdir(direct):
		currFile = open(direct+"/"+filename, "r")
		for line in currFile:
			words = line.split()
			if words[0] != "#" and words[0] != "@":
				for index, word in enumerate(words):
					if index > 1 and abs(float(word)) < 90:
						print "Cis bond found!", direct, filename, words[0], index - 2
						cisFound = True
if not cisFound:
	print "No cis angles found" 
