#Aidan Fike
#July 21, 2019

#File used to create random dihedral angles and random sequences

import random
import sys
import subprocess 
import os
from shutil import copyfile

random.seed(a=None)

#Conversion between ints and amino acid names (used for random sequence creation)
def getNumLetterConv():
    return {0: 'G', 1: 'A', 2: 'V', 3: 'F', 4: 'N', 5: 'S', 6: 'R', 7: 'D'}

#All letters of valid amino acids
def getAllAminos():
    return ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

#See if a sequence is made up of only valid amino acids
def isAllowedSeq(seq):
    aminos = getAllAminos()

    isAllowed = True
    for i in range(len(seq)):
        try:
            if not seq[i].upper() in aminos: 
                isAllowed = False
        except: #If anything goes wrong, the sequence is probably incorrect.
            print("Something weird happened when evaluating the sequence. Please double-check your input and try again.")
            isAllowed = False

    return isAllowed

#Create a new random sequence. Look through previously created sequences to make sure there is no
#redundancy
def createNewSeq(numAmino, readPrev=True):
    currSeq = ""

    #Find already-created sequences and their cyclic equivalents
    prevSeqs = []
    if readPrev:
        if os.path.exists("finSeqs.txt"):
            seqFile = open("finSeqs.txt", "r")
            for line in seqFile:
                currSeq = line
                if len(line.split()[0]) == numAmino:
                    for i in range(numAmino):
                        cyclicEquiv = ""
                        for j in range(numAmino):
                            cyclicEquiv += currSeq[(j + i) % numAmino] 
                        prevSeqs.append(cyclicEquiv)
            seqFile.close()

    numLetterConv = getNumLetterConv()    

    #Generate random sequences (making sure you aren't creating already-created structs)
    newSeq = ""
    badSeq = True
    while(badSeq):
        newSeq = ""
        doWhile = False
        badSeq = False
        for i in range(numAmino):
            randIndex = random.randrange(len(numLetterConv))
            newSeq += (numLetterConv[randIndex]) 
        for prevseq in prevSeqs:
            if newSeq == prevseq:
                badSeq = True

    print("\nThe sequence %s was choosen randomly\n\n" % newSeq)
    
    return newSeq

#Write random phi/psi angles to a file
def writeToDihedFile(seq, phi, psi, structNum):
    dihedFile = open("s" + str(structNum) + ".di", "w") 

    dihedFile.write(seq + "\n\n")
    dihedFile.write("s" + str(structNum) + " phi/psi:\n")
    for i in range(len(phi)):
	    dihedFile.write(str(phi[i]) + " " + str(psi[i]) + "\n")
    dihedFile.write("\n")   

    dihedFile.close()

#Create an array of random phi/psi angles
def createRandomPhiPsi(numDihedCombs, newSeq):
    numAminos = len(newSeq)
    phi = []
    psi = []
    for i in range(numDihedCombs):
        phi.append([])
        psi.append([])
        for j in range(numAminos):
            minPhi = -180
            maxPhi = 180
            minPsi = -180
            maxPsi = 180

            phi[i].append(random.randrange(minPhi, maxPhi))
            psi[i].append(random.randrange(minPsi, maxPsi))    
    return phi, psi
