#Aidan Fike
#July 21, 2019

#Edited and updated to version 2
#Kevin Schult, March/April 2020

#File with functions which measure omega angles and rmsd.
#Dependent on MDAnalysis library, which does not come standard in most Python installations

import os
import sys
import dihedDir
import shutil
import chimScriptMaker

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
import numpy as np
import networkx as nx

#Function that evaluates a structure to determine whether it is "valid,"
#that is, whether it meets the criteria of having the right peptide bonds
#and has a low head-to-tail distance.
#
#Returns boolean validPeptide, which is True when the structure is valid.
#Side effects: moves or deletes structure file depending on whether it is
#valid or invalid, respectively.
def checkPeptide(structNum, newSeq, bonds, verbose):
    flippedBondFound = checkPeptideBonds(structNum, newSeq, bonds, verbose)
    correctCyclization = checkCyclization(structNum, newSeq, verbose)
    chiralityError = checkChirality(structNum, newSeq, verbose)
    validPeptide = -1

    if correctCyclization and not flippedBondFound and not chiralityError:
        print("s%d is a valid structure.\n" % structNum)
        os.rename("s%d_int_noh.pdb" % structNum, "validStructs/s%d_int_noh.pdb" % structNum)
        validPeptide = True
    else:
        print("s%d is not a valid structure and will be discarded.\n" % structNum)
        os.rename("s%d_int_noh.pdb" % structNum, "s%d_invalid.pdb" % structNum)
        validPeptide = False

    return validPeptide

#Function to check if the chiralities are correct.
#Returns a boolean value. True if there's a problem, False otherwise.
def checkChirality(structNum, newSeq, verbose):
    chiralityTolerance = 10
    chiralityRange = 40

    u = mda.Universe("s%d_int_noh.pdb" % structNum)
    peptide = u.select_atoms('protein')

    length = len(newSeq)
    chiralityError = False
    if verbose:
        print("Chirality:")
    for i in range(length):
        #The idea here is to calculate the improper dihedral CA-N-C-CB
        #for each amino acid. If the dihedral is positive, residue is L,
        #otherwise residue is D. This is true even for cysteine, since
        #it's definted here by spatial orientation, not by R/S notation.
        #Glycine is excluded since it has no chirality to be checked.

        chirality = 'none' #correct initialization

        if newSeq[i].upper() == "G":
            if verbose:
                print(newSeq[i], "N/A")
            continue #skip to next residue.

        res = str(i+1)
        chiralDihedralGroup = peptide.select_atoms('resnum ' + res + ' and name CA')\
                            + peptide.select_atoms('resnum ' + res + ' and name N')\
                            + peptide.select_atoms('resnum ' + res + ' and name C')\
                            + peptide.select_atoms('resnum ' + res + ' and name CB')
        #print(chiralDihedralGroup, chiralDihedralGroup.atoms)
        chiralDihedral = chiralDihedralGroup.dihedral.value()

        if chiralityTolerance < chiralDihedral < chiralityTolerance+chiralityRange:
            chirality = "L"
        elif -1*chiralityTolerance > chiralDihedral > -1*(chiralityTolerance+chiralityRange):
            chirality = "D"

        chiralStr = formatDihedral(chiralDihedral) #returns a string with the right length
        #debugging
        if verbose:
            print(newSeq[i], chirality, chiralStr)

        if (chirality == "L" and newSeq[i].isupper()) or \
(chirality == "D" and newSeq[i].islower()):
            continue #This residue is okay, check the next one.
        else:
            chiralityError = True
            if verbose:
                print("Chirality error found at residue " + str(i+1))

    return chiralityError

#Function to measure omega angles of a given structure and determine
#if it has cis bond(s) in the right places.
#
#Params: structNum - int: The structure number being evaluated
#        newSeq - string: The sequence of the structure being evaluated 
#        bonds - string: Whether each peptide bond is cis or trans. "ttctt"
#                        is an example for a 5-mer.
#Be wary of the omega tolerance. A deviation larger than omegaTolerance from
#cis/trans omega angles will throw out the structure.
#
#Retruns boolean flippedBondFound, which is True when a bond is the opposite
#of its intended orientation.
def checkPeptideBonds(structNum, newSeq, bonds, verbose):
        omegaTolerance = 30.0
        u = mda.Universe("s%d_int_noh.pdb" % structNum)
        peptide = u.select_atoms('protein')
        omegas = [res.omega_selection() for res in peptide.residues]          
        #MDAnalysis can't naturally handle cyclization, so we have to spoon-feed it
        #the atom group for the omega angle between the last and first residues.
        #At the moment, the last element of omegas is <None>, so let's manually define it.
        omegas[-1] = peptide.select_atoms('resnum ' + str(len(newSeq)) + \
                                          ' and name CA') + \
                     peptide.select_atoms('resnum ' + str(len(newSeq)) + \
                                          ' and name C') + \
                     peptide.select_atoms('resnum 1' + \
                                          ' and name N') + \
                     peptide.select_atoms('resnum 1' + \
                                          ' and name CA')
    
        flippedBondFound = False
        for i in range(len(newSeq)):
            #debugging
            #print(omegas[i].dihedral.value(), i, newSeq[i])
            omegaValue = omegas[i].dihedral.value()
            omegaStr = formatDihedral(omegaValue) #returns a string
            bondType = "None"

            if -1*omegaTolerance < omegaValue < omegaTolerance:
                bondType = 'c' #cis
            elif (-180.0 < omegaValue < -180.0+omegaTolerance) or \
                 (180.0-omegaTolerance < omegaValue < 180.0):
                bondType = 't' #trans

            if bondType == bonds[i]:
                if verbose:
                    print("Omega %s at residue %d" % (omegaStr, i+1))
                continue #This bond checks out, check the next one.
            else:
                #Something is incorrect, and the structure will be discarded.
                flippedBondFound = True
                if verbose:
                    print("Omega %s at residue %d, incorrect bondtype %s"\
                          % (omegaStr, i+1, bondType))

        return flippedBondFound

def formatDihedral(angle):
    angleStr = "%.1f" % angle
    while len(angleStr) < 6: #max length is 6, e.g. -126.2 
        angleStr = " " + angleStr #formatting
    return angleStr

#Function to check cyclization via measuring head-to-tail distance.
#
#Returns boolean value correctCyclization which is False when the
#bond distance is greater than 1.4 A.
def checkCyclization(structNum, newSeq, verbose):
    maxBondLength = 1.4 #Angstroms, avg peptide bond is 1.3 A
    u = mda.Universe("s%d_int_noh.pdb" % structNum)
    peptide = u.select_atoms('protein')

    nTerm = peptide.select_atoms('resnum 1 and name N')
    cTerm = peptide.select_atoms('resnum ' + str(len(newSeq)) + ' and name C')

    distanceArray = mda.analysis.distances.dist(nTerm,cTerm)
    #distanceArray contains three objects:
    #resID of first atomGroup [0]
    #resID of second atomGroup [1]
    #distance array [2]
    headToTail = distanceArray[2][0]
    if verbose:
        print("s%d has a head-to-tail distance of %.2f A" % (structNum, headToTail))

    correctCyclization = (headToTail <= maxBondLength)
    if not correctCyclization:
        print("s%d did not cyclize correctly -- head-to-tail distance is too large" % structNum)
    
    return correctCyclization

def addToGraph(thresh, structNum, structs, includeOxygen):
    #Add new structure to graph of structures
    structs.add_node(structNum)
    #Calc rmsd between new structure and other nodes - add edge if over thresh
    structs = calcRmsd(structNum, thresh, structs, includeOxygen) 

    return structs

#Function to calulate the rmsd between all non-cis sequences
#
#Return: Graph including all structures as nodes and all RMSD>thresh as edges
def calcRmsd(structNum, thresh, structs, includeOxygen):
    u = {} #dictionary of validStructs structures being compared

    #MDAnalysis 'backbone' includes C,CA,N,O.
    #We usually want to just use C,CA,N to calculate rmsd.
    #includeOxygen flag switches to non-default option.
    atomGroup = 'backbone and not name O'
    if includeOxygen:
        atomGroup = 'backbone'

    u[0] = mda.Universe("validStructs/s%d_int_noh.pdb" % structNum)
    name2 = "s%d" % structNum
    #Only need to calculate rmsd with the newly generated structure, named u[0] here
    for strIndex in list(structs.nodes):
        if strIndex == structNum:
            continue #don't compare structure to itself
        u[strIndex] = mda.Universe("validStructs/s%d_int_noh.pdb" % strIndex)
        name1 = "s%d" % strIndex

        rmsValue = rms.rmsd(u[0].select_atoms(atomGroup).positions, \
                                        u[strIndex].select_atoms(atomGroup).positions, \
                                        center=True, \
                                        superposition=True)

        if rmsValue > thresh:
            structs.add_edge(structNum, strIndex, rms=rmsValue)

        print("\n%s and %s have a backbone rmsd of %f Angstroms" % (name1, \
                                name2, rmsValue))

    return structs
     
def checkSubGraphs(structNum, structs, setSize, verbose):
    foundSet = False
    maxClique = "Groups not evaluated, too few connected structures\n"
    #Sanity check - no need to check subgraphs if structNum isn't 
    #above RMSD threshold for enough structures
    if structs.degree(structNum) >= (setSize - 1):
        #Find largest set of structures that each meet the RMSD threshold
        completeSubgraphs = nx.make_clique_bipartite(structs)
        #print(structs.nodes)
        #print(completeSubgraphs.nodes.data())
        #Line above returns a new graph where each structure is connected
        #to the "maximal cliques" that include that structure.
        #(The "maximal clique" of a node /n/ is the largest set of nodes
        #containing /n/ where every node is adjacent to ever other node.)
        #We want to find the maximal clique of our new structure structNum
        #and check its size, but structNum can be involved in the
        #maximal cliques of other nodes, so we need to do more detailed selection.
        maxCliqueSize = 0
        for clique in list(completeSubgraphs[structNum]):
            #print(completeSubgraphs[clique])
            cliqueNodes = list(completeSubgraphs[clique])
            size = len(cliqueNodes)
            if size > setSize:
                #...what?
                #I don't _think_ it's possible, but the error message below is for the
                #case in which you find a complete set that is larger than expected.
                print("Interesting! Somehow a complete group of " + size +\
" structures was found before a complete group of your specified size.\n\
If you're reading this, let Kevin know and save your console input/output.\n\
The script should finish fine, but you might end up with more structures than\n\
you asked for. In this case, you can manually ignore one of the generated structures\n\
to get your specified set size.")
            if size > maxCliqueSize:
                maxClique = cliqueNodes
                maxCliqueSize = size
                if size >= setSize: #store node list and crash out
                    foundSet = True

    if verbose:
        print("Largest complete group so far: ", maxClique)
    return foundSet, maxClique

#Function to move relevant log and dihedral files to the current directory
#Params: 
#structs - At this point in the script, structs is a complete graph of size setSize
#whose nodes are the structures to be returned to user and whose edges are
#the RMSD values between these structures.
def saveSetFiles(newSeq, structs, completeSet):
    nodes = list(structs.nodes)
    for structNum in nodes:
        structName = "s" + str(structNum)
        shutil.copyfile("%s.dChimOut" % structName, structName + "DChim.log" )
        shutil.copyfile("%s.chimOut" % structName, structName + "Chim.log" )
        shutil.copyfile("%s.di" % structName, structName +".dihed")

        if completeSet:
            print("\n%s is one of your initial structures.\n\
It is written to %s_int_oneh.pdb." % (structName, structName)) 
        os.rename("validStructs/%s_int_noh.pdb" % structName, structName + "_int_oneh.pdb")

    if completeSet:           
        rmsFileName = newSeq + ".rmsd"
        rmsFile = open(rmsFileName, "w")
        print("")
        for u,v,rmsValue in structs.edges.data('rms'):
            rmsText = "RMSD between %s and %s is %f A\n" % ('s'+str(u), 's'+str(v), rmsValue)
            rmsFile.write(rmsText)
            print(rmsText)
        rmsFile.close()
       
#Move the final desired sequences and relevant log files to a new directory with the 
#name of the sequence you are creating.
#
#Params: directoryName - string: The name of the output directory, usually the name of the sequence
def organizeFiles(directoryName, debug, keepall):
    while os.path.exists(directoryName):
        replace = input("Your created sequence already has a directory made for it.\n\
Would you like to replace this directory, '" + directoryName + "'? \n\
If 'yes', the directory will be replaced. \n\
If 'no', you will be prompted for a name for your output directory. \n")
        if replace == "yes": 
            shutil.rmtree(directoryName) #breaks loop
        if replace == "no":
            directoryName = input("Please input a name for your output directory.\n")
    
    os.mkdir(directoryName)

    print("Moving your structure PDBs and log files to %s/" % directoryName)
            
    if (keepall or debug): #if keepall or debug are on, we don't want to delete anything
        os.system("mv *.chimOut *.dChimOut *.di %s/" % directoryName)
    else:
        os.system("rm *invalid.pdb *mut_min.pdb *inverted.pdb")
        os.system("rm *.dChimOut")
        os.system("rm *.chimOut")
        os.system("rm *.di")
    os.system("mv *.log *.pdb *.dihed *.rmsd %s/" % directoryName)
    os.system("mv validStructs %s/allValidStructs" % directoryName)

#Function to create the temporary structure directory, removing validStructs directories 
#if they already exist
def mkStructsDirectory():
    if os.path.exists("validStructs"):
        shutil.rmtree("validStructs")
    os.mkdir("validStructs")
