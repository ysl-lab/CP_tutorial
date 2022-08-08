#Aidan Fike
#July 21, 2019
#
#Updated by Kevin Schult March-July 2020 with contribution from Tim Ling
#File to create a chimera script which creates a random pdb structure

#NOTE: Errors that are produced by Chimera are sent to a log file, not to
#the console. Check everything and make sure that the logs make sense.

import os
import numpy as np

#File to create a chimera script to create a structure with amino acids in "seq", and with the
#passed phi/psi angles
def createChimScript(seq, phi, psi, bonds, structNum, debug):
    scriptName = "chimScript.py"
    chim = open(scriptName, "w") 
    chim.write("import chimera\n")
    chim.write("from chimera import runCommand\n\n")

    chim.write("runCommand('open code/gly.pdb')\n\n")

    chim.write("runCommand('rotation 1 :1@C :1@CA')\n")
    chim.write("runCommand('rotation 1 %s')\n" % phi[0])
    chim.write("runCommand('~select :1@C :1@CA')\n")

    for index in range(1,len(seq)):
        chim.write("runCommand('addaa gly,%s,%s,%s %s')\n" % (str(index+1), phi[index],\
                                                                psi[index],":"+str(index)+".a"))
        if bonds[index-1] == "c":
            chim.write("runCommand('rotation 2 %s %s')\n" % (":"+str(index)+".a@C",\
                                                             ":"+str(index+1)+".a@N"))
            chim.write("runCommand('rotation 2 180')\n")
            chim.write("runCommand('~rotation 2')\n")
            chim.write("runCommand('~select :.a')\n")

    if debug:
        #Debug: write linear glycine peptide
        chim.write("runCommand('write #0 s%d_linear.pdb ')\n" % structNum) 

    chim.write("runCommand('delete :%s@OXT')\n" % len(seq))
    chim.write("runCommand('bond :%s@C :1@N')\n" % len(seq))

    chim.write("runCommand('minimize nogui True nsteps 1000')\n")

    if debug:
        #Debug: write cyclic minimized glycine peptide
        chim.write("runCommand('write #0 s%d_cGly.pdb ')\n" % structNum) 

    for index, amino in enumerate(seq):
        toInvert = False
        isPro = False
        if amino.islower():
            toInvert = True
            amino = amino.upper()
        if amino == "P":
            amino = "A"
            isPro = True
        chim.write("runCommand('swapaa %s #:%s.a')\n" % (oneToThree(amino), index+1))
        if toInvert:
            chim.write("runCommand('invert :%s@ca')\n" % (index + 1))
            if amino == "V":
                chim.write("runCommand('invert :%s@cb')\n" % (index + 1))
        if isPro:
            chim.write("runCommand('swapaa pro #:%s.a')\n" % (index + 1))

    if debug:
        #Debug: write mutated, unminimized cyclic peptide
        chim.write("runCommand('write #0 s%d_mut_nomin.pdb ')\n" % structNum) 

    chim.write("runCommand('minimize nogui True nsteps 15000 cgsteps 1000')\n")

    for i in range(len(seq)):
        chim.write("runCommand('chirality :%d.A@CA')\n" % (i+1))

    chim.write("runCommand('select element.H')\n")
    chim.write("runCommand('~select :1.A@H')\n")
    chim.write("runCommand('delete sel')\n")

    chim.write("runCommand('write #0 s%d_mut_min.pdb ')\n" % structNum) 

    chim.close()

    return scriptName

#Function to invert non-proline residues with the incorrect chirality
def invertBadChiralities(structNum, newSeq, verbose):
    chimLog = open("s%d.chimOut" % structNum, "r")

    chirality = []
    length = len(newSeq)

    for line in chimLog:
        words = line.split()
        for i in range(length):
            if len(words) > 0 and words[0] == "#0:%d.A@CA" % (i+1):
                chirality.append(words[2])
    chimLog.close()

    scriptName = "dAminoScript.py"
    damino = open(scriptName, "w")

    damino.write("import chimera\n")
    damino.write("from chimera import runCommand\n\n")

    damino.write("runCommand('open s%d_mut_min.pdb')\n" % structNum)

    #Keeping H2T bond consistent
    damino.write("runCommand('bond :%s@C :1@N')\n" % str(length))

    errorFound = False
    for i in range(len(chirality)):
        if newSeq[i].upper() == 'P':
            #inversion doesn't work on proline for cyclic peptides.
            if verbose:
                print("Proline detected, chirality inversion not attempted.")
            continue #move to next residue
        amino = "L"
        if newSeq[i].islower():
            amino = "D"
        if amino == "L":
            if (chirality[i] == 'R' and newSeq[i].upper() != 'C') \
                    or (chirality[i] == 'S' and newSeq[i].upper() == 'C'):
                if verbose:
                    print("Incorrect D amino acid found at %d%s in struct %d. Now inverting it to an L amino acid" % (i+1,newSeq[i],structNum))
                damino.write("runCommand('invert :%d.A@CA')\n" % (i+1))
                errorFound = True

        if amino == "D":
            if (chirality[i] == 'S' and newSeq[i].upper() != 'C') \
                    or (chirality[i] == 'R' and newSeq[i].upper() == 'C'):
                if verbose:
                    print("Incorrect L amino acid found at %d%s in struct %d. Now inverting it to a D amino acid" % (i+1,newSeq[i],structNum))
                damino.write("runCommand('invert :%d.A@CA')\n" % (i+1))
                errorFound = True

    if errorFound:
        damino.write("runCommand('minimize nogui True, nsteps 1000')\n")

        damino.write("runCommand('select element.H')\n")
        damino.write("runCommand('~select :1.A@H')\n")
        damino.write("runCommand('delete sel')\n")

        damino.write("runCommand('write #0 s%d_mut_inverted.pdb ')\n" % structNum)
    elif verbose:
        print("No incorrect chiralities inverted in struct %d" % structNum)

    damino.close()

    if errorFound:
        # Planning to change path to chimera:
        os.system("/cluster/tufts/ylin12/tim/UCSF-Chimera64-1.14/bin/chimera --script " + scriptName + " --nogui &> s%d.dChimOut" % structNum)
    else:
        dChimOut = open("s%d.dChimOut" % structNum, "w")
        dChimOut.write("No chirality errors found by chimera")
        dChimOut.close()

    os.remove(scriptName)

    return errorFound


#Kevin's note: We're trying to figure out what this function does, and/or why we
#need it. If you need or can't have CONECT records, make changes accordingly.
#Despite what Aidan says, it probably doesn't have anything to do with dihedrals.
#Aidan's commentary:
#Insert into a pdb file the dihedral connections which are missed by chimera because of the 
#weird cyclic nature
def adjustStructFile(structNum, newSeq, verbose, debug):

    chiralityError = invertBadChiralities(structNum, newSeq, verbose)

    if chiralityError:
        oldPdb = open("s%d_mut_inverted.pdb" % structNum, "r")
    else:
        oldPdb = open("s%d_mut_min.pdb" % structNum, "r")
    newPdb = open("s%d_int_noh.pdb" % structNum, "w")

    N1 = -1
    CA1 = -1 
    H1 = -1
    CAx = -1
    Cx = -1
    Ox = -1

    #Find relevant atoms in the current pdb and write the new pdb
    for line in oldPdb:
        words = line.split()
        if words[0] != "END" and words[0] != "CONECT":
            if words[5] == "1":
                if words[2] == "N":
                    N1 = int(words[1])
                if words[2] == "CA":
                    CA1 = int(words[1])
                if words[2] == "H":
                    H1 = int(words[1])
            elif words[5] == str(len(newSeq)):
                if words[2] == "C":
                    Cx = int(words[1])
                if words[2] == "CA":
                    CAx = int(words[1])
                if words[2] == "O":
                    Ox = int(words[1])

            newPdb.write(line) 
        else:
            newPdb.write("CONECT %4d %4d %4d %4d\n" % (N1, CA1, Cx, H1))
            newPdb.write("CONECT %4d %4d %4d %4d\n" % (Cx, N1, CAx, Ox))
            newPdb.write("END\n")

    oldPdb.close()
    newPdb.close()
    if verbose:
        print("Successfully created s%d_int_noh.pdb" % structNum)
    
#Code to convert one letter amino acid codes to three-letter amino acid codes
def oneToThree(one):
    if one == "A":
        return "ala"
    if one == "C":
        return "cys"
    if one == "D":
        return "asp"
    if one == "E":
        return "glu"
    if one == "F":
        return "phe"
    if one == "G":
        return "gly"
    if one == "H":
        return "his"
    if one == "I":
        return "ile"
    if one == "K":
        return "lys"
    if one == "L":
        return "leu"
    if one == "M":
        return "met"
    if one == "N":
        return "asn"
    if one == "P":
        return "pro"
    if one == "Q":
        return "gln"
    if one == "R":
        return "arg"
    if one == "S":
        return "ser"
    if one == "T":
        return "thr"
    if one == "V":
        return "val"
    if one == "W":
        return "trp"
    if one == "Y":
        return "tyr"
