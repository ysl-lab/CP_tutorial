#Aidan Fike
#July 21, 2019
#
#Edited and updated to version 2
#Kevin Schult, March/April 2020
#
#Edit and revised
#Minh N. Ho, August 2022
#
#Main function. Able to create initial sequences for a given sequence passed by the user.
#See README for more information

import sys
sys.path.insert(1, 'code/') #allows for the files in code/ to be imported and used here.

import chooseStruct 
import parseCommandLine 
import dihedDir 
import chimScriptMaker 
import os
import shutil
import networkx as nx

# NOTE: make sure to specify a binary executable to chimera
# For example: chimeraBinaryExecPath=/cluster/tufts/ylin12/tim/UCSF-Chimera64-1.14/bin/chimera"
chimeraBinaryExecPath="/path/to/chimera_X.Y.Z/bin/chimera"

def main():
    #Parse command line arguments/obtain the sequence whose structures will be made
    newSeq, thresh, setSize, bonds, maxTries, includeOxygen, verbose, debug = parseCommandLine.getArgvInfo()

    #Mk directories
    chooseStruct.mkStructsDirectory()

    structs = nx.Graph()
    foundSet = False

    #For maxTries iterations, create a structure with random phi psi angles, then only keep it if it
    #does not have any incorrect peptide bonds
    structNum = 0
    while not foundSet:
        structNum += 1
        
        if structNum > maxTries:
            print("%d structures have been created and evaluated for rmsd/cis bonds but no complete set of\n\
structures with a backbone rmsd greater than %f A were found. Try rerunning the program \n\
with a lower backbone rmsd threshold value or more attempts for better success (Remember, the rmsd should \n\
be given in angstroms, not [nm]).\n" % (maxTries, thresh))
            response = ""
            while response != "continue" and response != "keepall" and response != "end":
                response = input("If you wish to continue and generate an additional " + str(maxTries) + " structures, enter 'continue'.\n\
If you wish to end and keep all valid and invalid structures, enter 'keepall'.\n\
If you wish to end and delete all structures, enter 'end'.\n")
            if response == "continue":
                maxTries = int(2*maxTries)
                print("Resuming...\n")
            elif response == "keepall":
                keepall = True
                chooseStruct.saveSetFiles(newSeq, structs, foundSet)
                chooseStruct.organizeFiles(newSeq + "_terminated", debug, keepall)
                sys.exit()
            elif response == "end":
                os.system("rm *.dChimOut") 
                os.system("rm *.chimOut") 
                os.system("rm *.di") 
                shutil.rmtree("validStructs")
                sys.exit()

        phi, psi = dihedDir.createRandomPhiPsi(1, newSeq)
        dihedDir.writeToDihedFile(newSeq, phi[0], psi[0], structNum)
        scriptName = chimScriptMaker.createChimScript(newSeq, phi[0], psi[0], bonds, structNum, debug)
        os.system(chimeraBinaryExecPath + " --script " + scriptName + " --nogui &> s%d.chimOut" % structNum)
        os.remove(scriptName)
        chimScriptMaker.adjustStructFile(structNum, newSeq, verbose, debug)

        validPeptide = chooseStruct.checkPeptide(structNum, newSeq, bonds, verbose)
        if not validPeptide:
            continue #restart loop without adding structure to structs
        elif structs.number_of_nodes() == 0:
            structs.add_node(structNum)
        else:
            structs = chooseStruct.addToGraph(thresh, structNum, structs, includeOxygen)
            #Check subgraphs containing structNum - might be one that's large enough
            foundSet, maxClique = chooseStruct.checkSubGraphs(structNum, structs, setSize, verbose)

        if foundSet:
            #Loop's about to end - prep for data reporting and file organization
            for node in list(structs.nodes):
                if not node in maxClique:
                    structs.remove_node(node)
        
    #Record that this structure has been generated
    seqFile = open("finSeqs.txt", "a+")
    seqFile.write(newSeq + "\n")
    seqFile.close()

    #Move the final structures and log files to a new directory named "newSeq"
    #and clean up irrelevant log files
    chooseStruct.saveSetFiles(newSeq, structs, foundSet)
    chooseStruct.organizeFiles(newSeq, debug, False)

    #os.chdir("..")
    #prefix = input("What do you want as your numeral prefix")
    #os.mkdir(prefix + "_" + newSeq)

main()
