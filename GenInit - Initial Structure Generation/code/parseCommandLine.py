#Aidan Fike
#July 21, 2019

#File with function to parse command line arguments and extract relevant information. Also
#informs the user if they make a mistake

import sys
import os
import dihedDir

usage = "\nUSAGE: To use this program run:\n\n\
\n\
    python main.py [-flag argument] [-flag argument] [-flags]\n\
\n\
    Where [-flag] is one of the below:\n\
\n\
        -seq | Takes an argument of a string of upper- and lowercase letters to be\n\
                the cyclic peptide sequence.\n\
\n\
        -randseq | Takes an integer argument. Instructs the program to randomly\n\
                generate a cyclic peptide sequence of that length.\n\
                Only (G, A, V, F, S, N, R, D) will be used for random sequences.\n\
\n\
        -n | Takes an integer argument that indicates the number of distinct initial\n\
                structures to generate. Set to 2 by default.\n\
\n\
        -bonds | Takes an argument of a string of 't's and 'c's that specify whether\n\
                the peptide bond at the corresponding location should be cis or trans.\n\
                See examples. All trans by default.\n\
\n\
        -thresh | Takes a floating point argument. Structures will be generated\n\
                until enough of them are -thresh angstroms different in RMSD.\n\
\n\
        -maxTries | Takes an integer argument. If -maxTries structures have been\n\
                generated without creating an acceptable output set, the script\n\
                will terminate.\n\
\n\
        -includeOxygen | Takes no argument. By default, this script\n\
                calculates RMSD between the (C,CA,N) atoms of each structure.\n\
                This optional flag will tell the script to use (C,CA,N,O) atoms\n\
                in calculating RMSD.\n\
\n\
        -h | Takes no argument. Displays help text, then terminates the script.\n\
\n\
        -v | Takes no argument. Tells the script to produce verbose output.\n\
\n\
        -debug | Takes no argument. Tells the script to produce all intermediate structures.\n\
\n\
    Note that the script requires -thresh and one of -seq or -randseq to run.\n\
    -maxTries is 20 by default. -n is 2 by default. When -bonds is not set, all\n\
    peptide bonds are trans.\n\
\n\
    A couple examples are:\n\
        For a random 8mer, generate 2 structures with a minimum rmsd of 1.9A \n\
        that includes oxygen in rmsd calculation:\n\
            python main.py -randseq 8 -thresh 1.9 -includeOxygen\n\
        For a desired sequence, GNSRVGGGGG, generate 4 strucutres in 30 attempts with a minimum rmsd of 3.0A:\n\
            python main.py -seq GNSRVGGGGG -n 4 -thresh 3.0 -maxTries 30\n\
        For a desired sequence, EasESLCGG with a cis bond between the 3rd and 4th residues,\n\
        generate 3 structures in 50 attempts with a minimum rmsd of 2.1A:\n\
            python main.py -seq EasESLCGG -bonds ttctttttt -n 3 -thresh 2.1 -maxTries 50\n\n"





print("\nNOTE: To successfully run this program, you will need to \n\
load a version of python 3 which includes MDAnalysis. \n\
See the README for detailed instructions on loading the modules.\n\
\n\
NOTE: This version of the script takes RMSD threshold in angstroms, not nm. \n")

#Function to parse through all of the command line arguments. 
#
#Return: newSeq: string - The sequence whose structures will be created in this program
#        thresh: float - The backbone rmsd threshold
def getArgvInfo():
    newSeq = ""
    includeOxygen = False
    thresh = -1.0
    maxTries = 20
    setSize = 2 #number of distinct structures to try for
    bonds = ""
    verbose = False
    debug = False

    argumentsSpecified = []
    flags = ['-h','-debug','-v','-includeOxygen','-randseq','-seq','-thresh','-maxTries','-n','-bonds']
    flagsWithValues = ['-randseq','-seq','-thresh','-maxTries','-n','-bonds']

    for i in range(1,len(sys.argv)): #sys.argv is a list of the arguments.
                                     #The 0th argument is the name of the program, main.py.

        argI = sys.argv[i]

        #Check that the argument is a flag or the value for a processed flag.
        if not (argI in flags): 
            if sys.argv[i-1] in flagsWithValues: #argI is a value for a valid flag. Move along.
                continue
            else: #We have no idea what the user is doing. Stop it now.
                exitWithError("Argument '" + argI + "' is not recognized.")
        
        #Check for repeated arguments
        if argI in argumentsSpecified: #User has used the same flag twice.
            exitWithError("You have used the flag " + argI + " too many times.")
        else: #Mark this flag to not be repeated.
            argumentsSpecified.append(argI)

        #Check for arguments that don't have an additional value
        if argI == '-h': #User is requesting the script usage readme.
            print(usage, "\n")
            sys.exit()
        elif argI == '-includeOxygen': #User is asking for C,CA,N,O backbone definition for rmsd.
            includeOxygen = True
        elif argI == '-v': #User is asking for verbose output.
            verbose = True
        elif argI == '-debug':
            debug = True
            verbose = True #if they're debugging, they will also want verbose console output.

        #Check for arguments that preceed an input value
        if i != len(sys.argv): 
            
            if argI == '-randseq': #User wants random sequence, use following integer to generate one
                numAminos = sys.argv[i+1]
                if not str.isdigit(numAminos):
                    exitWithError("You didn't correctly specify the number of amino acids \
you wanted in your new random sequence.")
                numAminos = int(numAminos)
                newSeq = dihedDir.createNewSeq(numAminos)

            elif argI == '-seq': #User wants sequence as defined in following argument
                newSeq = sys.argv[i+1]
                if len(newSeq) == 1:
                    exitWithError("This program does not work for sequences of only 1 amino acid.")
                if not (newSeq.isalpha() and dihedDir.isAllowedSeq(sys.argv[i+1])):
                    exitWithError("The sequence you wanted to create is not possible with this program. \n\
Next time, please input a sequence with letters representing biological amino acids or their D-variants.")

            elif argI == '-thresh': #User is specifying threshold rmsd as following float argument.
                thresh = sys.argv[i+1]
                if not isFloat(thresh):
                    exitWithError("'%s' is not a valid decimal-formatted number. \n\
Next time you run, this argument should be used to input the threshold rmsd value you want. \n\
This will become the minimium backbone-aligned backbone rmsd between your two structures." % thresh)

            elif argI == '-maxTries': #User is manually setting the maximum number of attempted structures.
                maxTries = sys.argv[i+1]
                if not str.isdigit(maxTries):
                    exitWithError("You didn't correctly specify the maximum number of attempted structures.")
                maxTries = int(maxTries)

            elif argI == '-bonds': #User wants peptide bond omega angles as defined in following argument
                bonds = sys.argv[i+1]
                if not isBondsString(bonds):
                    exitWithError("The peptide bond orientation wasn't correctly specified.")

            elif argI == '-n': #User is manually setting the number of distinct structures needed.
                setSize = sys.argv[i+1]
                if not str.isdigit(setSize):
                    exitWithError("You didn't correctly specify the number of needed structures.")
                setSize = int(setSize)






    #All inputs gathered. Make sure that they make sense.

    if maxTries < setSize:
        exitWithError("maxTries must be greater than or equal to the size of the set you wish to produce.")

    if newSeq == "" or thresh == -1.0:
        exitWithError("You need to specify both your sequence (-seq or -randseq) \n\
and your rmsd threshold in angstroms (-thresh). One or both of these were \n\
not correctly defined.")

    #generate a -bonds value if none exists - all trans by default
    if bonds == "":
        for i in range(len(newSeq)):
            bonds = bonds + "t"
    else:
        print("NOTE: You have manually specified your peptide bond omega angles as cis or trans.\n\
Be aware that peptides with cis bonds may take more attempts to generate, \n\
as trans bonds are generated more often. You may need to consider whether your \n\
structure can be reasonably generated from random phi-psi angles and energy minimization.\n\
The verbose flag -v may be used to examine the pattern of cis bond generation in your sequence.\n")


    if len(newSeq) != len(bonds):
        exitWithError("Your sequence and your bond specification don't match in length.")

    #debugging
    #print(newSeq, float(thresh), setSize)
    return newSeq, float(thresh), setSize, bonds, maxTries, includeOxygen, verbose, debug

#Check if the sequence has already been specified.
def checkSeqDuplicate(seqSpecified):
    if seqSpecified:
        exitWithError("Please only input one -seq or -randseq flag.")
    return

#Check if a given user input is a float
def isFloat(x):
    try:
        float(x)
        return True
    except:
        return False

#Check if the "bonds" user input meets requirements
def isBondsString(bonds):
    try:
        if not (bonds.isalpha() and bonds.islower()):
        #    print("error A")
            return False
        for letter in bonds:
            if letter != "t" and letter != "c":
        #        print("error B, " + letter)
                return False
        return True
    except:
        #print("error C")
        return False

def exitWithError(error):
    print(usage, "\n")
    print("ERROR: " + error + "\n")    
    sys.exit()
