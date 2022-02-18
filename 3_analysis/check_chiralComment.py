#This script will create an index (.ndx) file that contains the atoms needed to calculate improper angles for chirality.
#Then the index file is input for the GROMACS command "gmx angle", which will calulate the improper angles 
#and output a .xvg file containing the angles for each residue for each frame in the trajectory.
#Finally, the script will read this .xvg file and check if any angle is below a threshold angle, such 
#that the angle would be classified as "cis".
#Feb 2022, This script runs after loading Gromacs and required packages (gcc/4.9.2, openmpi/2.1.2, and cuda/8.0.44, as well as loading "module load python/3.6.0", I have found issues with "module load python", which loads python 2.7.3 or more directly "module load python/2.7.3", which appears to disrupt Gromacs.

import optparse
import numpy as np
import os
#This function searches the input .gro file for the atom numbers relevant to the improper dihedrals to be added.
def get_improper_atoms(gro):
    fi = open(gro)
    lines = fi.readlines()
    fi.close()
    #Reads the .gro file line by line. Looks at the second to last line for length of peptide.
    length = lines[-2][0:5].strip()
    target = ["CA", "C", "N","CB","CH3"]
    mem = {}
    for line in lines:
        #Collect the residue id in the current line
        resID = line[0:5].strip()
        #collect the residue name in the current line
        resname = line[5:10].strip()
        #Collect the atom in the current line
        atomName = line[10:15].strip()
        #Ignore all atoms that are not a CA, C, N, or CB.
        if atomName not in target:
            continue
        #Determine the atom number for the current line, and relevant atom
        atomID = line[15:20].strip()
        #Store both the residue number and atom number in a dictionary to index later.
        entryID = atomName + resID
        mem[entryID] = atomID
    return length, mem

#Identify if the peptide has glycines
def get_glycines(gro):
    glycine_indicies=[]
    fi = open(gro)
    lines = fi.readlines()
    fi.close()
    peptide_length = lines[-2][0:5].strip()
    Glycine_checker=np.zeros(int(peptide_length))
    target = ["GLY"]
    for line in lines:
        #Collect the residue id in the current line
        resID = line[0:5].strip()
        #collect the residue name in the current line
        resname = line[5:10].strip()
        #Ignore all atoms that are not a CA, C, N, H or O.
        if resname in target:
            glycine_indicies.append(resID)
    #Determine the unique indicies
    unique_indices=set(glycine_indicies)
    #Create a sorted list of the unique indices and make them integers
    sorted_indices = sorted([int(i) for i in list (unique_indices)])
    #Create a vector that indicates where the glycines are located. 1='L' and 0='D'
    for i in sorted_indices:
        Glycine_checker[i-1] = 1
    for i in list(unique_indices):
        print("Residue",i,"is glycine")
    return Glycine_checker

#Identify if the peptide is capped
def get_caps(gro):
    fi = open(gro)
    lines = fi.readlines()
    fi.close()
    #Collect the length of the peptide, initialize a indicator of caps
    peptide_length = lines[-2][0:5].strip()
    cap_checker=np.zeros(int(peptide_length))

    N_cap=type(False)
    C_cap=type(False)
    fi = open(gro)
    lines = fi.readlines()
    fi.close()
    peptide_length = lines[-2][0:5].strip()
    target = ["ACE","NME"]
    #We focus on only the first and last residues when looking for caps
    first_residue=(lines[+2][5:10].strip())
    last_residue=(lines[-2][5:10].strip())
    if first_residue in target:
        cap_checker[0]=1
        print("Residue 1 is either Ace or NMe")
    if last_residue in target:
        cap_checker[int(peptide_length)-1]=1
        print("Residue",peptide_length,"is either Ace or NMe")
    return cap_checker

#A function to write the improper angle defintions into a .ndx file.
def write_impropers(impropers,filename):#,index):
    f= open(filename,"w+")
    f.write('[Impropers]\n')
    for i in range(len(impropers)):
     f.write(impropers[i]+'\n')
    f.close()
    return

#A quick function to turn "True" into True the boolean and "False" to False the boolean.
def get_parity(arg):
    return True if arg[0].upper() == "T" else False

#The main function collects the commandline arguments and uses the other functions in the program to determine the correct improper dihedrals' definitions.
if __name__ == "__main__":
    #Define commandline arguments/flags
    parser = optparse.OptionParser()
    parser.add_option('--gro', dest = 'gro')
    parser.add_option('--trj', dest = 'trj')
    parser.add_option('--seq', dest = 'seq', default = '')
  
    #Collect the arguments given by the user
    (options, args) = parser.parse_args()
    gro = options.gro
    trj = options.trj
    seq = options.seq

    #Collect sequence if --seq flag not given
    while not seq:
        seq = input("Enter sequence:\n")

    #Translate sequence to vector where L=1 and D=0. L is a uppercase letter in user-inputed sequence and D is a lowercase letter in user-input
    test_vec=np.zeros(len(seq))
    for i in range(len(seq)):
        if seq[i].isupper():
            test_vec[i]=1

    #Determine the atom numbers that define each of the target atoms in the .gro file
    length, atom_dic = get_improper_atoms(gro)

    #Identify if there are glycines or capped residues present
    chirality_gly= get_glycines(gro)
    chirality_caps= get_caps(gro)

    #Create a vector that has 1="L" set for all Glycines and Caps
    chirality=np.zeros(int(length))
    for i in range(int(length)):
        if chirality_gly[i] ==1 or chirality_caps[i] ==1:
            chirality[i] =1
    print(chirality)
#For each residue, find its' atoms that define its' improper related to chirality. 
#For example, residue i would need the atom indices corresponding to CAi, Ni, Ci, and CBi.
#These definitions don't exist for Glycine and capped residues that don't have beta C.
#To deal with this, we create a dummy defintion, and those angles get's ignored in analysis/ 
#automatically set to L-form chirality.
    impropers = []
    length= int(length)
    for i in range(1,length+1):
        if chirality_gly[i-1]==1:
            dihedral="%5s %5s %5s %5s" %(atom_dic["CA" + str(i)] ,atom_dic["N" + str(i)], \
                atom_dic["C" + str(i)],atom_dic["CA" + str(i) ])
        elif chirality_caps[i-1]==1:
            dihedral="%5s %5s %5s %5s" %(atom_dic["CH3" + str(i)] ,atom_dic["CH3" + str(i)], \
                atom_dic["CH3" + str(i)],atom_dic["CH3" + str(i) ])
        else:       
            dihedral="%5s %5s %5s %5s" %(atom_dic["CA" + str(i)] ,atom_dic["N" + str(i)], \
                atom_dic["C" + str(i)],atom_dic["CB" + str(i) ])
        impropers.append(dihedral)

    #Naming of input .ndx files. 
    baseName = trj.split(".")
    indexName = baseName[0]+"_impropers_chirality.ndx"
    print(baseName[0])
    angleXVGFile = baseName[0]+"_impropers_chirality.xvg"
    angdistFile = baseName[0]+"_angdist_chirality.xvg"
    write_impropers(impropers,indexName)

    #Execute GROMACS command using the generated index file. Output are 2 .xvg files. the XXX_omega.xvg contains the improper angles for each residue at each frame 
    #and the XXX_angdist.xvg contains a distribution fo the angles (which is not what we need for chirality, but gromacs will output it by default).
    os.system("gmx_mpi angle -f "+trj+" -n "+indexName+" -ov "+angleXVGFile+" -od "+angdistFile+" -all -xvg none -type dihedral  &> "+baseName[0]+"_GromacsChiral.log")
    print(list(range(2,length+2)))
    Angles = np.loadtxt(angleXVGFile,usecols=list(range(2,length+2)))
    #NOTE THAT BADANGLES IS ZERO INDEXED
    BadAngles = []
    #This adds an extra empty dimension if you have only the gro file to analyze.
    if gro == trj:
        Angles = Angles[np.newaxis]
    print(len(Angles[:,0]))
    #Scan through each angle in every frame for any that do not match the correct chirality.
    for i in range(len(Angles[:,0])):
        for j in range(length):
            #Here we have a residue that is achiral
            if chirality[j] == 1:
                pass
            #For L chirality
            elif test_vec[j] == 1:
                if Angles[i,j] < 0:
                    BadAngles.append([i,j])
                if abs(Angles[i,j]) == 0:
                    print("Error: Inconclusive Chirality for frame:",str(i),", residue:",str(j+1), ", angle:",Angles[i,j],".")
                    
            #For D chirality
            else:
                if Angles[i,j] > 0:
                    BadAngles.append([i,j])
                if abs(Angles[i,j]) == 0:
                    print("Error: Inconclusive Chirality for frame:",str(i),", residue:",str(j+1), ", angle:",Angles[i,j],".")
            
    #print(BadAngles)
    #Show the user which frames and residues have bad chirality
    for Angle in BadAngles:
        print("Residue", Angle[1]+1, "in frame",Angle[0],"has incorrect chirality")

    #Calculate the percent of frames where at least one residue has an incorrect chirality.
    print("Percent of frames with bad chirality is ", 100*float(len(set(list(i[0] for i in BadAngles)))/len(Angles[:,0])))
    #Separate the frames into those with all good chiralities, and those with at least one bad chirality.
    GoodIndex = []
    BadIndex = list(set(list(i[0] for i in BadAngles)))
    for i in range(1,len(Angles[:,0])+1):
        if i not in BadIndex:
            GoodIndex.append(i)
    #Save to an index file the frames that are not problematic.
    with open(baseName[0]+"_good_chirality.ndx",'w') as outfile:
        outfile.write(" [frame] \n")
        for i in GoodIndex:
            outfile.write(str(i)+"\n")
    #Save to an index file the frames that have at least one bad chirality.
    with open(baseName[0]+"_bad_chirality.ndx",'w') as outfile:
        outfile.write(" [frame] \n")
        for i in BadIndex:
            outfile.write(str(i)+"\n")

    print("The analyzed sequence had", sum(chirality_gly), "glycines, and",sum(chirality_caps),"caps identified") 
