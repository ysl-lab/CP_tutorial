#This script will create an index (.ndx) file that contains the atoms needed to calculate omega angles.
#Then the index file is input for the GROMACS command "gmx angle", which will calulate the omega angles 
#and output a .xvg file containing the omega angles for each residue for each frame in the trajectory.
#Finally, the script will read this .xvg file and check if any angle is below a threshold angle, such 
#that the angle would be classified as "cis".

import optparse
import numpy as np
import os

#This function searches the input .gro file for the atom numbers relevant to the omega dihedrals to be added.
def get_omega_atoms(gro):
    fi = open(gro)
    lines = fi.readlines()
    fi.close()
    #Reads the .gro file line by line. Looks at the second to last line for length of peptide.
    length = lines[-2][0:5].strip()
    #Atoms needed to define omega.
    target = ["CA", "C", "N","CH3"]
    mem = {}
    for line in lines:
        #Collect the residue id in the current line
        resID = line[0:5].strip()
        #collect the residue name in the current line
        resname = line[5:10].strip()
        #Collect the atom in the current line
        atomName = line[10:15].strip()
        #Ignore all atoms that are not a CA, C, N, H or O.
        if atomName not in target:
            continue
        #Determine the atom number for the current line, and relevant atom
        atomID = line[15:20].strip()
        #Store both the residue number and atom number in a dictionary to index later.
        entryID = atomName + resID
        mem[entryID] = atomID
    return length, mem

def get_caps(gro):
    N_cap=type(False)
    C_cap=type(False)
    fi = open(gro)
    lines = fi.readlines()
    fi.close()
    peptide_length = lines[-2][0:5].strip()
    target = ["ACE","NME"]
    first_residue=(lines[+2][5:10].strip())
    last_residue=(lines[-2][5:10].strip())
    if first_residue in target:
        N_cap=True
    if last_residue in target:
        C_cap=True
    return N_cap,C_cap

#A quick function to turn "True" into True the boolean and "False" to False the boolean.
def get_parity(arg):
    return True if arg[0].upper() == "T" else False

#A function to write the omega defintions into a .ndx file. 
def write_omega(omegas,filename):
    f= open(filename,"w+")
    f.write('[Omegas]\n')
    for i in range(len(omegas)):
     f.write(omegas[i]+'\n')
    f.close()
    return

#The main function collects the commandline arguments and uses the other functions in the program to determine the presence of cis and trans bonds.
if __name__ == "__main__":
    #Define commandline arguments/flags
    parser = optparse.OptionParser()
    parser.add_option('--gro', dest = 'gro')
    parser.add_option('--trj', dest = 'trj')
    parser.add_option('--cyclic', dest = 'cyclic', default = 'false')
    parser.add_option('--cutoff', dest = 'cutoff', default = 90)
    #Collect the arguments given by the user
    (options, args) = parser.parse_args()
    gro = options.gro
    trj = options.trj
    cutoff = options.cutoff
    cyclic = get_parity(options.cyclic)

length, atom_dic = get_omega_atoms(gro)
omegas = []
length= int(length)
#First check if the peptide is capped with ACE and/or NME.
N_cap,C_cap= get_caps(gro)

#For each residue, find its' atoms that define its' omega angle. 
#For example, residue i would need the atom indices corresponding to CAi, Ci, Ni+1, and CAi+2.
#These definitions are a little different for capped peptides.

#Linear capped case: 
if N_cap ==True:
    dihedral ="%5s %5s %5s %5s" %(atom_dic["CH3" + str(1)] ,atom_dic["C" + str(1)], \
    atom_dic["N" + str(2)],atom_dic["CA" + str(2) ])
    omegas.append(dihedral)

    for i in range(2,length-1):
        dihedral="%5s %5s %5s %5s" %(atom_dic["CA" + str(i)] ,atom_dic["C" + str(i)], \
        atom_dic["N" + str(i+1)],atom_dic["CA" + str(i+1) ])
        omegas.append(dihedral)
    if C_cap==True:
        dihedral="%5s %5s %5s %5s" %(atom_dic["CA" + str(length-1)] ,atom_dic["C" + str(length-1)], \
        atom_dic["N" + str(length)],atom_dic["CA" + str(length) ])
        omegas.append(dihedral)
    if C_cap==False:
        dihedral="%5s %5s %5s %5s" %(atom_dic["CA" + str(length-1)] ,atom_dic["C" + str(length-1)], \
        atom_dic["N" + str(length)],atom_dic["CA" + str(length)])
        omegas.append(dihedral)

if (N_cap ==False) and (C_cap==True):
    for i in range(1,length-1):
        dihedral="%5s %5s %5s %5s" %(atom_dic["CA" + str(i)] ,atom_dic["C" + str(i)], \
        atom_dic["N" + str(i+1)],atom_dic["CA" + str(i+1) ])
        omegas.append(dihedral)

    dihedral="%5s %5s %5s %5s" %(atom_dic["CA" + str(length-1)] ,atom_dic["C" + str(length-1)], \
    atom_dic["N" + str(length)],atom_dic["CA" + str(length) ])
    omegas.append(dihedral)


#Uncapped peptide cases (linear and cyclic):
if (N_cap == False) and (C_cap==False):
    for i in range(1,length):
      dihedral="%5s %5s %5s %5s" %(atom_dic["CA" + str(i)] ,atom_dic["C" + str(i)], \
            atom_dic["N" + str(i+1)],atom_dic["CA" + str(i+1) ])
      omegas.append(dihedral)

    #For the cyclic peptide case, we need to add the head_tail dihedral. 
    #For example, a linear pentapeptide has 4 omega angles, but a cyclic pentapeptide has 5 omega angles.
    if cyclic:
      head_tail_dihedral= "%5s %5s %5s %5s" %(atom_dic["CA" + str(length)] ,atom_dic["C" + str(length)], \
            atom_dic["N" + str(1)],atom_dic["CA" + str(1) ])
      omegas.append(head_tail_dihedral)

#Naming of input .ndx files. 
baseName = trj.split(".")
indexName = baseName[0]+"_omega.ndx"
print(baseName[0])
angleXVGFile = baseName[0]+"_omega.xvg"
angdistFile = baseName[0]+"_angdist.xvg"
write_omega(omegas,indexName)

#Execute GROMACS command using the generated index file. Output are 2 .xvg files. the XXX_omega.xvg contains the omega angles for each residue at each frame 
#and the XXX_angdist.xvg contains a distribution fo the omega angles (which is not what we need for cis/trans, but gromacs will output it by default).
os.system("gmx_mpi angle -f "+trj+" -n "+indexName+" -ov "+angleXVGFile+" -od "+angdistFile+" -all -xvg none -type dihedral") #this is where you need .trj .xtc
#print(list(range(2,length+2)))

#Read only the angles that we want from the output of gmx angle
if cyclic:
    Angles = np.loadtxt(angleXVGFile,usecols=list(range(2,length+2)))
else:
    Angles = np.loadtxt(angleXVGFile,usecols=list(range(2,length+1)))
#Define an array that holds the frame number and angle number that do not satisfy our criteria for cis/trans
BadAngles = []
#print(len(Angles[:,0]))
#Fix an issue with "not enough dimensions" if dealing with a gro file
if gro == trj:
    Angles = Angles[np.newaxis]
print(len(Angles[:,0]))
#Scan through every angle and frame in the collected omega angles array and check each one against the cutoff.
for i in range(len(Angles[:,0])):
    for j in range(length):
        if abs(Angles[i,j]) < int(cutoff):
            BadAngles.append([i,j])

#print(BadAngles)
#Ensure the user sees any problematic angles
for Angle in BadAngles:
    print("Angle", Angle[1], "in frame",Angle[0],"is cis")

#Calculate the percent of frames one or more cis bonds
print("Percent of frames with cis bonds is:", 100*float(len(set(list(i[0] for i in BadAngles)))/len(Angles[:,0])),"%")

#Now collect any frames that have at least one cis bond in CisAngles, and any frames that are not in CisAngles in TransAngles
TransIndex = []
CisIndex = list(set(list(i[0] for i in BadAngles)))
for i in range(1,len(Angles[:,0])+1):
    if i not in CisIndex:
        TransIndex.append(i)
#Write the frames to a .ndx file so that a user could filter their trajectory. This is only really helpful for .xtc analyses.
with open(baseName[0]+"_trans.ndx",'w') as outfile:
    outfile.write(" [frame] \n")
    for i in TransIndex:
        outfile.write(str(i)+"\n")
with open(baseName[0]+"_cis.ndx",'w') as outfile:
    outfile.write(" [frame] \n")
    for i in CisIndex:
        outfile.write(str(i)+"\n")
