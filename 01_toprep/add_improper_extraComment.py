#This program adds the two improper dihedrals in a cyclic peptide that were excluded during the process of creating a cyclic peptide topology with Gromacs according to lab protocol.
#It also adds one improper dihedral along each peptide bond in the cyclic peptide in order to help prevent cis bonds during the production simulation.
#To use this program, do: python --gro [.gro file] --ori [Topology_without_Imp.top] --out [Topology_with_imp.top]
#If not specified the program searches for the default filenames: "prot.gro", "cx_amber99sbMod_tip3p_temp.top", "cx_amber99sbMod_tip3p.top", respectively for the .gro file, the input topology, and the output topology.

import optparse
#This function searches the input .gro file for the atom numbers relevant to the dihedrals to be added.
def get_improper_atom(gro):
    fi = open(gro)
    lines = fi.readlines()
    fi.close()
    #Collect the number of atoms in the .gro file
    length = lines[-2][0:5].strip()
    
    target = ["CA", "C", "N", "H", "O"]
    mem = {}
    for line in lines:
        #Collect the residue id in the current line
        resID = line[0:5].strip()
        #Ignore all atoms in the middle residues of the .gro file
        if resID != "1" and resID != length:
            continue
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

#This function adds to the input topology improper dihedrals that are missing. The resulting text is written to the output topology file.
def add_improper(topOri, topOut, improper1, improper2,extraImps):
    #Read the input topology
    fi = open(topOri)
    lines = fi.readlines()
    fi.close()
    with open(topOut, "w+") as fo:
        for i in range(1, len(lines)):
            #Add the additional dihedrals
            if lines[i] == "; Include Position restraint file\n":
                fo.write(improper1)
                fo.write("\n")
                fo.write(improper2)
                fo.write("\n")
		for Imp in extraImps:
		    fo.write(Imp)
		    fo.write("\n")
            #Write the original topology to the new file
            fo.write(lines[i - 1])
        fo.write(lines[-1])
#This function scans the gro file for relevant atoms to the peptide bonds connecting each residue.
def get_extra_improper_atom(gro):
    #Read the input topology
    fi = open(gro)
    lines = fi.readlines()
    fi.close()
    #Collect the number of atoms in the .gro file
    length = lines[-2][0:5].strip()
    
    target = ["C", "N", "H", "O"]
    mem = {}
    
    for line in lines:
        #Collect the residue id in the current line
        resID = line[0:5].strip()
        #Collect the atom in the current line
        atomName = line[10:15].strip()
            #Ignore all atoms that are not C, N, H or O.
        if atomName not in target:
            continue
        #Determine the atom number for the current line, and relevant atom
        atomID = line[15:20].strip()
        #Store both the residue number and atom number in a dictionary to index later.
        entryID = atomName + resID
        mem[entryID] = atomID
    return length, mem
    
#The main function collects the commandline arguments and uses the other functions in the program to determine the correct improper dihedrals' definitions.
if __name__ == "__main__":
    #Define commandline arguments/flags
    parser = optparse.OptionParser()
    parser.add_option('--gro', dest = 'gro', default = 'prot.gro')
    parser.add_option('--ori', dest = 'topOri', default = 'cx_amber99sbMod_tip3p_temp.top')
    parser.add_option('--out', dest = 'topOut', default = 'cx_amber99sbMod_tip3p.top')
    #Collect the arguments given by the user
    (options, args) = parser.parse_args()
    gro = options.gro
    topOri = options.topOri
    topOut = options.topOut

    #Determine the atom numbers that define each of the target atoms in the .gro file
    length, improper_dic = get_improper_atom(gro)
    #Use the atom names as a key for the dictionary storing the atom numbers for the two improper dihedrals to be added.
    improper1 = "%5s %5s %5s %5s %5s     ; added by MD" % (improper_dic["CA" + length],
                                                           improper_dic["N1"],
                                                           improper_dic["C" + length],
                                                           improper_dic["O" + length],
                                                           4)
    improper2 = "%5s %5s %5s %5s %5s     ; added by MD" % (improper_dic["C" + length],
                                                           improper_dic["CA1"],
                                                           improper_dic["N1"],
                                                           improper_dic["H1"],
                                                           4)
    
    #Add one improper dihedral along each peptide bond in the protein
    extraImps = []
    extraLength,extra_dic = get_extra_improper_atom(gro)
    #Do the 1 to n-1th dihedrals
    for i in range(1,int(length)):
        currentImp = "%5s %5s %5s %5s %5s     ; added by MD" % (extra_dic["H" + str(i+1)],
                                                           extra_dic["N"+str(i+1)],
                                                           extra_dic["C" + str(i)],
                                                           extra_dic["O" + str(i)],
                                                           4)
        extraImps.append(currentImp)
    #Do the nth dihedral
    currentImp = "%5s %5s %5s %5s %5s     ; added by MD" % (extra_dic["H" + "1"],
                                                           extra_dic["N"+"1"],
                                                           extra_dic["C" + str(length)],
                                                           extra_dic["O" + str(length)],
                                                           4)
    extraImps.append(currentImp)
    
    #Add the dihedrals to the output topology based on the input topology
    add_improper(topOri, topOut, improper1, improper2,extraImps)
