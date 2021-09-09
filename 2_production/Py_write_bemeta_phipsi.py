"""
python Py_write_bemeta_phipsi.py prot.gro

Created by Diana Slough
Created on 06-26-2016
Modified on 10-3-2016

This program writes out 2D-BEMETA files that use phi [C N CA C] and psi [N CA C N]
NOTE: assumes CP and currently supports all residues



"""

import sys
from get_dihedrals import get_psi, get_phi

NNEUTRAL = 5    # number of neutral replicas


###############################################################
# read_lines(f)
#
# Loads all lines of file and returns list containing the lines
#
# Input:
#   f: input file that has already been opened
###############################################################
def read_lines(f):
    lines = []
    for l in f: 
        lines.append(l)
    return lines


###############################################################
# read_gro(lines)
#
# Takes the lines of a file and reads the gromacs data
#
# Input:
#   lines: list containing all lines of the file
################################################################
def read_gro(lines):
    num_res = 0
    resids = []         # Residue ID numbers
    residues = []       # residue name
    atoms = []          # atom name
    indexes = []        # atom index

    # Loop over all lines of the file
    for l in range(len(lines)):
        if l != 0 and l != 1:               # ignore the first two lines
            if l < (len(lines) - 1):        
                if int(lines[l].strip().split()[0][0]) > num_res: 
                    num_res = int(lines[l].strip().split()[0][0])

                # get the residue information
                resids.append(int(lines[l].strip().split()[0][0])) 
                residues.append(str(lines[l].strip().split()[0][1:]))
                atoms.append(str(lines[l].strip().split()[1]))
                indexes.append(int(lines[l].strip().split()[2]))
            else:  # box line
                break

    # Return data
    return num_res, resids, residues, atoms, indexes


################################################################
# write_bemeta_files(num_res, phi, psi)
#
# Writes the information to a bemeta file
#
# Input:
#   num_res: number of residues in this structure
#   phi: phi indexes
#   psi: psi indexes
################################################################
def write_bemeta_files(num_res, phi, psi):

    num_replicas = num_res * 2 + NNEUTRAL

    # Open the file
    with open('bemeta.dat', 'w') as f:
        f.write('RANDOM_EXCHANGES\n\n')
        counter = 0
        rep_counter = 0

        # Print phi[i] - psi[i] data for each residue
        for j in range(num_res):
            f.write('#Rep%d; Res%d:phi/psi\n' % (rep_counter, j + 1))
            f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n' % (counter, phi[j][0], phi[j][1], phi[j][2], phi[j][3]))
            f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n\n' % (counter + 1, psi[j][0], psi[j][1], psi[j][2], psi[j][3]))
            counter += 2
            rep_counter += 1 
        
        # Print psi[i] - phi[i+] data
        for j in range(num_res):
            if j != num_res - 1:
                f.write('#Rep%d; Res%d:psi/Res%d:phi\n' % (rep_counter, j + 1, j + 2))
                f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n' % (counter, psi[j][0], psi[j][1], psi[j][2], psi[j][3]))
                f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n\n' % (counter + 1, phi[j + 1][0], phi[j + 1][1], phi[j + 1][2], phi[j + 1][3]))
            else:
                f.write('#Rep%d; Res%d:psi/Res%d:phi\n' % (rep_counter, j + 1, 1))
                f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n' % (counter, psi[j][0], psi[j][1], psi[j][2], psi[j][3]))
                f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n\n' % (counter + 1, phi[0][0], phi[0][1], phi[0][2], phi[0][3]))
            counter += 2
            rep_counter += 1 
        
        # Print the rest of the bemeta.dat
        cv_str = "{"
        height_str = "{"
        for i in range(num_res):
            cv_str += '{cv%d,cv%d},' % (i * 4, i * 4 + 1)
            cv_str += '{cv%d,cv%d},' % (i * 4 + 2, i * 4 + 3)
            height_str += '0.1,0.1,'
        for i in range(NNEUTRAL):
            cv_str += '{cv0,cv1},'
            height_str += '0.0,'
        cv_str = cv_str[:-1] + "}"
        height_str = height_str[:-1] + "}"
        f.write('metad: METAD ...\n')
        f.write('   ARG=@replicas:%s\n' % cv_str)
        f.write('   HEIGHT=@replicas:%s\n' % height_str)
        f.write('   SIGMA=0.31416,0.31416\n')
        f.write('   PACE=2000\n')
        f.write('   FILE=HILLS\n')
        f.write('...\n\n')
        f.write('PRINT ARG=@replicas:%s STRIDE=500 FILE=COLVAR\n\n' % cv_str) 
        f.write('ENDPLUMED\n')
    return


#########################################################################
# main()
# loads *.gro file, gets phi and psi indexes and saves to bemeta file
#########################################################################
def main():
    with open(sys.argv[1], 'r') as f:
        lines = read_lines(f)
   
    # Read the gromacs file
    num_res, resids, residues, atoms, indexes = read_gro(lines)
   
    # Get phi and psi indexes
    phi = get_phi(num_res, resids, residues, atoms, indexes)
    psi = get_psi(num_res, resids, residues, atoms, indexes)

    # Write to file
    write_bemeta_files(num_res, phi, psi)


if __name__ == '__main__':
    main()
