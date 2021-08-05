#!/usr/bin/env bash
# the next line restart using vmd \
exec vmd "-dispdev text -e  $0"

################################################################################
#                                                                              #
# Generate the index file for calculating the Phi/Psi angles of the residue    #
# by using Gromacs g_dihedral programs.                                        #
#                                                                              #
# This script is designed for cyclic peptide, since the g_rama fails to work   #
# for the calculation of N-terminal Phi angle and C-terminal Psi angle.        #
#                                                                              #
################################################################################


proc getAtomID {molid resid name} {
  set sel [atomselect $molid "resid $resid and $name"]
  set atomid [$sel get serial]
#  set atomid [$sel get index]
  return $atomid
}

#############################################
set numRes 6

set GRO "solute.gro"
set molid [ mol new $GRO type gro ]

#set PDB "prot.pdb"
#set molid [mol new $PDB type pdb]

# Set up residue list based on CA atoms
set calpha [atomselect $molid "name CA"]
set resids [$calpha get resid]


set ofp [open "IntGroups.ndx" w]

foreach resi $resids {
  puts $ofp "\[ NH$resi \]"
  puts $ofp " [getAtomID $molid $resi "name N H"]"
  
  puts $ofp "\[ CO$resi \]"
  puts $ofp " [getAtomID $molid $resi "name C O"]"
}

close $ofp

exit 0
