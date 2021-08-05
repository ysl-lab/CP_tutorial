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
  set sel [atomselect $molid "resid $resid and name $name"]
  set atomid [$sel get serial]
  return $atomid
}


proc genPhiPsiAtomList {molid resid0 resid1 resid2} {
  set atomidList [list]

  lappend atomidList [getAtomID $molid $resid0 "C"]
  lappend atomidList [getAtomID $molid $resid1 "N"]
  lappend atomidList [getAtomID $molid $resid1 "CA"]
  lappend atomidList [getAtomID $molid $resid1 "C"]
  lappend atomidList [getAtomID $molid $resid2 "N"]
  return $atomidList
}


proc writeAtomList {ofp atomidList} {
  lassign $atomidList a0 a1 a2 a3 a4 
  puts $ofp "$a0 $a1 $a2 $a3"
  puts $ofp "$a1 $a2 $a3 $a4"
}



#############################################
set numRes 6
set GRO "equil_npt.gro"

set molid [ mol new $GRO type gro ]

# Set up residue list based on CA atoms
set calpha [atomselect $molid "name CA"]
set resids [$calpha get resid]


set ofp [open "index.ndx" w]
puts $ofp {[ PhiPsi ]}

foreach resi $resids {
  if { $resi == 1} {
    set atomidList [ genPhiPsiAtomList $molid $numRes 1 2 ]
  } elseif { $resi == $numRes} {
    set atomidList [ genPhiPsiAtomList $molid [expr $numRes-1] $numRes 1 ]
  } else {
    set atomidList [ genPhiPsiAtomList $molid [expr $resi -1] $resi [expr $resi+1]]
  }

  writeAtomList $ofp $atomidList

}


close $ofp

exit 0
