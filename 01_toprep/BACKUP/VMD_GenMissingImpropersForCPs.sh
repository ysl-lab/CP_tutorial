#!/usr/bin/env bash
# the next line restart using vmd \
exec vmd "-dispdev text -e  $0"


set NRES 6
set IRES 1

################################################################################
#                                                                              #
# Generate the index file for calculating the Phi/Psi angles of the residue    #
# by using Gromacs g_dihedral programs.                                        #
#                                                                              #
# This script is designed for cyclic peptide, since the g_rama fails to work   #
# for the calculation of N-terminal Phi angle and C-terminal Psi angle.        #
#                                                                              #
################################################################################


proc getNextResID { resi } {
  global NRES
  global IRES
  if { $resi == $NRES } {
    return $IRES 
  } else {
    return [expr $resi + 1 ] 
  }
}

proc getPreviousResID {resi } {
  global NRES
  global IRES
  if { $resi == $IRES } {
    return $NRES
  } else {
    return [expr $resi - 1 ] 
  }
}


proc getAtomID {molid resid name} {
  set sel [atomselect $molid "resid $resid and name $name"]
  set atomid [$sel get serial]
#  set atomid [$sel get index]
  return $atomid
}


proc getOmega1Atoms {molid resi} {
  set resi0 [getPreviousResID $resi ]
  set atomidList [list]
  lappend atomidList [getAtomID $molid $resi0 "C"]
  lappend atomidList [getAtomID $molid $resi  "CA"]
  lappend atomidList [getAtomID $molid $resi "N"]
  lappend atomidList [getAtomID $molid $resi "H"]
  return $atomidList
}


proc getOmega2Atoms {molid resi} {
  global NRES 
  set resi1 [getNextResID $resi ]
  set atomidList [list]
  lappend atomidList [getAtomID $molid $resi "CA"]
  lappend atomidList [getAtomID $molid $resi1 "N"]
  lappend atomidList [getAtomID $molid $resi "C"]
  lappend atomidList [getAtomID $molid $resi "O"]
  return $atomidList
}




#############################################
set GRO "02_prot.gro"
set molid [ mol new $GRO type gro ]

set ofp [open "Impropers.dat" w]
set Omega2  [getOmega2Atoms $molid $NRES ]
lassign $Omega2 a1 a2 a3 a4
puts $ofp "$a1 $a2 $a3 $a4"
set Omega1  [getOmega1Atoms $molid $IRES ]
lassign $Omega1 a1 a2 a3 a4
puts $ofp "$a1 $a2 $a3 $a4"

close $ofp
exit 0
